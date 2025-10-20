#### Write model evaluation
library(dplyr)
results_folder <- "/data/scps/v5/results"
model_acro <- "mpaeu"
done <- as.integer(extract_storr())

# Get evaluation results
results <- read.csv(file.path("internal/sdm_expert_eval", "Evaluation-2025-08-29.csv"))
results_c <- results |>
    filter(!user_code %in% c("s.principe@unesco.org"))

user_info <- readxl::read_xlsx("internal/sdm_expert_eval/original_users.xlsx")
user_info <- user_info |>
    mutate(Name = paste(stringr::str_to_title(Name), stringr::str_to_title(Surname))) |>
    rename(user_code = `E-mail`) |>
    distinct(user_code, .keep_all = T) |>
    select(Name, user_code)

unique_keys <- unique(results_c$species_key)

for (k in done) {
    message("Processing taxonid=", k)
    if (k %in% unique_keys) {
        sp_dat <- results_c |>
            select(-id) |>
            filter(species_key == k) |>
            mutate(question_key = case_when(
                question_key == "Q1" ~ "1. How confident are you in your knowledge of this species’ distribution?",
                question_key == "Q2" ~ "2. Do you think the distribution of species records reflects your knowledge of the species’ actual distribution?",
                question_key == "Q3" ~ "3. To what extent do the spatial predictions align with your understanding of the species’ distribution?",
                question_key == "Q4" ~ "4. How well do the additional model predictions for the current period agree with the main spatial prediction?",
                question_key == "Q5A" ~ "5. Impact and importance for ecosystem processes: Is it an alien or invasive species?",
                question_key == "Q5B" ~ "6. Impact and importance for ecosystem processes: Relevance for the marine environment",
                question_key == "Q6" ~ "7. Notes or comments",
                .default = question_key
            )) |>
            tidyr::pivot_wider(
                names_from = question_key,
                values_from = answer,
                names_sort = TRUE
            ) |>
            left_join(user_info, by = "user_code") |>
            relocate(Name) |>
            select(-user_code) |>
            rename(Evaluator = Name, taxonID = species_key) 

        best_score <- sp_dat$`3. To what extent do the spatial predictions align with your understanding of the species’ distribution?`
        best_score <- best_score[order(best_score)][1]
        best_score_n <- as.integer(strtrim(best_score, 1))

        species_general <- list(
            model = model_acro,
            taxonID = k,
            best_score = best_score,
            best_score_n = best_score_n,
            cbi_score = NA,
            average_score = NA,
            evaluators = unique(sp_dat$Evaluator)
        )
        eval_status <- "evaluated"
    } else {
        eval_status <- "not_evaluated"
        species_general <- NULL
        sp_dat <- NULL
    }
    jsonlite::write_json(
        list(
            status = eval_status,
            summary = species_general,
            evaluations = sp_dat
        ),
        file.path(
            results_folder,
            glue::glue("taxonid={k}/model={model_acro}/taxonid={k}_model={model_acro}_what=experteval.json")
        ),
        pretty = TRUE
    )
}