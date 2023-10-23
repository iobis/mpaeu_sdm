####

library(glue)
library(terra)
library(dismo)
library(tidyverse)
library(lme4)
library(parallel)

# Set species
sp <- c(101, 103)

plot_theme <- theme_light() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = 8),
        legend.key.width = unit(4, "pt"),
        legend.key.height = unit(4, "pt"))

outfolder <- "results/vsp_testing_old"

# Function
extract_metrics <- function(file) {
  metrics <- read.csv(file)
  metrics <- metrics %>%
    pivot_longer(1:ncol(metrics), values_to = "value", names_to = "metric") %>%
    mutate(model = gsub(
      "_spatial_gridmetrics.csv|_spatial_latmetrics.csv|_randommetrics.csv",
      "",
      gsub(
        "bias_high_rep[[:digit:]]*_|bias_low_rep[[:digit:]]*_|nobias_high_rep[[:digit:]]*_|nobias_low_rep[[:digit:]]*_", "",
        gsub(
          glue("{outfolder}/[[:digit:]]*[[:digit:]]/sp[[:digit:]]*[[:digit:]]_"),
          "",
          file
        )))) %>%
    mutate(replicate = str_extract(file, "rep[[:digit:]]")) %>%
    mutate(bias = ifelse(grepl("nobias", file), "no", "yes")) %>%
    mutate(nsamp = ifelse(grepl("low", file), "low", "high"))
    
  return(metrics)
}


for (i in 1:length(sp)) {
  for (ty in c("spatial_gridmetrics", "spatial_latmetrics", "randommetrics")) {
    m_files <- list.files(glue("{outfolder}/{sp[i]}/"),
                             pattern = ty)
    m_files <- glue("{outfolder}/{sp[i]}/{m_files}")
    eval_metrics <- bind_rows(lapply(m_files, extract_metrics))
    eval_metrics$cv_type <- ty
    if (ty == "spatial_gridmetrics") {
      all_types <- eval_metrics
    } else {
      all_types <- bind_rows(all_types, eval_metrics)
    }
  }
  all_types$species <- sp[i]
  if (i == 1) {
    all_metrics <- all_types
  } else {
    all_metrics <- bind_rows(all_metrics, all_types)
  }
}


# Get only summaries
summ_metrics <- all_metrics %>%
  group_by(species, cv_type, replicate, model, metric, bias, nsamp) %>%
  summarise(value = mean(value),
            n_models = n())


### AUC metrics ----
auc_metrics <- all_metrics %>%
  filter(metric == "auc",
         cv_type == "spatial_gridmetrics") %>%
  mutate(model = toupper(gsub("_", " ", model)),
         bias = ifelse(bias == "yes", "With bias", "No bias"),
         nsamp = ifelse(nsamp == "high", "High sample", "Low sample")) %>%
  mutate(model = as_factor(model),
         cv_type = as_factor(cv_type),
         species = as_factor(species),
         bias = as_factor(bias),
         nsamp = as_factor(nsamp))

ggplot(data = auc_metrics, aes(x = model, y = value, color = model))+
  geom_boxplot(width = .42, outlier.shape = NA) +
  geom_jitter(alpha = .2, width = .12, size = .5) +
  scale_color_discrete("Model") +
  ylab("Mean AUC") + xlab(NULL) +
  facet_wrap(~ bias + nsamp) + plot_theme

ggsave("auc_metrics_vsp.png", width = 1280*2, height = 720*2, unit = "px")

### CBI metrics ----
cbi_metrics <- all_metrics %>%
  filter(metric == "cbi",
         cv_type == "spatial_gridmetrics") %>%
  mutate(model = toupper(gsub("_", " ", model)),
         bias = ifelse(bias == "yes", "With bias", "No bias"),
         nsamp = ifelse(nsamp == "high", "High sample", "Low sample")) %>%
  mutate(model = as_factor(model),
         cv_type = as_factor(cv_type),
         species = as_factor(species),
         bias = as_factor(bias),
         nsamp = as_factor(nsamp))

ggplot(data = cbi_metrics, aes(x = model, y = value, color = model))+
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_boxplot(width = .42, outlier.shape = NA) +
  geom_jitter(alpha = .3, width = .12, size = .5) +
  scale_color_discrete("Model") +
  ylab("Mean CBI") + xlab(NULL) +
  facet_wrap(~ bias + nsamp) + plot_theme

ggsave("cbi_metrics_vsp.png", width = 1280*2, height = 720*2, unit = "px")


### TSS (p10) metrics ----
tss_metrics <- all_metrics %>%
  filter(metric == "tss_p10",
         cv_type == "spatial_gridmetrics") %>%
  mutate(model = toupper(gsub("_", " ", model)),
         bias = ifelse(bias == "yes", "With bias", "No bias"),
         nsamp = ifelse(nsamp == "high", "High sample", "Low sample")) %>%
  mutate(model = as_factor(model),
         cv_type = as_factor(cv_type),
         species = as_factor(species),
         bias = as_factor(bias),
         nsamp = as_factor(nsamp))

ggplot(data = tss_metrics, aes(x = model, y = value, color = model))+
  geom_boxplot(width = .42, outlier.shape = NA) +
  geom_jitter(alpha = .3, width = .12, size = .5) +
  scale_color_discrete("Model") +
  ylab("Mean TSS (p10 threshold)") + xlab(NULL) +
  facet_wrap(~ bias + nsamp) + plot_theme

ggsave("tss_metrics_vsp.png", width = 1280*2, height = 720*2, unit = "px")

auc_bias <- auc_metrics %>% filter(bias == "With bias")
auc_nobias <- auc_metrics %>% filter(bias == "No bias")


m_auc_bias <- lmer(value ~ model + (1 | nsamp) + (1 | species),
                   data = auc_bias)

(re.effects <- sjPlot::plot_model(m_auc_bias, type = "re", terms = c("nsamp", "species"), show.values = TRUE))


#### Compare predictions

# Function
extract_map_metric <- function(index, files) {
  
  file <- files[index]
  
  results <- data.frame(
    i_metric = NA,
    jacc_mtp = NA,
    jacc_p10 = NA
  )
  
  bias <- ifelse(grepl("nobias", file), "", "bias_")
  nsamp <- ifelse(grepl("high", file), "high", "low")
  
  sp <- gsub("sp", "", str_extract(file, "sp[:digit:]*[:digit:]"))
  replicate <- str_extract(file, "rep[:digit:]*[:digit:]")
  
  base <- glue("{outfolder}/{sp}")
  pred <- rast(glue("{base}/{file}"))
  
  pred <- crop(pred, current)
  pred <- mask(pred, current)
  
  pred <- (pred - global(pred, min, na.rm = T)[,1])/(global(pred, max, na.rm = T)[,1] - global(pred, min, na.rm = T)[,1])
  
  jacc <- function(rast1, rast2) {
    comb <- rast1 + rast2
    inter <- comb
    union <- comb
    inter[] <- ifelse(union[] >= 2, 1, 0)
    union[] <- ifelse(comb[] >= 1, 1, 0)
    
    cinter <- freq(inter)
    cunion <- freq(union)
    
    return((cinter$count[cinter$value == 1] / cunion$count[cunion$value == 1]))
  }
  
  if (grepl("current", file)) {
    
    points <- arrow::read_parquet(glue("data/virtual_species/key={sp}/occurrences_{nsamp}_{bias}{replicate}.parquet"))
    points$pred <- terra::extract(pred, points[,1:2], ID = F)[,1]
    
    mtp <- modEvA::getThreshold(obs = points$presence, pred = points$pred, threshMethod = "MTP")
    p10 <- modEvA::getThreshold(obs = points$presence, pred = points$pred, threshMethod = "MTP", quant = 0.1)
    
    binary_pred_mtp <- terra::app(pred, function(x) ifelse(x >= mtp, 1, 0))
    binary_pred_p10 <- terra::app(pred, function(x) ifelse(x >= p10, 1, 0))
    
    results$i_metric <- dismo::nicheOverlap(raster(pred), raster(current),
                                             mask = F, checkNegatives = F)
    results$jacc_mtp <- jacc(current_binary, binary_pred_mtp)
    results$jacc_p10 <- jacc(current_binary, binary_pred_p10)
    
    results$type <- "current"
  } else {
    if (grepl("ssp1", file)) {
      results$i_metric <- dismo::nicheOverlap(raster(pred), raster(ssp1),
                                              mask = F, checkNegatives = F)
      results$type <- "ssp1"
    } else {
      results$i_metric <- dismo::nicheOverlap(raster(pred), raster(ssp5),
                                              mask = F, checkNegatives = F)
      results$type <- "ssp5"
    }
  }
  
  results$bias <- ifelse(grepl("nobias", file), "no", "yes")
  results$nsamp <- nsamp
  results$model <- gsub(
    "_current.tif|_ssp1.tif|_ssp5.tif",
    "",
    gsub(
      "bias_high_rep[[:digit:]]*_|bias_low_rep[[:digit:]]*_|nobias_high_rep[[:digit:]]*_|nobias_low_rep[[:digit:]]*_", "",
      gsub(
        glue("sp[[:digit:]]*[[:digit:]]_"),
        "",
        file
      )))
  results$sp <- sp
  results$rep <- replicate
  
  return(results)
}

for (i in 1:length(sp)) {
  cat("Running", sp[i], "\n")
  
  # Get original suitability
  current <- rast(glue("data/virtual_species/key={sp[i]}/suitability_key{sp[i]}.tif"))
  current_binary <- rast(glue("data/virtual_species/key={sp[i]}/pa_key{sp[i]}.tif"))
  ssp1 <- rast(glue("data/virtual_species/key={sp[i]}/ssp1_suitability_key{sp[i]}.tif"))
  ssp5 <- rast(glue("data/virtual_species/key={sp[i]}/ssp5_suitability_key{sp[i]}.tif"))
  
  current <- crop(current, ext(c(-41, 47, 20, 89)))
  current_binary <- crop(current_binary, ext(c(-41, 47, 20, 89)))
  ssp1 <- crop(ssp1, ext(c(-41, 47, 20, 89)))
  ssp5 <- crop(ssp5, ext(c(-41, 47, 20, 89)))
  
  m_files <- list.files(glue("{outfolder}/{sp[i]}"),
                        pattern = ".tif")
  m_files <- m_files[-grep(".json", m_files)]
  
  # all_types <- lapply(cli::cli_progress_along(m_files), extract_map_metric,
  #                     files = m_files)
  all_types <- mclapply(1:length(m_files), extract_map_metric, files = m_files,
                        mc.cores = 6)
  
  all_types <- bind_rows(all_types)

  if (i == 1) {
    all_metrics <- all_types
  } else {
    all_metrics <- bind_rows(all_metrics, all_types)
  }
}

all_metrics

teste <- all_metrics %>%
  filter(type == "current") 

ggplot(teste) +
  geom_boxplot(aes(x = model, y = i_metric)) +
  facet_wrap(~ bias + nsamp) + plot_theme

ggplot(teste) +
  geom_boxplot(aes(x = model, y = jacc_mtp)) +
  facet_wrap(~ bias + nsamp) + plot_theme

teste <- all_metrics %>%
  filter(type == "ssp5") 

ggplot(teste) +
  geom_boxplot(aes(x = model, y = i_metric)) +
  facet_wrap(~ bias + nsamp) + plot_theme




for (i in 1:length(sp)) {
  time_f <- list.files(glue("results/vsp_testing_old/{sp[i]}"), pattern = "timings",
                       full.names = T)
  
  sp_timings <- lapply(time_f, function(x){
    tim <- read.csv(x)
    colnames(tim) <- c("step", "time_minutes")
    
    bias <- ifelse(grepl("nobias", x), "no", "yes")
    nsamp <- ifelse(grepl("high", x), "high", "low")
    
    sp <- gsub("sp", "", str_extract(x, "sp[:digit:]*[:digit:]"))
    replicate <- str_extract(x, "rep[:digit:]*[:digit:]")
    
    tim$model <- gsub(
      "_timings.csv",
      "",
      gsub(
        "bias_high_rep[[:digit:]]*_|bias_low_rep[[:digit:]]*_|nobias_high_rep[[:digit:]]*_|nobias_low_rep[[:digit:]]*_", "",
        gsub(
          glue("sp[[:digit:]]*[[:digit:]]_"),
          "",
          gsub("results/vsp_testing_old/[[:digit:]]*[[:digit:]]/", "", x)
        )))
    tim$bias <- bias
    tim$nsamp <- nsamp
    tim$sp <- sp
    tim$rep <- replicate
    
    return(tim)
  })
  
  sp_timings <- bind_rows(sp_timings)
  
  if (i == 1) {
    all_timings <- sp_timings
  } else {
    all_timings <- bind_rows(all_timings, sp_timings)
  }
}

timings_summ <- all_timings %>%
  group_by(step, model) %>%
  summarise(mean_t = mean(time_minutes),
          uppr_t = mean(time_minutes) + sd(time_minutes),
          lwr_t = mean(time_minutes) - sd(time_minutes)) %>%
  filter(step != "response curve") %>%
  mutate(model = toupper(gsub("_", " ", model)))


ggplot(timings_summ, aes(color = model)) +
  geom_line(aes(y = mean_t, x = step, group = model)) +
  #geom_ribbon(aes(ymin = lwr_t, ymax = uppr_t, x = step, group = model))
  xlab(NULL) + ylab("Cummulative mean time (in minutes, from the start)") +
  scale_color_discrete("Model") +
  scale_x_discrete(limits = c("tuning", "cv", "evaluate dataset", "evaluate final"),
                   labels = c("Tuning", "Cross-validation", "Evaluate dataset", "Evaluate final model")) +
  plot_theme + theme(axis.text.x = element_text(angle = 0))

ggsave("timings_vsp.png", width = 1280*2, height = 720*2, unit = "px")




### Plot example maps for reference

ori_curr <- rast("data/virtual_species/key=101/suitability_key101.tif")
ori_ssp1 <- rast("data/virtual_species/key=101/ssp1_suitability_key101.tif")

max_curr <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_maxnet_current.tif")
max_ssp1 <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_maxnet_ssp1.tif")

brt_curr <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_brt_naive_current.tif")
brt_ssp1 <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_brt_naive_ssp1.tif")

las_curr <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_lasso_current.tif")
las_ssp1 <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_lasso_ssp1.tif")

rfd_curr <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_rf_downsampled_current.tif")
rfd_ssp1 <- rast("results/vsp_testing_old/101/sp101_bias_high_rep1_rf_downsampled_ssp1.tif")

europe <- vect("data/shapefiles/mpa_europe_starea_v2.shp")

ori_curr <- (ori_curr - global(ori_curr, min, na.rm = T)[,1])/(global(ori_curr, max, na.rm = T)[,1] - global(ori_curr, min, na.rm = T)[,1])
ori_ssp1 <- (ori_ssp1 - global(ori_ssp1, min, na.rm = T)[,1])/(global(ori_ssp1, max, na.rm = T)[,1] - global(ori_ssp1, min, na.rm = T)[,1])
max_curr <- (max_curr - global(max_curr, min, na.rm = T)[,1])/(global(max_curr, max, na.rm = T)[,1] - global(max_curr, min, na.rm = T)[,1])
max_ssp1 <- (max_ssp1 - global(max_ssp1, min, na.rm = T)[,1])/(global(max_ssp1, max, na.rm = T)[,1] - global(max_ssp1, min, na.rm = T)[,1])
brt_curr <- (brt_curr - global(brt_curr, min, na.rm = T)[,1])/(global(brt_curr, max, na.rm = T)[,1] - global(brt_curr, min, na.rm = T)[,1])
brt_ssp1 <- (brt_ssp1 - global(brt_ssp1, min, na.rm = T)[,1])/(global(brt_ssp1, max, na.rm = T)[,1] - global(brt_ssp1, min, na.rm = T)[,1])
las_curr <- (las_curr - global(las_curr, min, na.rm = T)[,1])/(global(las_curr, max, na.rm = T)[,1] - global(las_curr, min, na.rm = T)[,1])
las_ssp1 <- (las_ssp1 - global(las_ssp1, min, na.rm = T)[,1])/(global(las_ssp1, max, na.rm = T)[,1] - global(las_ssp1, min, na.rm = T)[,1])
rfd_curr <- (rfd_curr - global(rfd_curr, min, na.rm = T)[,1])/(global(rfd_curr, max, na.rm = T)[,1] - global(rfd_curr, min, na.rm = T)[,1])
rfd_ssp1 <- (rfd_ssp1 - global(rfd_ssp1, min, na.rm = T)[,1])/(global(rfd_ssp1, max, na.rm = T)[,1] - global(rfd_ssp1, min, na.rm = T)[,1])

ori_curr <- crop(ori_curr, europe)
ori_ssp1 <- crop(ori_ssp1, europe)
max_curr <- crop(max_curr, europe)
max_ssp1 <- crop(max_ssp1, europe)
brt_curr <- crop(brt_curr, europe)
brt_ssp1 <- crop(brt_ssp1, europe)
las_curr <- crop(las_curr, europe)
las_ssp1 <- crop(las_ssp1, europe)
rfd_curr <- crop(rfd_curr, europe)
rfd_ssp1 <- crop(rfd_ssp1, europe)

max_curr <- mask(max_curr, ori_curr)
max_ssp1 <- mask(max_ssp1, ori_curr)
brt_curr <- mask(brt_curr, ori_curr)
brt_ssp1 <- mask(brt_ssp1, ori_curr)
las_curr <- mask(las_curr, ori_curr)
las_ssp1 <- mask(las_ssp1, ori_curr)
rfd_curr <- mask(rfd_curr, ori_curr)
rfd_ssp1 <- mask(rfd_ssp1, ori_curr)


ori_curr <- as.data.frame(ori_curr, xy = T)
ori_ssp1 <- as.data.frame(ori_ssp1, xy = T)
max_curr <- as.data.frame(max_curr, xy = T)
max_ssp1 <- as.data.frame(max_ssp1, xy = T)
brt_curr <- as.data.frame(brt_curr, xy = T)
brt_ssp1 <- as.data.frame(brt_ssp1, xy = T)
las_curr <- as.data.frame(las_curr, xy = T)
las_ssp1 <- as.data.frame(las_ssp1, xy = T)
rfd_curr <- as.data.frame(rfd_curr, xy = T)
rfd_ssp1 <- as.data.frame(rfd_ssp1, xy = T)

ori_curr$model <- "ori"
ori_ssp1$model <- "ori"
max_curr$model <- "max"
max_ssp1$model <- "max"
brt_curr$model <- "brt"
brt_ssp1$model <- "brt"
las_curr$model <- "las"
las_ssp1$model <- "las"
rfd_curr$model <- "rfd"
rfd_ssp1$model <- "rfd"

colnames(ori_curr)[3] <- "value"
colnames(ori_ssp1)[3] <- "value"
colnames(max_curr)[3] <- "value"
colnames(max_ssp1)[3] <- "value"
colnames(brt_curr)[3] <- "value"
colnames(brt_ssp1)[3] <- "value"
colnames(las_curr)[3] <- "value"
colnames(las_ssp1)[3] <- "value"
colnames(rfd_curr)[3] <- "value"
colnames(rfd_ssp1)[3] <- "value"

current <- data.frame(rbind(max_curr, brt_curr, las_curr, rfd_curr))
future <- data.frame(rbind(max_curr, brt_curr, las_curr, rfd_curr))


base <- rnaturalearth::ne_countries(returnclass = "sf")

# New facet label names for dose variable
new_labs <- c("Maxent", "BRT", "LASSO", "RF down-sampled")
names(new_labs) <- c("max", "brt", "las", "rfd")

p1 <- ggplot()+
  geom_sf(data = base, fill = "grey80", color = "grey80") +
  geom_raster(data = current, aes(x = x, y = y, fill = value)) +
  coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
  scale_fill_distiller("", direction = 1) +
  xlab(NULL) + ylab(NULL) +
  theme_light() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "grey10"),
        text = element_text(size = 6)) +
  facet_wrap(~ model, labeller = labeller(model = new_labs))

(p2 <- ggplot()+
  geom_sf(data = base, fill = "grey80", color = "grey80") +
  geom_raster(data = ori_curr, aes(x = x, y = y, fill = value)) +
  coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
  scale_fill_distiller("", direction = 1) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Virtual species 101") +
  theme_light() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "grey10"),
        text = element_text(size = 6)))

library(patchwork)

p2 + p1

ggsave("map_current.png", width = 1280*2, height = 720*2, unit = "px")


p1 <- ggplot()+
  geom_sf(data = base, fill = "grey80", color = "grey80") +
  geom_raster(data = future, aes(x = x, y = y, fill = value)) +
  coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
  scale_fill_distiller("", direction = 1) +
  xlab(NULL) + ylab(NULL) +
  theme_light() +
  theme(panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color = "grey10"),
        text = element_text(size = 6)) +
  facet_wrap(~ model, labeller = labeller(model = new_labs))

(p2 <- ggplot()+
    geom_sf(data = base, fill = "grey80", color = "grey80") +
    geom_raster(data = ori_ssp1, aes(x = x, y = y, fill = value)) +
    coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
    scale_fill_distiller("", direction = 1) +
    xlab(NULL) + ylab(NULL) +
    ggtitle("Virtual species 101") +
    theme_light() +
    theme(panel.border = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(color = "grey10"),
          text = element_text(size = 6)))


p2 + p1

ggsave("map_future.png", width = 1280*2, height = 720*2, unit = "px")
