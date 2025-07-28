############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ SDM - bootstrap functions #########################

source("functions/components_model_species.R")

bootstrap_model <- function(species, iterations = 20,
                            model_acro = "mpaeu",
                            results_folder = "results",
                            target = "all",
                            n_back_max = NULL, verbose = FALSE, retry_max = 10) {
    temp_outf <- glue::glue("{results_folder}/taxonid={species}/model={model_acro}/predictions/temp_pred/")
    fs::dir_create(temp_outf)

    log_file <- jsonlite::read_json(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=log.json"))

    best_conf <- log_file[["model_bestparams"]]
    best_variables <- unlist(log_file[["model_details"]][["variables"]])
    good_models <- unlist(log_file[["model_good"]])
    metric_type <- log_file[["model_good_metric"]][[1]]
    group <- log_file[["group"]][[1]]

    if (target == "all") {
        to_do <- good_models
        if (is.numeric(to_do)) {
          to_do <- "esm"
        }
    } else if (target == "best") {
        metrics <- lapply(good_models, function(mn) {
            mfiles <- list.files(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/metrics/"), full.names = T)
            mf <- mfiles[grepl("what=cvmetr", mfiles)]
            mf <- mf[grepl(paste0("method=", substr(mn, 1, 2)), mf)]
            mo <- arrow::read_parquet(mf)
            return(mean(mo[metric_type][[1]], na.rm = T))
        })
        to_do <- good_models[which.max(metrics)]
    } else {
        to_do <- base::intersect(target, good_models)
    }

    if (length(to_do) < 1) {
        cli::cli_abort("No good model corresponding to target='{target}'. Available: {.val {good_models}} model{?s}")
    }

    # For back compatibility
    if ("hab_depth" %in% names(log_file)) {
        hab_depth <- log_file[["hab_depth"]][[1]]
    } else {
        eco_info <- arrow::open_csv_dataset("data/species_ecoinfo.csv") %>%
            filter(taxonID == species) %>%
            collect()

        if (nrow(eco_info) < 1) {
            eco_info <- obissdm::mp_get_ecoinfo(
                species_list = species,
                outfile = NULL,
                return_table = TRUE,
                show_progress = FALSE
            )
        }

        # Select ecological information for the species
        hab <- eco_info$mode_life
        hab_depth <- hab_to_depth(hab)
    }

    # Load env layers
    env <- obissdm::get_envofgroup(unlist(group),
        depth = hab_depth, load_all = T,
        conf_file = "sdm_conf.yml", verbose = verbose
    )
    env <- terra::subset(env$layers, best_variables)

    mask_l <- terra::rast(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/predictions/taxonid={species}_model={model_acro}_mask_cog.tif"))
    mask_l <- terra::subset(mask_l, "fit_region")
    mask_l[mask_l == 0] <- NA

    env_masked <- terra::mask(env, mask_l)

    # Load points
    pts <- arrow::read_parquet(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=fitocc.parquet"))

    quad_n <- unlist(log_file[["model_details"]][["background_size"]][[1]])

    if (!is.null(n_back_max)) {
        quad_n <- ifelse(quad_n > n_back_max, n_back_max, quad_n)
    }

    env_xy <- terra::as.data.frame(prod(env_masked), xy = T)[, 1:2]
    env_xy_samp <- sample(seq_len(nrow(env_xy)), size = quad_n)
    env_xy <- env_xy[env_xy_samp, ]
    colnames(env_xy) <- c("decimalLongitude", "decimalLatitude")

    sampled_quad <- cbind(
        env_xy,
        presence = 0,
        terra::extract(env_masked, env_xy, ID = F)
    )

    for (i in seq_len(iterations)) {
        if (verbose) cli::cat_line(cli::bg_cyan(glue::glue("Running iteration {i} out of {iterations}")))
        samp_pts <- sample(seq_len(nrow(pts)), nrow(pts), replace = T)
        samp_quad <- sample(seq_len(nrow(sampled_quad)), nrow(sampled_quad), replace = T)

        sp_data <- obissdm::mp_prepare_data(pts[samp_pts, ],
            eval_data = NULL,
            species_id = unlist(log_file[["scientificName"]]),
            env_layers = env_masked,
            pred_quad = sampled_quad[samp_quad, ],
            verbose = verbose
        )

        for (z in seq_len(length(to_do))) {
            if (to_do[z] == "maxent") {
                bp <- best_conf[["maxent"]][[1]]
            } else if (to_do[z] == "esm") {
                bp <- best_conf[[to_do[z]]]
                bp$individual_parameters <- as.data.frame(dplyr::bind_rows(bp$individual_parameters))
                # Deal with exception, when the resulting object contain lists
                if (class(bp$individual_parameters[,1]) == "list") {
                  new_df <- data.frame(
                    remult = unlist(bp$individual_parameters$remult),
                    features = unlist(bp$individual_parameters$features),
                    part = unlist(bp$individual_parameters$part)
                  )
                  bp$individual_parameters <- new_df
                }
                bp$variable_combinations <- lapply(bp$variable_combinations, function(x){
                  unlist(x)
                })
                bp$scores <- unlist(bp$scores)
            } else {
                bp <- best_conf[[to_do[z]]]
            }

            fit <- try(obissdm::model_bootstrap(sp_data, algo = to_do[z], params = bp, verbose = verbose), silent = T)
            
            if (inherits(fit, "try-error") && !is.null(retry_max)) {
              initc <- 0
              while (initc <= retry_max & inherits(fit, "try-error")) {
                samp_pts <- sample(seq_len(nrow(pts)), nrow(pts), replace = T)
                samp_quad <- sample(seq_len(nrow(sampled_quad)), nrow(sampled_quad), replace = T)
                
                sp_data <- obissdm::mp_prepare_data(pts[samp_pts, ],
                                                    eval_data = NULL,
                                                    species_id = unlist(log_file[["scientificName"]]),
                                                    env_layers = env_masked,
                                                    pred_quad = sampled_quad[samp_quad, ],
                                                    verbose = verbose
                )
                fit <- try(obissdm::model_bootstrap(sp_data, algo = to_do[z], params = bp, verbose = verbose), silent = T)
                initc <- initc + 1
              }
              if (inherits(fit, "try-error")) stop("Error when trying refit of ", to_do[z])
            }

            pred <- predict_bootstrap(fit, sp_data, species,
                                      group = group, hab_depth = hab_depth,
                                      pred_out = temp_outf, outacro = model_acro,
                                      iteration = i,
                                      type = ifelse(to_do[z] == "esm", "esm", "normal"),
                                      aggregate = 4,
                                      verbose = verbose)
        }
    }

    return(invisible(NULL))
}

predict_bootstrap <- function(fit, sdm_data, species, group, hab_depth,
                              pred_out, outacro, iteration, type = "normal",
                              aggregate, verbose = FALSE) {

    scenarios <- data.frame(
        scenario = c("current", rep(c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
            each = 2
        )),
        year = c(NA, rep(c("dec50", "dec100"), 5))
    )

    if (type == "esm") {
        model_name <- "esm"
    } else {
        model_name <- fit$name
    }

    for (k in seq_len(nrow(scenarios))) {
        if (is.na(scenarios$year[k])) {
            period <- NULL
        } else {
            period <- scenarios$year[k]
        }

        if (verbose) cli::cli_alert_info("Predicting scenario {k} of {nrow(scenarios)}.")
        outpred <- paste0(
            pred_out, "taxonid=", species, "_model=", outacro,
            "_method=", model_name, "_scen=", scenarios$scenario[k],
            ifelse(is.null(period), "", paste0("_", period)), "_it=", iteration, ".tif"
        )

        env_to_pred <- obissdm::get_envofgroup(unlist(group),
            depth = unlist(hab_depth), load_all = T,
            scenario = scenarios$scenario[k],
            period = period,
            env_folder = "data/env",
            conf_file = "sdm_conf.yml",
            verbose = verbose
        )
        env_to_pred <- terra::subset(env_to_pred$layers, colnames(sdm_data$training)[-1])
        if (aggregate > 0) {
            env_to_pred <- terra::aggregate(env_to_pred, fact = aggregate, na.rm = T)
        }

        if (type == "normal") {
            pred <- predict(fit, env_to_pred)
        } else {
            pred_multi <- lapply(fit, function(x) {
                predict(x, env_to_pred)
            })
            pred_multi <- terra::rast(pred_multi)
            scores_ens <- attr(fit, "scores")
            pred <- weighted.mean(pred_multi, scores_ens, na.rm = T)
        }

        names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))

        pred <- pred * 100
        pred <- terra::as.int(pred)
        terra::writeRaster(pred, outpred, overwrite = T, datatype = "INT1U")
    }

    return(invisible(NULL))
}

bootstrap_sp <- function(species, target = "all",
                         acro = "mpaeu", results_folder = "results",
                         iterations = 20, n_back_max = 50000, retry_max = 10) {

    bt_result <- try(bootstrap_model(
        species,
        iterations = iterations,
        model_acro = acro,
        results_folder = results_folder,
        target = target,
        n_back_max = n_back_max,
        retry_max = retry_max
    ), silent = F)

    if (inherits(bt_result, "try-error")) return(list("failed", bt_result))

    all_preds <- list.files(glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions/temp_pred/"),
        full.names = T
    )

    un_models <- unique(gsub(".*method=([a-zA-Z0-9_]+)_scen.*", "\\1", all_preds))

    for (un in un_models) {
        m_pred <- all_preds[grepl(un, all_preds)]

        scenarios <- gsub(".*scen=([a-zA-Z0-9_]+)_it=.*", "\\1", m_pred)
        scenarios <- unique(scenarios)

        for (sc in scenarios) {
            f_preds <- list.files(glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions"),
                full.names = T
            )
            f_preds <- f_preds[grepl(paste0("method=", un, "_scen=", sc), f_preds)]
            f_preds <- f_preds[!grepl("bootcv", f_preds)]

            base_layer <- terra::rast(f_preds)

            m_pred_sc <- m_pred[grepl(sc, m_pred)]

            pred_rast <- terra::rast(m_pred_sc)

            #cv_rast <- (terra::app(pred_rast, "sd") / mean(pred_rast))
            sd_rast <- terra::app(pred_rast, "sd")
            mean_rast <- mean(pred_rast)
            cv_rast <- c(sd_rast, mean_rast)

            cv_rast <- terra::disagg(cv_rast, fact = 4)
            
            # To handle special case of older wavefetch (back compatibility)
            if (terra::ext(cv_rast) != terra::ext(base_layer)) {
              cv_rast <- terra::crop(cv_rast, base_layer)
            }
            
            cv_rast <- terra::mask(cv_rast, base_layer)

            #cv_rast <- cv_rast * 100
            cv_rast <- terra::as.int(cv_rast)

            outf <- basename(m_pred_sc[1])
            outf <- gsub("_it=*.*", "", outf)
            outf <- glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions/{outf}_what=bootcv.tif")

            terra::writeRaster(cv_rast, outf, overwrite = T, datatype = "INT2U")
            cogeo_optim(outf)
        }
    }

    if (length(un_models) > 1) {
        for (sc in scenarios) {
            f_preds <- list.files(glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions"),
                full.names = T
            )
            f_preds <- f_preds[grepl(paste0("_scen=", sc), f_preds)]
            f_preds <- f_preds[!grepl("bootcv", f_preds)]
            f_preds <- f_preds[grepl(paste0(un_models, collapse = "|"), f_preds)]

            base_layer <- terra::rast(f_preds)
            base_layer <- mean(base_layer)

            m_pred <- all_preds[grepl(paste0("_scen=", sc), all_preds)]
            m_pred <- m_pred[grepl(paste0(un_models, collapse = "|"), m_pred)]

            pred_rast <- lapply(seq_len(iterations), function(x) NULL)

            for (it in seq_len(iterations)) {
                pred_rast[[it]] <- mean(terra::rast(m_pred[grepl(paste0("it=", it, ".tif"), m_pred)]))
            }

            pred_rast <- terra::rast(pred_rast)

            #cv_rast <- (terra::app(pred_rast, "sd") / mean(pred_rast))
            sd_rast <- terra::app(pred_rast, "sd")
            mean_rast <- mean(pred_rast)
            cv_rast <- c(sd_rast, mean_rast)

            cv_rast <- terra::disagg(cv_rast, fact = 4)
            
            # To handle special case of older wavefetch (back compatibility)
            if (terra::ext(cv_rast) != terra::ext(base_layer)) {
              cv_rast <- terra::crop(cv_rast, base_layer)
            }
            
            cv_rast <- terra::mask(cv_rast, base_layer)

            #cv_rast <- cv_rast * 100
            cv_rast <- terra::as.int(cv_rast)

            outf <- glue::glue("taxonid={species}_model={acro}_method=ensemble_scen={sc}")
            outf <- glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions/{outf}_what=bootcv.tif")

            terra::writeRaster(cv_rast, outf, overwrite = T, datatype = "INT2U")
            cogeo_optim(outf)
        }
    }

    jsonf <- glue::glue("{results_folder}/taxonid={species}/model={acro}/taxonid={species}_model={acro}_what=log.json")
    log_file <- jsonlite::read_json(jsonf)
    log_file$model_uncertainty$bootstrap_status <- "done"
    log_file$model_uncertainty$bootstrap_iterations <- iterations
    if (length(un_models) > 1) {
        log_file$model_uncertainty$bootstrap_models <- c(un_models, "ensemble")
    } else {
        log_file$model_uncertainty$bootstrap_models <- un_models
    }
    log_file$model_uncertainty$bootstrap_max_n <- n_back_max
    jsonlite::write_json(log_file, path = jsonf, pretty = TRUE)

    fs::dir_delete(glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions/temp_pred/"))

    return("done")
}

# A problem was identified that, even if the bootstrap was done for all models
# it would only save one of those in the log file. This does not affect any results
# but we fixed the json files using the following function. Any new bootstrap realization
# is already fixed by defaul.
fix_boot_json <- function(species, results_folder, acro) {
    preds <- list.files(glue::glue("{results_folder}/taxonid={species}/model={acro}/predictions"))
    bootf <- preds[grepl("what=bootcv", preds)]
    bootf <- gsub("_scen*.*", "", bootf)
    bootf <- gsub("*.*_method=", "", bootf)
    spd <- "not_changed"
    if (length(unique(bootf)) > 1) {
        counts <- table(bootf)
        if (length(unique(counts)) == 1) {
            # Needs correction
            models_done <- unique(bootf)
            jsonf <- glue::glue("{results_folder}/taxonid={species}/model={acro}/taxonid={species}_model={acro}_what=log.json")
            log_file <- jsonlite::read_json(jsonf)
            log_file$model_uncertainty$bootstrap_models <- models_done
            jsonlite::write_json(log_file, path = jsonf, pretty = TRUE)
            spd <- "changed"
        }
    }
    return(spd)
}
# END
