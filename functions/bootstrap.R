############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ SDM - bootstrap functions #########################

bootstrap_model <- function(species, iterations = 20,
                            model_acro = "mpaeu",
                            results_folder = "results",
                            target = "all",
                            n_back_max = NULL, verbose = FALSE) {
    temp_outf <- glue::glue("{results_folder}/taxonid={species}/model={model_acro}/predictions/temp_pred/")
    fs::dir_create(temp_outf)

    source("functions/components_model_species.R")
    log_file <- jsonlite::read_json(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=log.json"))

    best_conf <- log_file[["model_bestparams"]]
    best_variables <- unlist(log_file[["model_details"]][["variables"]])
    good_models <- unlist(log_file[["model_good"]])
    metric_type <- log_file[["model_good_metric"]][[1]]
    group <- log_file[["group"]][[1]]

    if (target == "all") {
        to_do <- good_models
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
    env <- obissdm::get_envofgroup(group,
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

    quad_n <- log_file[["model_details"]][["background_size"]][[1]]

    if (!is.null(n_back_max)) {
        quad_n <- ifelse(quad_n > n_back_max, n_back_max, quad_n)
    }

    env_xy <- terra::as.data.frame(prod(env_masked), xy = T)[, 1:2]
    env_xy_samp <- sample(1:nrow(env_xy), size = quad_n)
    env_xy <- env_xy[env_xy_samp, ]
    colnames(env_xy) <- c("decimalLongitude", "decimalLatitude")

    sampled_quad <- cbind(
        env_xy,
        presence = 0,
        terra::extract(env_masked, env_xy, ID = F)
    )

    for (i in 1:iterations) {
        if (verbose) cli::cat_line(cli::bg_cyan(glue::glue("Running iteration {i} out of {iterations}")))
        samp_pts <- sample(1:nrow(pts), nrow(pts), replace = T)
        samp_quad <- sample(1:nrow(sampled_quad), nrow(sampled_quad), replace = T)

        sp_data <- obissdm::mp_prepare_data(pts[samp_pts, ],
            eval_data = NULL,
            species_id = unlist(log_file[["scientificName"]]),
            env_layers = env_masked,
            pred_quad = sampled_quad[samp_quad, ],
            verbose = verbose
        )

        for (z in 1:length(to_do)) {
            if (to_do[z] == "maxent") {
                bp <- best_conf[["maxent"]][[1]]
            } else {
                bp <- best_conf[[to_do[z]]]
            }

            fit <- obissdm::model_bootstrap(sp_data, algo = to_do[z], params = bp, verbose = verbose)

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

    for (k in 1:nrow(scenarios)) {
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

        env_to_pred <- obissdm::get_envofgroup(group,
            depth = hab_depth, load_all = T,
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
            pred_multi <- rast(pred_multi)
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

bootstrap_sp <- function(species, target = "all") {

    bt_result <- try(bootstrap_model(
        species,
        iterations = 20,
        model_acro = acro,
        results_folder = "results",
        target = target,
        n_back_max = 50000
    ), silent = F)

    if (inherits(bt_result, "try-error")) return(list("failed", bt_result))

    all_preds <- list.files(glue::glue("results/taxonid={species}/model={acro}/predictions/temp_pred/"),
        full.names = T
    )

    un_models <- unique(gsub(".*method=([a-zA-Z0-9_]+)_scen.*", "\\1", all_preds))

    for (un in un_models) {
        m_pred <- all_preds[grepl(un, all_preds)]

        scenarios <- gsub(".*scen=([a-zA-Z0-9_]+)_it=.*", "\\1", m_pred)
        scenarios <- unique(scenarios)

        for (sc in scenarios) {
            f_preds <- list.files(glue::glue("results/taxonid={species}/model={acro}/predictions"),
                full.names = T
            )
            f_preds <- f_preds[grepl(paste0("method=", un, "_scen=", sc), f_preds)]
            f_preds <- f_preds[!grepl("bootcv", f_preds)]

            base_layer <- rast(f_preds)

            m_pred_sc <- m_pred[grepl(sc, m_pred)]

            pred_rast <- rast(m_pred_sc)

            cv_rast <- (app(pred_rast, "sd") / mean(pred_rast))

            # To handle cases in that 0/0 == NA
            cv_rast[is.na(cv_rast)] <- 0
            cv_rast <- terra::mask(cv_rast, pred_rast[[1]])

            cv_rast <- disagg(cv_rast, fact = 4)
            cv_rast <- terra::mask(cv_rast, base_layer)

            cv_rast <- cv_rast * 100
            cv_rast <- terra::as.int(cv_rast)

            outf <- basename(m_pred_sc[1])
            outf <- gsub("_it=*.*", "", outf)
            outf <- glue::glue("results/taxonid={species}/model={acro}/predictions/{outf}_what=bootcv.tif")

            writeRaster(cv_rast, outf, overwrite = T, datatype = "INT2U")
            cogeo_optim(outf)
        }
    }

    if (length(un_models) > 1) {
        for (sc in scenarios) {
            f_preds <- list.files(glue::glue("results/taxonid={species}/model={acro}/predictions"),
                full.names = T
            )
            f_preds <- f_preds[grepl(paste0("_scen=", sc), f_preds)]
            f_preds <- f_preds[!grepl("bootcv", f_preds)]
            f_preds <- f_preds[grepl(paste0(un_models, collapse = "|"), f_preds)]

            base_layer <- rast(f_preds)
            base_layer <- mean(base_layer)

            m_pred <- all_preds[grepl(paste0("_scen=", sc), all_preds)]
            m_pred <- m_pred[grepl(paste0(un_models, collapse = "|"), m_pred)]

            pred_rast <- lapply(1:20, function(x) NULL)

            for (it in 1:20) {
                pred_rast[[it]] <- mean(rast(m_pred[grepl(paste0("it=", it, ".tif"), m_pred)]))
            }

            pred_rast <- rast(pred_rast)

            cv_rast <- (app(pred_rast, "sd") / mean(pred_rast))

            # To handle cases in that 0/0 == NA
            cv_rast[is.na(cv_rast)] <- 0
            cv_rast <- terra::mask(cv_rast, pred_rast[[1]])

            cv_rast <- disagg(cv_rast, fact = 4)
            cv_rast <- terra::mask(cv_rast, base_layer)

            cv_rast <- cv_rast * 100
            cv_rast <- terra::as.int(cv_rast)

            outf <- glue::glue("taxonid={species}_model={acro}_method=ensemble_scen={sc}")
            outf <- glue::glue("results/taxonid={species}/model={acro}/predictions/{outf}_what=bootcv.tif")

            writeRaster(cv_rast, outf, overwrite = T, datatype = "INT2U")
            cogeo_optim(outf)
        }
    }

    fs::dir_delete(glue::glue("results/taxonid={species}/model={acro}/predictions/temp_pred/"))

    return("done")
}

# END
