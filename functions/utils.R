# Storr utils
extract_storr <- function(
    storr_path = "mpaeu_storr",
    what = "succeeded"
) {
    require(storr)
    st <- storr_rds(storr_path)
    keys <- st$list()
    status <- st$mget(keys)
    status <- unlist(lapply(status, \(x) x[[1]]), use.names = F)
    keys <- keys[status == what]
    cli::cli_alert_success("{.val {length(keys)}} key{?s} found for {.var {what}}")
    return(keys)
}

clean_storr <- function(
    storr_path = "mpaeu_storr",
    what = "running"
) {
    require(storr)
    st <- storr_rds(storr_path)
    keys <- st$list()
    status <- st$mget(keys)
    status <- unlist(lapply(status, \(x) x[[1]]), use.names = F)
    keys <- keys[status == what]
    cli::cli_alert_danger("Are you sure you want to remove {.val {length(keys)}} key{?s} for {.var {what}}?")
    input <- readline("Type Y/n\n")
    if (input == "Y") {
        st$del(keys)
        cli::cli_alert_success("{.val {length(keys)}} key{?s} for {.var {what}} removed.")
    } else {
        cli::cli_alert_success("No keys removed")
    }
    return(invisible(NULL))
}

extract_errors_storr <- function(
    storr_path = "mpaeu_storr",
    full_error = FALSE
) {
    require(storr)
    st <- storr_rds(storr_path)
    keys <- st$list()
    status <- st$mget(keys)
    status <- unlist(lapply(status, \(x) x[[1]]), use.names = F)
    keys <- keys[status == "failed"]
    if (length(keys) > 0) {
        status <- st$mget(keys)
        if (full_error) {
            errors <- lapply(status, \(x){
                x[[2]]
            })
        } else {
            errors <- lapply(status, \(x){
                x[[2]][[1]]
            })
            errors <- unlist(errors, use.names = F)
        }
        names(errors) <- keys
        return(errors)
    } else {
        cli::cli_alert_success("No errors found")
        return(invisible(NULL))
    }
}
