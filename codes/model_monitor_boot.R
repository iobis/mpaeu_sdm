library(storr)
outacro <- "mpaeu"
st <- storr_rds(paste0(outacro, "_boot_storr"))

cli::cli_h1("Monitoring model {cli::bg_cyan(outacro)}")

all <- st$mget(st$list())
all_read <- st$list()
all <- unlist(lapply(all, function(x) x[[1]]))

good <- all_read[which(all %in% c("done", "succeeded"))]
failed <- all_read[which(all %in% c("failed","no_good_model"))]
runing <- all_read[which(all == "running")]
low <- all_read[which(all %in% c("low_data"))]

cli::cli_alert_success("Success: {length(good)}")
cli::cli_alert_info("Runing: {runing} (total of {length(runing)})")
cli::cli_alert_warning("Low number (skipped): {length(low)}")
cli::cli_alert_danger("Failed: {length(failed)}")

i = 1

while (i < 99999999999999) {
  
  Sys.sleep(10)
   
  all <- st$mget(st$list())
  all_read <- st$list()
  all <- unlist(lapply(all, function(x) x[[1]]))
  
  good <- all_read[which(all %in% c("done", "succeeded"))]
  failed <- all_read[which(all %in% c("failed","no_good_model"))]
  runing <- all_read[which(all == "running")]
  low <- all_read[which(all %in% c("low_data"))]
  
  if (interactive()) {
    cat("\014") 
  } else {
    system("clear")
  }
  cli::cli_h1("Monitoring model {cli::bg_cyan(outacro)}")
  cli::cli_alert_success("Success: {length(good)}")
  cli::cli_alert_info("Runing: {runing} (total of {length(runing)})")
  cli::cli_alert_warning("Low number (skipped): {length(low)}")
  cli::cli_alert_danger("Failed: {length(failed)}")
  
  i <- i + 1
}