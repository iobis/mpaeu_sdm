monitor_model_shiny <- function(model_acro, total_sp, update_time = 60) {
  
  require(shiny)
  require(gridlayout)
  require(bslib)
  require(ggplot2)
  
  runApp(list(
    ui = grid_page(
      layout = c(
        "header  header   area4   ",
        "sidebar bluePlot bluePlot",
        "area3   area3    area3   "
      ),
      row_sizes = c(
        "125px",
        "1.40fr",
        "0.6fr"
      ),
      col_sizes = c(
        "310px",
        "1.12fr",
        "0.88fr"
      ),
      gap_size = "1rem",
      grid_card(
        area = "sidebar",
        card_header("Models"),
        card_body(
          "Number of species:",
          textOutput(outputId = "numberOutput"),
          "In execution:",
          textOutput(outputId = "executingOutput"),
          "Failed/Skipped:",
          textOutput(outputId = "failedOutput"),
          "Done:",
          textOutput(outputId = "doneOutput")
        )
      ),
      grid_card_text(
        area = "header",
        content = "MPA Europe - SDM modeling tracker",
        alignment = "start",
        is_title = FALSE
      ),
      grid_card_plot(area = "bluePlot"),
      grid_card(
        area = "area3",
        full_screen = TRUE,
        card_header("Errors"),
        card_body(verbatimTextOutput("errorsOutput", placeholder = T), padding = 0)
      ),
      grid_card(
        area = "area4",
        card_body(
          value_box(
            title = "Time running",
            value = textOutput("timeOutput"),
            showcase = bsicons::bs_icon("clock-fill")
          )
        )
      )
    ),
    server = function(input, output, session) {
      
      require(storr)
      outacro <- model_acro
      st <- storr_rds(paste0(outacro, "_storr"))
      
      total <- total_sp
      
      next_whole <- lubridate::ceiling_date(Sys.time(), "10 seconds")
      go_signal <- reactiveVal(FALSE)
      
      first_check <- observe({
        invalidateLater(10)
        req(next_whole - Sys.time() < 0)
        go_signal(TRUE)
        first_check$destroy()
      })
      
      
      dat <- reactiveValues(data = data.frame(type = c("Total", "Remaining", "Executing", "Failed/Skipped", "Done"),
                                              value = c(0, 0, 0, 0, 0),
                                              class = c(rep("Modeling", 2), rep("Results", 3))),
                            time = 0, time_old = Sys.time(),
                            errors = "Wait...")
      
      observe({
        if(!go_signal()) {return(format(Sys.time(), "%H:%M:%S"))}
        isolate({
          
          all <- st$mget(st$list())
          all_read <- st$list()
          all <- unlist(lapply(all, function(x) x[[1]]))
          
          good <- length(all_read[which(all %in% c("done", "succeeded"))])
          failed <- length(all_read[which(all %in% c("failed","no_good_model"))])
          runing <- length(all_read[which(all == "running")])
          low <- length(all_read[which(all %in% c("low_data"))])
          
          dat$data$value[dat$data$type == "Executing"] <- runing
          dat$data$value[dat$data$type == "Failed/Skipped"] <- failed + low
          dat$data$value[dat$data$type == "Done"] <- good
          
          dat$data$value[dat$data$type == "Total"] <- total
          dat$data$value[dat$data$type == "Remaining"] <- total - (good + failed + low)
          
          dif <- round(difftime(Sys.time(), dat$time_old, units = "mins"), 2)
          
          if (dif < 60) {
            dif <- paste(round(difftime(Sys.time(), dat$time_old, units = "mins"), 1), "mins")
          } else if (dif < 1440) {
            dif <- paste(round(difftime(Sys.time(), dat$time_old, units = "hours"), 1), "hours")
          } else {
            dif <- paste(round(difftime(Sys.time(), dat$time_old, units = "hours"), 1), "days")
          }
          
          dat$time <- dif
          
          errors_n <- all_read[which(all %in% c("failed"))]
          errors <- st$mget(rev(errors_n))
          
          dat$errors <- paste("key =", rev(errors_n), "|", unlist(lapply(errors, function(x) x[[2]][1])))
          
        })
        invalidateLater(update_time*1000, session)
        format(Sys.time(), "%H:%M:%S")
      })
      
      output$bluePlot <- renderPlot({
        ggplot(dat$data) +
          geom_bar(aes(x = type, y = value, fill = type), stat = "identity") +
          scale_fill_manual(values = c("#218380", "#ffbc42", "#8f2d56", "#16697a", "#489fb5")) +
          #facet_wrap(~ class, scales = "free") +
          xlab(NULL) + ylab(NULL) +
          theme_light()
      })
      
      output$numberOutput <- renderText({total})
      output$executingOutput <- renderText({dat$data$value[dat$data$type == "Executing"]})
      output$failedOutput <- renderText({dat$data$value[dat$data$type == "Failed/Skipped"]})
      output$doneOutput <- renderText({dat$data$value[dat$data$type == "Done"]})
      
      output$timeOutput <- renderText({dat$time})
      
      output$errorsOutput <- renderText({dat$errors})
    }
  ))
  
}


#monitor_model_shiny("inteval", 1000)
