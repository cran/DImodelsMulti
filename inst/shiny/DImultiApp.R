ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),

  tags$head(
    tags$style(HTML("
    .shiny-file-input-progress {
      display: none !important;
    }
  "))
  ),

  # --- App Title ---
  titlePanel("Repeated measures and multivariate Diversity-Interactions models using DImodelsMulti"),

  # --- Two-Column Layout ---
  fluidRow(
    # Control board for specifying inputs ----------------------------------------------------
    column(
      width = 4,
      div(
        style = "background-color:#f8f9fa; padding:15px; border-radius:10px; height: 95vh; overflow-y: auto; box-shadow: 0 0 5px rgba(0,0,0,0.1);",

        h4("Control Board", class = "text-primary fw-bold mb-3"),

        # Data input
        h5("Data", class = "text-secondary"),
        fileInput("data_file", "Choose data file", accept = ".csv"),
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",", Tab = "\t", Semicolon = ";"),
                     selected = ","),
        helpText(
          tags$em(
            "Note: The loaded data is used only within this session for analysis.",
            "It is not stored, shared, or saved and all data is automatically deleted when the session ends."
          )),
        tags$div(style = "margin-top: 10px;"),
        textOutput("upload_status"),
        tags$div(style = "margin-top: 10px;"),
        actionButton("load_data", "Load Data", icon = icon("file-import")),
        hr(),

        # Model Configuration
        h5("Model configuration", class = "text-secondary"),
        uiOutput("input_selections"),
        hr(),

        # Action Buttons
        # h5("Run Model", class = "text-secondary"),
        actionButton("run_model", "Fit Model(s)", icon = icon("play"), class = "btn-success w-100"),
        br(), br(),
        actionButton("reset_all", "Reset All", icon = icon("redo"), class = "btn-outline-danger w-100")
      )
    ),

    # RIGHT: Outputs ----------------------------------------------------------
    column(
      width = 8,
      div(
        style = "padding:15px;",

        h4("Outputs", class = "text-primary fw-bold mb-3"),

        tabsetPanel(
          id = "output_tabs",
          type = "tabs",

          # Tab 1: Data Preview
          shiny::tabPanel("Data Preview",
                   shinydashboard::box(width = 12, title = " ", status = "primary", solidHeader = TRUE,
                       shiny::tableOutput("data_preview"))
          ),

          # Tab 2: Model Summary
          # shiny::tabPanel("Model Summary",
          #          shinydashboard::box(width = 12, title = "Summary", status = "info", solidHeader = TRUE,
          #              verbatimTextOutput("model_summary"))
          # ),
          shiny::tabPanel("Model Summary",
                   uiOutput("model_summary_ui"),
                   uiOutput("model_error_ui")
          ),

          # Tab 3: Diagnostics
          shiny::tabPanel("Diagnostics",
                   tags$div(style = "margin-top: 20px;"),
                   shinydashboard::box(width = 12, title = "Diagnostic Plots", status = "warning", solidHeader = TRUE,
                       plotOutput("diagnostics", height = "600px"))
          ),

          # Tab 4: Model interpretation
          shiny::tabPanel("Results Visualization",
                   tags$div(style = "margin-top: 20px;"),
                   shinydashboard::box(width = 12, title = "Model predictions", status = "success", solidHeader = TRUE,
                       plotOutput("predictions_plot", height = "650px"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # Reactive components
  data_reactive <- eventReactive(input$load_data, {
    req(input$data_file)
    read.csv(input$data_file$datapath, sep = input$sep)
  })

  data_loaded <- reactiveVal(FALSE)
  model_ready <- reactiveVal(FALSE)
  best_model <- reactiveVal(NULL)
  best_model_name <- reactiveVal(NULL)

  # UI Placeholders
  observe({
    req(input$data_file)
    df <- read.csv(input$data_file$datapath, sep = input$sep)

    updateSelectInput(session, "response", choices = names(df))
    updateSelectInput(session, "props", choices = names(df))
    updateSelectInput(session, "time_col",
                      choices = c("Select column..." = "", names(df)),
                      selected = "")
    updateSelectInput(session, "response_id_col",
                      choices = c("Select column..." = "", names(df)),
                      selected = "")
    updateSelectInput(session, "plots",
                      choices = c("Select column..." = "NA", names(df)),
                      selected = "")
  })

  output$input_selections <- renderUI({
    output$param_selectors <- renderUI({
      tagList(

        h6(tags$b("Model specification")),
        # Model type
        checkboxGroupInput(
          "model_type",
          "Select model type(s) to fit",
          choices = c("Repeated measures" = "RM",
                      "Multivariate" = "MV")
          # no default selected
        ),

        # Variance–Covariance Structure
        tags$div(
          title = "Select variance–covariance structure",
          style = "cursor: help;",

          # RM options
          conditionalPanel(
            condition = "input.model_type.indexOf('RM') >= 0",
            radioButtons(
              "vc_structure_RM",
              "Variance–Covariance Structure (RM)",
              choices = c(
                "Compound Symmetry (CS)" = "CS",
                "Autoregressive (AR1)" = "AR1",
                "Unstructured (UN)" = "UN"
              ),
              selected = "UN"
            )
          ),

          # MV options
          conditionalPanel(
            condition = "input.model_type.indexOf('MV') >= 0",
            radioButtons(
              "vc_structure_MV",
              "Variance–Covariance Structure (MV)",
              choices = c(
                "Compound Symmetry (CS)" = "CS",
                "Unstructured (UN)" = "UN"
              ),
              selected = "UN"
            )
          )
        ),

        tags$div(style = "margin-top: 30px;"),

        h6(tags$b("Data inputs")),
        # Prop, y and plot_IDS
        tags$div(
          title = "For wide data: choose all relevant columns | For long data: choose name–value pair",
          style = "cursor: help;",   # shows help cursor on hover
          selectInput(
            "response",
            "Response(s)",
            choices = c("Select column(s)..." = ""),
            multiple = TRUE
          )
        ),

        tags$div(
          title = "For repeated measures models, select the column indicating time.",
          style = "cursor: help;",   # shows help cursor on hover
          conditionalPanel(
            condition = "input.model_type.indexOf('RM') >= 0",
              selectInput(
              "time_col",
              "Time column (required for Repeated Measures)",
              choices = c("Select column..." = ""),
              multiple = FALSE
            )
          )
        ),

        tags$div(
          title = "For multivariate models, select the column indicating response names if data is in long-format.",
          style = "cursor: help;",
          conditionalPanel(
            condition = "input.model_type.indexOf('MV') >= 0",
            selectInput(
              "response_id_col",
              "Response indicator column (optional, only need if data in long format)",
              choices = c("Select column..." = ""),
              multiple = FALSE
            )
          )
        ),
        tags$div(
          title = "Columns containing species proportions",
          style = "cursor: help;",   # shows help cursor on hover
          selectInput(
            "props",
            "Species proportions",
            choices = c("Select column(s)..." = ""),
            multiple = TRUE
          )
        ),
        tags$div(
          title = "Column containing identifier for experimental unit",
          style = "cursor: help;",   # shows help cursor on hover
          selectInput(
            "plots",
            "Experimental unit identifier",
            choices = c("Select column..." = ""),
            multiple = FALSE
          )
        ),
        tags$div(style = "margin-top: 30px;"),

        h6(tags$b("DI model type")),
        # DI interaction structures
        tags$div(
          title = "Diversity-Interactions model structure",
          style = "cursor: help;",   # shows help cursor on hover
          checkboxGroupInput(
            "int_structures",
            "Choose DI models to fit",
            choices = c("Intercept only (STR)",
                        "Identities only (ID)",
                        "Average interaction (AV)",
                        "Functional group interactions (FG)",
                        "Additive species interactions (ADD)",
                        "Full pairwise interactions (FULL)"),
            inline = FALSE
          )
        ),
        tags$div(style = "margin-top: 30px;"),

        h6(tags$b("Additional customisation")),
        # FG and ID values
        textInput("FG_input", "Functional groups (optional; needed only for FG model)", placeholder = 'e.g., c("G", "G", "L", "L")'),
        textInput("ID_input", "Identity effect groups (optional)", placeholder = 'e.g., c("G1", "G1", "G2", "G2")'),
        # Extra formula
        textInput("formula_input", "Additional fixed effects (optional)", placeholder = "Specified as a formula e.g., ~ x1 + x2"),
        tags$div(style = "margin-top: 30px;"),

        h6(tags$b("Theta")),
        radioButtons(
          "theta_type",
          "Theta value",
          choices = c("Fixed numeric" = "num", "Estimate" = "est"),
          inline = FALSE
        ),

        conditionalPanel(
          condition = "input.theta_type == 'num'",
          numericInput("theta_val", "Numeric value (same across all models)", value = 1, min = 0.01, max = 1.5, step = 0.01)
        ),

        conditionalPanel(
          condition = "input.theta_type == 'est'",
          helpText("Theta will estimated using profile likelihood. If multiple DI models are selected, it will be estimated for the simplest interaction structure and reused across others.")
        ),

        # # JS logic to toggle visibility of AR1 option
        # tags$script(HTML("
        #   $(document).on('shiny:inputchanged', function(event) {
        #     if (event.name === 'model_type') {
        #       if (event.value === 'MV') {
        #         // Hide AR1 when Multivariate is selected
        #         $('input[name=\"vc_structure\"][value=\"AR1\"]').closest('div.radio').hide();
        #       } else {
        #         // Show AR1 again for Repeated measures
        #         $('input[name=\"vc_structure\"][value=\"AR1\"]').closest('div.radio').show();
        #       }
        #     }
        #   });
        # "))
      )
    })

  })

  # Input validation
  validate_inputs <- function(data, inputs) {

    errors <- list()

    # Ensure at least one model type selected
    if (is.null(inputs$model_type) || length(inputs$model_type) == 0) {
      errors$model_type <- "Please select at least one model type (Repeated Measures and/or Multivariate)."
    }

    # Ensure at least one DI model is selected
    if (is.null(inputs$int_structures) || length(inputs$int_structures) == 0) {
      errors$int_structures <- "Please select at least one DI model to fit."
    }

    # Ensure things aren't empty
    if (is.null(inputs$response) || length(inputs$response) == 0) {
      errors$response <- "No response variable(s) selected."
    }

    if (is.null(inputs$props) || length(inputs$props) == 0) {
      errors$props <- "No species proportion column(s) selected."
    }

    if (inputs$plots == "" || length(inputs$plots) == 0) {
      errors$plots <- "No experimental unit identifier selected."
    }

    # Ensure time column is specified if RM is selected
    if ("RM" %in% inputs$model_type) {
      if (inputs$time_col == "" || length(inputs$time_col) == 0) {
        errors$time_col <- "For repeated measures models, a time column should be selected."
      }
    }

    # Ensure time column is specified if MV is selected and only one response column selected
    if ("MV" %in% inputs$model_type && length(inputs$response) == 1) {
      if (inputs$response_id_col == "" || length(inputs$response_id_col) == 0) {
        errors$response_id_col <- "For multivariate models, if specifying a single response column, then a column containing response names should be selected."
      }
    }

    # Ensure columns exist
    # all_selected <- unlist(c(inputs$response, inputs$props, inputs$plots))
    # missing_cols <- setdiff(all_selected, names(data))
    # if (length(missing_cols) > 0) {
    #   errors$missing <- paste("The following columns are missing in the data:",
    #                           paste(missing_cols, collapse = ", "))
    # }

    # Ensure responses are numeric
    if (length(inputs$response) > 0) {
      # if(length(inputs$response) > 2) {
      non_numeric <- sapply(inputs$response, function(col) !is.numeric(data[[col]]))
      if (any(isTRUE(non_numeric))) {
        errors$numeric_res <- paste("For data in wide format (i.e., specifying more than two columns), all response columns must be numeric. These columns aren't: ",
                                    paste(names(non_numeric), collapse = ", "))
      # }
      }
      # else if (length(inputs$response) == 2) {
      #   non_numeric <- sapply(inputs$response, function(col) !is.numeric(data[[col]]))
      #   if (all(isTRUE(non_numeric))) {
      #     errors$numeric_res <- paste("For data in long format (i.e., specifying name-value pairs), both response columns can't be characters.")
      #   }
      # }
      # else {
      #   errors$only_one <- paste("Specify additional responses")
      # }
    }

    # Ensure proportions are numeric
    if (length(inputs$props) > 0) {
      non_numeric <- sapply(inputs$props, function(col) !is.numeric(data[[col]]))
      if (any(isTRUE(non_numeric))) {
        errors$numeric_prop <- paste("Species proportions must be numeric. These columns aren't:",
                                     paste(names(non_numeric), collapse = ", "))
      }
    }

    # Ensure proportions sum to 1
    if (length(inputs$props) > 0) {
      oob <- sapply(inputs$props, function(col) {
        any(data[[col]] < 0 | data[[col]] > 1, na.rm = TRUE)
      })
      if (any(oob)) {
        errors$bounds <- paste("Species proportions should be between 0 and 1. Currently these rows have values beyond these bounds.",
                               paste(names(oob)[oob], collapse = ", "))
      }

      if(!is.null(errors$numeric_prop)){
        sums <- rowSums(data[, inputs$props])

        if (!all((near(sums, 1, tol = 1e-3) | sums==0))) {
          errors$sums <- paste("Species proportions should sum to 1. Currently certain rows do not sum to 1.")
        }
      }
    }

    # Ensure FG and ID are characters
    if (!is.null(inputs$FG_input) && all(inputs$FG_input != "")) {
      if(length(inputs$FG_input) != length(inputs$props)){
        errors$FG_err <- "Functional grouping should be specified as a character vector with same length as number of species."
      }
    }

    if (!is.null(inputs$ID_input) && all(inputs$ID_input != "")) {
      if(length(inputs$ID_input) != length(inputs$props)){
        errors$ID_err <- "Identity grouping should be specified as a character vector with same length as number of species."
      }
    }

    # Ensure extra-formula contains columns in data
    if (!is.null(inputs$formula_input) && inputs$formula_input != "") {
      # Try to parse the formula
      formula_result <- tryCatch(
        as.formula(inputs$formula_input),
        error = function(e) return(e)
      )

      if (inherits(formula_result, "error")) {
        errors$formula <- paste0("Invalid formula syntax due to following error: ", formula_result$message,
                                 ".\nPlease specify formula with correct syntax, i.e., ~ x1 + x2")
      } else {
        # Formula valid so now check to ensure terms are present
        vars_in_formula <- all.vars(terms(formula_result))
        missing_in_formula <- setdiff(vars_in_formula, names(data))
        if (length(missing_in_formula) > 0) {
          errors$formula_vars <- paste(
            "The following variables in the formula are not in the data:",
            paste(missing_in_formula, collapse = ", ")
          )
        }
      }
    }


    # Return results
    if (length(errors) > 0) {
      return(list(valid = FALSE, messages = errors))
    } else {
      return(list(valid = TRUE))
    }
  }

  # Fitted models
  models_fitted <- eventReactive(input$run_model, {
    data <- data_reactive()
    inputs <- list(
      response = input$response,
      time_col = input$time_col,
      response_id_col = input$response_id_col,
      props = input$props,
      plots = input$plots,
      model_type = input$model_type,
      vc_structure_RM = input$vc_structure_RM,
      vc_structure_MV = input$vc_structure_MV,
      int_structures = input$int_structures,
      FG_input = eval(parse(text = input$FG_input)),
      ID_input = eval(parse(text = input$ID_input)),
      formula_input = input$formula_input,
      theta_type = input$theta_type,
      theta_val = input$theta_val
    )

    # Input validation
    validation <- validate_inputs(data, inputs)
    if (!validation$valid) {
      showModal(modalDialog(
        title = "Input Validation Error",
        easyClose = TRUE,
        footer = modalButton("OK"),
        tagList(
          lapply(validation$messages, function(msg) {
            tags$p(style = "color:red;", paste("•", msg))
          })
        )
      ))
      return(NULL)  # Stop model fitting
    }

    is_discrete <- function(x) {
      is.numeric(x) && all(!is.na(x)) && all(x == floor(x))
    }

    # Model fitting if inputs were okay
    withProgress(message = "Fitting selected models...", value = 0, {
      results <- list()
      errors <- list()
      DI_codes <- sub(".*\\(([^)]+)\\).*", "\\1", inputs$int_structures)

      total_models <- length(inputs$int_structures)
      for (i in seq_along(inputs$int_structures)) {
        DI_mod <- DI_codes[i]
        incProgress(1 / total_models, detail = paste("Fitting", DI_mod))

        # Prepare inputs as needed by DImulti
        # Responses
        # long_flag <- FALSE
        # long_name <- "NA"
        # # Wide format
        # if(length(inputs$response) == 2){
        #   # Wide format but with only two responses
        #   if(all(sapply(data[, inputs$response], is.numeric)) &
        #      all(sapply(data[, inputs$response], function(x) is_discrete(x)))){
        #     ys <- inputs$response
        #   # Long format
        #   } else {
        #     long_flag <- TRUE
        #     ys <- inputs$response[sapply(data[, inputs$response], function(x) !(is.character(x) | is_discrete(x)))]
        #     long_name <- inputs$response[inputs$response != ys]
        #   }
        # } else {
          ys <- inputs$response
        # }

        # MV or RM flags
        time_val <- c("NA", "NA")
        func_val <- c("NA", "NA")
        if ("RM" %in% inputs$model_type) {
          time_val <- c(inputs$time_col, inputs$vc_structure_RM)
        }
        if ("MV" %in% inputs$model_type) {
          resp_id <- ifelse(inputs$response_id_col == "", "NA", inputs$response_id_col)
          func_val <- c(resp_id, inputs$vc_structure_MV)
        }

        # Safe model fitting
        fit_result <- tryCatch({
          suppressWarnings(DImulti(
            y = ys,
            eco_func = func_val,
            time = time_val,
            unit_IDs = inputs$plots,
            prop = inputs$props,
            data = data,
            DImodel = DI_mod,
            FG = if (all(inputs$FG_input == "")) NULL else inputs$FG_input,
            ID = if (all(inputs$ID_input == "")) NULL else inputs$ID_input,
            extra_fixed = if (inputs$formula_input == "") NULL else as.formula(inputs$formula_input),
            estimate_theta = ifelse(inputs$theta_type == "num", FALSE, TRUE),
            theta = ifelse(!is.null(inputs$theta_val), as.numeric(inputs$theta_val), 1),
            method = "ML",
            method_theta = "univariate"
          ))
        }, error = function(e) {
          e
        })

        if (inherits(fit_result, "error")) {
          errors[[DI_mod]] <- fit_result$message
        } else {
          results[[DI_mod]] <- fit_result
        }
      }

      model_ready(TRUE)
      # Store results
      output$model_error_ui <- renderUI({
        if (length(errors) > 0) {
          shinydashboard::box(
            width = 12,
            title = "Model Fitting Errors",
            status = "danger",
            solidHeader = TRUE,
            div(
              style = "background-color: #fdecea; color: #b71c1c; padding: 10px; border-radius: 8px;",
              tags$ul(
                lapply(names(errors), function(name) {
                  tags$li(HTML(paste0("<b>", name, ":</b> ", errors[[name]])))
                })
              )
            )
          )
        } else {
          NULL
        }
      })

      # Display outputs
      if (length(results) == 0) {
        output$model_comparison_ui <- renderUI({
          tags$p("No models were successfully fitted.", style = "color:red;")
        })
        output$model_summary <- renderText("")
        return(NULL)
      }

      # browser()
      # If multiple models fitted
      if (length(results) > 1) {
        results_df <- eval(parse(text=paste("anova(",
                                            paste("results[[",1:length(results),"]]",
                                                  sep="",collapse=","),")")))
        results_df$AICc <- sapply(results, AICc)
        rownames(results_df) <- names(results)
        results_df <- dplyr::select(as.data.frame(results_df), Model, df, AIC, BIC, AICc, logLik, Test, L.Ratio, `p-value`)

        # browser()
        test_names <- rep("", nrow(results_df))
        fit_mod_codes <- names(results)
        for(i in 1:length(fit_mod_codes)){
          if (i == 1){
            test_names[i] <- ""
          } else {
            test_names[i] <- paste0(fit_mod_codes[i-1], " vs ", fit_mod_codes[i])
          }
        }
        test_names[which(test_names == "FG vs ADD")] <- "Not nested with FG"

        results_df <- dplyr::mutate(results_df, Test = test_names)
        results_df <- dplyr::mutate(results_df, Model = names(results),
                 AIC = round(AIC, 2),
                 BIC = round(BIC, 2),
                 AICc = round(AICc, 2),
                 logLik = round(logLik, 2),
                 L.Ratio = round(L.Ratio, 2),
                 `p-value` = round(`p-value`, 2))

        # Mark best model (lowest AICc)
        best_model <- results_df$Model[which.min(results_df$AICc)]
        # Best model name
        best_model_name(best_model)

        # Re-fit best model with REML
        # Safe model fitting
        best_REML <- tryCatch({
          suppressWarnings(DImulti(
            y = ys,
            eco_func = func_val,
            time = time_val,
            unit_IDs = inputs$plots,
            prop = inputs$props,
            data = data,
            DImodel = best_model,
            FG = if (all(inputs$FG_input == "")) NULL else inputs$FG_input,
            ID = if (all(inputs$ID_input == "")) NULL else inputs$ID_input,
            extra_fixed = if (inputs$formula_input == "") NULL else as.formula(inputs$formula_input),
            estimate_theta = ifelse(inputs$theta_type == "num", FALSE, TRUE),
            theta = ifelse(!is.null(inputs$theta_val), as.numeric(inputs$theta_val), 1),
            method = "REML",
            method_theta = "univariate"
          ))
        }, error = function(e) {
          e
        })

        # Save best model
        best_model(best_REML)

        output$model_comparison_ui <- renderUI({
          best_row <- which.min(results_df$AICc)

          # results_df <- lapply(seq_along(results_df), function(i) {
          #   vals <- as.character(results_df[[i]])
          #   vals[best_row] <- paste0("<b>", vals[best_row], "</b>")
          #   vals
          # }) %>% do.call(cbind, .) %>%
          #   `colnames<-`(colnames(results_df))

          colnames <- colnames(results_df)
          results_df <- lapply(seq_along(results_df), function(i) {
            vals <- as.character(results_df[[i]])
            vals[best_row] <- paste0("<b>", vals[best_row], "</b>")
            vals
          })
          results_df <- do.call(cbind, results_df)
          colnames(results_df) <- colnames


          # browser()
          HTML(paste0(
            "<div style='max-height:400px; overflow-y:auto; margin-bottom: 10px;'>",
            "<table style='border-collapse:collapse; width:100%; text-align:center;'>",
            "<thead style='background-color:#f8f9fa; font-weight:bold;'>",
            paste0("<tr>", paste0("<th>", colnames(results_df), "</th>", collapse = ""), "</tr>"),
            "</thead><tbody>",
            paste(apply(results_df, 1, function(row)
              paste0("<tr>", paste0("<td style='padding:5px; border-bottom:1px solid #ddd;'>", row, "</td>", collapse = ""), "</tr>")
            ), collapse = ""),
            "</tbody></table></div>",
            "<i style = 'margin-bottom: 5px;'>Maximum Likelihood (ML) estimation was used for model comparison as these models had different fixed effects.</i>",
            "<br>",
            "<i>Model with lowest AICc value highlighted in bold.</i>"
          ))
        })

        output$model_summary <- renderPrint({
          print(best_REML)
        })

        output$selected_model_name <- renderText({
          req(best_model_name())  # best_model_name() is a reactiveVal containing the model name
          paste("Selected model:", best_model_name())
        })

        output$REML_text <- renderUI({
          tags$i("Note: The final selected model has been re-estimated using Restricted maximum likelihood (REML) to obtain unbiased variance estimates.")
        })
      } else {
        # Only one model
        single_name <- names(results)
        best_model_name(single_name)             # store name
        # Re-fit best model with REML
        # Safe model fitting
        best_REML <- tryCatch({
          suppressWarnings(DImulti(
            y = ys,
            eco_func = func_val,
            time = time_val,
            unit_IDs = inputs$plots,
            prop = inputs$props,
            data = data,
            DImodel = single_name,
            FG = if (all(inputs$FG_input == "")) NULL else inputs$FG_input,
            ID = if (all(inputs$ID_input == "")) NULL else inputs$ID_input,
            extra_fixed = if (inputs$formula_input == "") NULL else as.formula(inputs$formula_input),
            estimate_theta = ifelse(inputs$theta_type == "num", FALSE, TRUE),
            theta = ifelse(!is.null(inputs$theta_val), as.numeric(inputs$theta_val), 1),
            method = "REML",
            method_theta = "univariate"
          ))
        }, error = function(e) {
          e
        })

        # Save best model
        best_model(best_REML)                 # set as best model

        output$model_comparison_ui <- renderUI({
          tags$p(paste("Only one model fitted: ", single_name))
        })
        output$model_summary <- renderPrint({
          print(best_REML)
        })
        output$selected_model_name <- renderText({
          req(best_model_name())  # best_model_name() is a reactiveVal containing the model name
          paste("Selected model:", best_model_name())
        })
        output$REML_text <- renderUI({
          tags$i("Note: This model has been estimated using Restricted maximum likelihood (REML) to obtain unbiased variance estimates.")
        })

      }
    })
  })


  # Events
  observeEvent(input$load_data, {
    req(input$data_file)
    data_loaded(TRUE)
    model_ready(FALSE)  # reset model fitted status when new data loaded
  })

  observeEvent(input$run_model, {
    req(data_loaded())
    models_fitted()
  })

  observeEvent(input$data_file, {
    data_loaded(FALSE)
    model_ready(FALSE)
  })


  # Outputs
  output$model_summary_ui <- renderUI({

    if (!(model_ready())) {
      # Show an informative message in the UI (not console)
      tags$div(
        style = "padding:10px; color: #666; font-style: italic;",
        "Model summary will appear here after fitting.",
        tags$br(),
        "Click 'Fit Model(s)' to start."
      )
    } else {
      # req(results_available())  # only render when results exist

      tagList(
        tags$div(style = "margin-top: 20px;"),
        shinydashboard::box(
          width = 12, title = "Model Comparison", status = "info", solidHeader = TRUE,
          uiOutput("model_comparison_ui")
        ),
        tags$div(style = "margin-top: 20px;"),
        shinydashboard::box(
          width = 12, title = "Selected Model Summary", status = "primary", solidHeader = TRUE,
          tags$h5(textOutput("selected_model_name"), style = "color:#2c3e50; font-weight:bold;"),
          tags$p(
            htmlOutput("REML_text"),
            style = "margin-bottom:5px; font-style: italic"
          ),
          verbatimTextOutput("model_summary")
        )
      )
    }
  })

  output$upload_status <- renderText({
    if (is.null(input$data_file)) {
      return("No file selected yet.")
    } else if (!data_loaded()) {   # before clicking 'Load Data'
      return(paste0("File selected: ", input$data_file$name, ". Click 'Load Data' to continue."))
    } else {
      return(paste("Data successfully loaded from", input$data_file$name))
    }
  })

  output$data_preview <- renderTable({
    if (!data_loaded()) {
      return(dplyr::tibble(
        " " = "Upload a data file (CSV), adjust the separator if needed, and click 'Load Data' to see a preview here."
      ))
    }
    head(data_reactive())
  })

  output$diagnostics <- renderPlot({
    if (!model_ready()) {
      pl <- ggplot2::ggplot(data = data.frame(x = 0.5, y = 0.5,
                                     label = "Diagnostics plots for the selected model will appear here after fitting model(s).")) +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = label), size = 6) +
        ggplot2::theme_void()
    } else {
      model <- best_model()
      req(model)
      withProgress(message = "Rendering diagnostic plot...", {
        pl <- DImodelsVis::model_diagnostics(model = model, which = 1:2,
                                only_extremes = TRUE)
      })
    }
    return(pl)
  })

  output$predictions_plot <- renderPlot({
    if (!model_ready()) {
      pl <- ggplot2::ggplot(data = data.frame(x = 0.5, y = 0.5,
                                     label = "Prediction plots for the selected model will appear here after fitting model(s).")) +
        ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = label), size = 6) +
        ggplot2::theme_void()
    } else {
      model <- best_model()
      req(model)
      withProgress(message = "Rendering predictions plot...", {
        pred_data <- DImodelsVis::gradient_change(model = model, plot = FALSE)
        # browser()
        attr(pred_data, "prop") <- NULL
        pl <- DImodelsVis::gradient_change_plot(data = pred_data, prop = NULL)
      })
    }

    return(pl)
  })

}

shinyApp(ui, server)
