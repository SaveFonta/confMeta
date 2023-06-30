# remotes::install_github("felix-hof/hMean", ref = "dev")
library(shiny)
library(hMean)
library(patchwork)

# Define a function that does the plots
make_plot <- function(pars) {
    with(
        pars,
        {
            pval <- ggPvalueFunction(
                thetahat = thetahat,
                se = se,
                level = level,
                heterogeneity = het,
                distr = distr,
                pValueFUN = pValueFUN,
                pValueFUN_args = list(
                    check_inputs = FALSE
                )
            )
            forest <- ForestPlot(
                thetahat = thetahat,
                se = se,
                level = level,
                distr = distr,
                pValueFUN = pValueFUN,
                heterogeneity = het
            )
            pval / forest
        }
    )
}

ui <- fluidPage(
    # Title bar
    titlePanel("Package hMean"),
    # Inputs
    fluidRow(
        column(
            12,
            selectInput(
                "fct1", "P-value function(s)",
                choices = c(
                    "hMean",
                    "k-Trials",
                    "Pearson",
                    "Edgington",
                    "Fisher"
                ),
                multiple = TRUE,
                selected = "hMean"
            ),
            conditionalPanel(
                'input.fct1.includes("hMean")',
                selectInput(
                    "distr1", "Harmonic Mean distribution",
                    choices = c("chisq", "f"),
                    multiple = TRUE,
                    selected = "chisq"
                )
            ),
            numericInput(
                "n1", "Number of studies", value = 3, min = 2, max = 50,
                step = 1
            ),
            numericInput("seed1", "RNG seed", value = 42),
            numericInput("mean1", "Effect: Mean", value = 0, step = 0.1),
            numericInput("sd1", "Effect: SD", value = 1, step = 0.1),
            numericInput("shape1", "SE: shape", value = 5, step = 0.1),
            numericInput("rate1", "SE: rate", value = 5, step = 0.1),
            selectInput(
                "het1", "Heterogeneity",
                choices = c("none", "additive", "multiplicative"),
                selected = "none",
                multiple = TRUE
            ),
            numericInput(
                "level1", "Confidence level",
                value = 0.95, min = 0.001, max = 0.999, step = 0.001
            ),
            sliderInput("xlim1", "Range of x-axis", value = c(-3, 3))
        ),
        # column(
        #     6,
        #     selectInput(
        #         "fct2", "P-value function(s)",
        #         choices = c(
        #             "hMean",
        #             "k-Trials",
        #             "Pearson",
        #             "Edgington",
        #             "Fisher"
        #         ),
        #         multiple = TRUE,
        #         selected = "hMean"
        #
        #     ),
        #     conditionalPanel(
        #         'input.fct2.includes("hMean")',
        #         selectInput(
        #             "distr2", "Harmonic Mean distribution",
        #             choices = c("chisq", "f"),
        #             multiple = TRUE,
        #             selected = "chisq"
        #         )
        #     ),
        #     numericInput(
        #         "n2", "Number of studies", value = 3, min = 2, max = 50,
        #         step = 1
        #     ),
        #     numericInput("seed2", "RNG seed", value = 42),
        #     numericInput("mean2", "Effect: Mean", value = 0, step = 0.1),
        #     numericInput("sd2", "Effect: SD", value = 1, step = 0.1),
        #     numericInput("shape2", "SE: shape", value = 5, step = 0.1),
        #     numericInput("rate2", "SE: rate", value = 5, step = 0.1),
        #     selectInput(
        #         "het2", "Heterogeneity",
        #         choices = c("none", "additive", "multiplicative"),
        #         selected = "none",
        #         multiple = TRUE
        #     ),
        #     numericInput(
        #         "level2", "Confidence level",
        #         value = 0.95, min = 0.001, max = 0.999, step = 0.001
        #     )
        # )
    ),

    # Action Button
    br(),
    actionButton("go", "Make the plots!"),
    br(),
    br(),
    br(),

    # Show outputs
    fluidRow(
        column(
            12,
            plotOutput("plot1")
        )#,
        # column(
        #     6,
        #     plotOutput("plot2")
        # )
    )
)

server <- function(input, output, session) {

    pars1 <- eventReactive(
        input$go,
        {
            set.seed(input$seed1)
            thetahat <- rnorm(
                n = input$n1,
                mean = input$mean1,
                sd = input$sd1
            )
            se <- rgamma(
                n = input$n1,
                shape = input$shape1,
                rate = input$rate1
            )
            het <- input$het1
            pValueFUN <- input$fct1
            d <- input$distr1
            distr <- if (!is.null(d)) "chisq" else d
            list(
                thetahat = thetahat,
                se = se,
                het = het,
                pValueFUN = pValueFUN,
                level = input$level1,
                distr = distr,
                xlim = input$xlim1
            )
        }
    )

    # pars2 <- eventReactive(
    #     input$go,
    #     {
    #         set.seed(input$seed2)
    #         thetahat <- rnorm(
    #             n = input$n2,
    #             mean = input$mean2,
    #             sd = input$sd2
    #         )
    #         se <- rgamma(
    #             n = input$n2,
    #             shape = input$shape2,
    #             rate = input$rate2
    #         )
    #         het <- input$het2
    #         pValueFUN <- input$fct2
    #         d <- input$distr2
    #         distr <- if (!is.null(d)) "chisq" else d
    #         list(
    #             thetahat = thetahat,
    #             se = se,
    #             het = het,
    #             pValueFUN = pValueFUN,
    #             level = input$level2,
    #             distr = distr
    #         )
    #     }
    # )

    # pars <- list(
    #     thetahat = rnorm(3),
    #     se = rgamma(3, 5, 5),
    #     level = 0.95,
    #     het = c("none", "additive"),
    #     pValueFUN = c("hMean"),
    #     distr = c("chisq", "f")
    # )
    # make_plot(pars)

    output$plot1 <- renderPlot(make_plot(pars1()))
    # output$plot2 <- renderPlot(make_plot(pars2()))
}

shinyApp(ui, server)
