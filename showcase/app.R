# remotes::install_github("felix-hof/hMean", ref = "dev")
library(shiny)
library(hMean)

ui <- fluidPage(
    # Title bar
    titlePanel("Package hMean"),
    # Inputs
    fluidRow(
        column(
            6,
            selectInput(
                "fct1", "P-value function(s)",
                choices = c("hMean", "k-Trials", "Pearson", "Edgington"),
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
            numericInput("mean1", "Effect: Mean", value = 0, step = 0.1),
            numericInput("sd1", "Effect: SD", value = 1, step = 0.1),
            numericInput("shape1", "SE: shape", value = 1, step = 0.1),
            numericInput("rate1", "SE: rate", value = 1, step = 0.1),
            numericInput("seed1", "RNG seed", value = 42),
            selectInput(
                "het1", "Heterogeneity",
                choices = c("none", "additive", "multiplicative"),
                selected = "none",
                multiple = TRUE
            ),
            numericInput(
                "level1", "Confidence level",
                value = 0.05, min = 0.001, max = 0.999, step = 0.001
            )

        ),
        column(
            6,
            selectInput(
                "fct2", "P-value function(s)",
                choices = c("hMean", "k-Trials", "Pearson", "Edgington"),
                multiple = TRUE,
                selected = "hMean"

            ),
            conditionalPanel(
                'input.fct2.includes("hMean")',
                selectInput(
                    "distr2", "Harmonic Mean distribution",
                    choices = c("chisq", "f"),
                    multiple = TRUE,
                    selected = "chisq"
                )
            ),
            numericInput(
                "n2", "Number of studies", value = 3, min = 2, max = 50,
                step = 1
            ),
            numericInput("mean2", "Effect: Mean", value = 0, step = 0.1),
            numericInput("sd2", "Effect: SD", value = 1, step = 0.1),
            numericInput("shape2", "SE: shape", value = 1, step = 0.1),
            numericInput("rate2", "SE: rate", value = 1, step = 0.1),
            numericInput("seed2", "RNG seed", value = 42),
            selectInput(
                "het2", "Heterogeneity",
                choices = c("none", "additive", "multiplicative"),
                selected = "none",
                multiple = TRUE
            ),
            numericInput(
                "level2", "Confidence level",
                value = 0.05, min = 0.001, max = 0.999, step = 0.001
            )
        )
    ),

    # Action Button
    br(),
    actionButton("go", "Make the plots!"),
    br(),

    # Show outputs
    fluidRow(
        column(
            6,
            verbatimTextOutput("pars1")
        ),
        column(
            6,
            verbatimTextOutput("pars2")
        )
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
            list(
                thetahat = thetahat,
                se = se,
                het = het,
                pValueFUN = pValueFUN,
                level = input$level1
            )
        }
    )

    pars2 <- eventReactive(
        input$go,
        {
            set.seed(input$seed2)
            thetahat <- rnorm(
                n = input$n2,
                mean = input$mean2,
                sd = input$sd2
            )
            se <- rgamma(
                n = input$n2,
                shape = input$shape2,
                rate = input$rate2
            )
            het <- input$het2
            pValueFUN <- input$fct2
            list(
                thetahat = thetahat,
                se = se,
                het = het,
                pValueFUN = pValueFUN,
                input = input$level2
            )
        }
    )

    output$p1 <- renderPlot({
        with(
            pars1(),

        )
    })



}

shinyApp(ui, server)


    # output$pars1 <- renderPrint({
    #     pars1()
    # })
    # output$pars2 <- renderPrint({
    #     pars2()
    # })

    # phi <- if (heterogeneity == "multiplicative") {
    #     estimatePhi(thetahat = thetahat, se = se)
    # } else {
    #     NULL
    # }
    # tau2 <- if (heterogeneity == "additive") {
    #     estimateTau2(thetahat = thetahat, se = se)
    # } else {
    #     NULL
    # }
    # mu <- seq(
    #     min(thetahat) - 0.5 * max(se),
    #     max(thetahat) + 0.5 * max(se),
    #     length.out = 1e5
    # )
    # alpha <- 0.05
    # funs <- list(
    #     "pearson" = hMean::pPearsonMu,
    #     "hMean" = hMean::hMeanChiSqMu,
    #     "k-Trials" = hMean::kTRMu,
    #     "edgington" = hMean::pEdgingtonMu,
    #     "fisher" = hMean::pFisherMu
    # )
    # p_vals <- do.call(
    #     "cbind",
    #     lapply(
    #         funs,
    #         function(f) {
    #             f(
    #                 thetahat = thetahat,
    #                 mu = mu,
    #                 se = se,
    #                 heterogeneity = heterogeneity,
    #                 phi = phi,
    #                 tau2 = tau2,
    #                 check_inputs = FALSE
    #             )
    #         }
    #     )
    # )
    # cis <- lapply(
    #     funs,
    #     function(f) {
    #         hMeanChiSqCI(
    #             thetahat = thetahat,
    #             se = se,
    #             level = 1 - alpha,
    #             alternative = "none",
    #             pValueFUN = f,
    #             pValueFUN_args = list(
    #                 check_inputs = FALSE,
    #                 heterogeneity = heterogeneity,
    #                 phi = phi,
    #                 tau2 = tau2
    #             )
    #         )
    #     }
    # )
    # plot_res <- function(
    #     mu,
    #     p_vals,
    #     cis = NULL,
    #     barheight = 0.05
    # ) {
    #     opar <- par(no.readonly = TRUE)
    #     par(las = 1)
    #     matplot(
    #         mu,
    #         p_vals,
    #         type = "l", lty = 1, lwd = 3,
    #         ylab = "p-value function", xlab = expression(mu)
    #     )
    #     legend("topleft",
    #         col = c(1, 2, 3, 4, 5),
    #         lwd = 3,
    #         lty = 1,
    #         legend = c("Pearson", "hMean", "k-Trials", "Edgington", "Fisher"),
    #         bty = "n",
    #         cex = 2
    #     )
    #     abline(h = 0.05, lty = 2)
    #     if (!is.null(cis)) {
    #         cis <- lapply(cis, "[[", i = "CI")
    #         jitter_inc <- 0.001
    #         jitter <- jitter_inc
    #         for (j in seq_along(cis)) {
    #             x <- cis[[j]]
    #             y_horiz <- alpha + (-1)^j * jitter
    #             if (j %% 2 == 0L) jitter <- jitter + jitter_inc
    #             for (i in seq_len(nrow(x))) {
    #                 l <- x[i, "lower"]
    #                 u <- x[i, "upper"]
    #                 lty <- 1
    #                 lwd <- 2
    #                 segments( # horizontal
    #                     x0 = l,
    #                     x1 = u,
    #                     y0 = y_horiz,
    #                     y1 = y_horiz,
    #                     lty = lty,
    #                     lwd = lwd,
    #                     col = j
    #                 )
    #                 segments( # error bar left
    #                     x0 = l,
    #                     x1 = l,
    #                     y0 = y_horiz - barheight / 2,
    #                     y1 = y_horiz + barheight / 2,
    #                     lty = lty,
    #                     lwd = lwd,
    #                     col = j
    #                 )
    #                 segments( # error bar left
    #                     x0 = u,
    #                     x1 = u,
    #                     y0 = y_horiz - barheight / 2,
    #                     y1 = y_horiz + barheight / 2,
    #                     lty = lty,
    #                     lwd = lwd,
    #                     col = j
    #                 )
    #             }
    #         }
    #     }
    #     par(opar)
    # }
    # plot_res(mu = mu, p_vals = p_vals, cis = cis)
