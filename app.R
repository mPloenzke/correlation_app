library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(tidyverse)

shinyApp(
    ui = dashboardPagePlus(
        header = dashboardHeaderPlus(
            enable_rightsidebar = TRUE,
            rightSidebarIcon = "gears",
            left_menu = tagList(
                dropdownButton(
                    label = "Controls",
                    icon = icon("sliders"),
                    status = "primary",
                    circle = FALSE,
                    sliderInput("p","Number in sensitive mixture:",min = 0,max = 100,value = 3,step=1),
                    sliderInput("var_y1","Covariance of sensitive mixture:",min = 0,max = 1,value = 0, step=.05),
                    sliderInput("var_y0","Covariance of resistant mixture:",min = 0,max = 1,value = 0, step=.05),
                    sliderInput("var_epsilon1","Noise of sensitive mixture:",min = 0,max = 1,value = .35, step=.05),
                    sliderInput("var_epsilon0","Noise of resistant mixture:",min = 0,max = 1,value = .65, step=.05),
                    sliderInput("max_signal","Maximum signal:",min = 1,max = 10,value = 5, step=1),
                    sliderInput("repetitions","Simulation repetitions:",min = 1,max = 25,value = 5, step=1)
                )
            )
        ),
        sidebar = dashboardSidebar(disable=TRUE),
        body = dashboardBody(
            titlePanel('Bivariate mixture of normals and correlation'),
            setShadow(class = "dropdown-menu"),
            fluidRow(
                splitLayout(cellWidths = c("49%", "49%"), 
                            tabsetPanel(type='tabs',
                                        tabPanel("Correlation measures for varying signal-to-noise regimes",plotOutput("correlationsPlot"))),
                            tabsetPanel(type='tabs',
                                        tabPanel("Simulation realizations", plotOutput("distributionsPlot")),
                                        tabPanel("Generating distributions", plotOutput("densitiesPlot"))))
            )
        ),
        rightsidebar = rightSidebar(),
        skin='purple',
    ),
    server = function(input, output) {
        twomixtures <- reactive({
            signal <- seq(0,input$max_signal,by=.25)
            repetitions <- input$repetitions
            p <- input$p/100
            mu_y0 <- 0
            var_y0 <- input$var_y0
            var_y1 <- input$var_y1
            var_epsilon0 <- input$var_epsilon0
            var_epsilon1 <- input$var_epsilon1
            noise <- var_epsilon0 + var_epsilon1
            fix_num_sensitive <- TRUE
            twomixtures <- lapply(1:length(signal), function(s2n) {
                mu_y1 <- signal[s2n]
                varcovar_matrix0 <- matrix(data=c(var_epsilon0+var_y0,var_y0,var_y0,var_epsilon0+var_y0),nrow=2,byrow=TRUE)
                varcovar_matrix1 <- matrix(data=c(var_epsilon1+var_y1,var_y1,var_y1,var_epsilon1+var_y1),nrow=2,byrow=TRUE)
                reprez <- tibble()
                for (rep in 1:repetitions) {
                    if (fix_num_sensitive) {
                        z <- c(rep(0,(1-p)*100),rep(1,p*100))
                    } else {
                        z <- rbinom(100,1,p)
                    }
                    y0 <- as.data.frame(mvtnorm::rmvnorm(100,mean=rep(mu_y0,2),sigma=varcovar_matrix0))
                    y1 <- as.data.frame(mvtnorm::rmvnorm(100,mean=rep(mu_y1,2),sigma=varcovar_matrix1))
                    z <- matrix(data=z,nrow=100,ncol=2,byrow=FALSE)
                    res <- (1-z)*y0 + z*y1
                    names(res) <- c('x','y')
                    q <- 1-p
                    covar_xx <- q*var_y0+q*mu_y0^2+p*mu_y1^2+p*var_y1-2*(q*p*mu_y0*mu_y1)-((q*mu_y0)^2)-((p*mu_y1)^2)
                    var_x <- q*((var_y0+var_epsilon0+mu_y0^2)-(var_y1+var_epsilon1+mu_y1^2)) + (var_y1+var_epsilon1+mu_y1^2) - (q*(mu_y0-mu_y1)+mu_y1)^2
                    res <- as_tibble(res) %>% 
                        mutate(signal=mu_y1-mu_y0,
                               mu_y1 = mu_y1,
                               noise=noise,
                               signal_to_noise = signal/noise,
                               repetition=rep,
                               mixture=z[,1],
                               dist = x + y,
                               #p_dist = dist>=nth(dist,sum(mixture),descending=TRUE),
                               p_dist = dist/2>=mu_y1/2,
                               spearman=cor(x,y,method='spearman'),
                               matthews = mltools::mcc(preds=p_dist,actuals=mixture),
                               pearson = covar_xx/(var_x))
                    reprez <- bind_rows(reprez,res)
                }
                reprez
            })
            bind_rows(twomixtures)
        })
        
        output$distributionsPlot <- renderPlot({
            rep2plot <- 1
            s2n_to_plot <- twomixtures() %>%
                distinct(signal_to_noise) %>%
                arrange() %>%
                pull() %>%
                quantile(probs=c(.2,.4,.6,.8)) %>% 
                as.vector()
            tdat <- twomixtures() %>%
                filter(repetition==rep2plot,
                       signal_to_noise %in% s2n_to_plot) %>%
                mutate(signal_to_noise = paste('Signal-to-noise: ',round(signal_to_noise,digits=2),sep='')) %>%
                tidyr::unite('called',c('mixture','p_dist'),sep='_',remove=FALSE)
            tdat_annotation <- tdat %>%
                group_by(signal_to_noise) %>%
                mutate(pearson = cor(x,y,method='pearson')) %>%
                ungroup() %>%
                mutate(maxx = -Inf, 
                       maxy = Inf,
                       pearson = as.character(round(pearson,digits=2))) %>%
                mutate(pearson = latex2exp::TeX(paste("$\\rho$ = ",pearson,sep=' '),output='character')) %>%
                distinct(signal_to_noise,maxx, maxy, pearson) %>%
                rename(x=maxx, y=maxy)
            p <- tdat %>%
                ggplot(aes(x=x,y=y)) + 
                geom_point(aes(color=as.factor(mixture)),show.legend = FALSE, alpha=.8) + 
                geom_vline(aes(xintercept=mu_y1/2),lty=2) +
                geom_hline(aes(yintercept=mu_y1/2),lty=2) +
                geom_text(data=tdat_annotation, aes(label=pearson,vjust=1.25,hjust=-.25), parse=TRUE) + 
                facet_wrap(vars(signal_to_noise),ncol=2) + 
                theme_bw() + 
                geom_smooth(method='lm',se=FALSE) + 
                xlab(latex2exp::TeX("$X_1$")) + 
                ylab(latex2exp::TeX("$X_2$")) + 
                scale_color_manual(values=c("black", "red","red","green")) + 
                theme(legend.position = 'none',text=element_text(size=18)) 
            p
        })
        
        output$correlationsPlot <- renderPlot({
            tdat <- twomixtures() %>%
                tidyr::gather('method','value',spearman:pearson) %>%
                mutate(method=as.factor(method))
            p2 <- tdat %>% 
                filter(method=='pearson') %>%
                distinct(signal_to_noise, value, method) %>%
                ggplot(aes(x=signal_to_noise,y=value, color=method)) + 
                geom_line(lwd=2,alpha=.75)+
                stat_smooth(data = tdat %>% filter(method != 'pearson'),aes(x=signal_to_noise,y=value, color=method), 
                            se=FALSE, lwd=2, method = 'gam', formula = y ~ s(x, bs = "cr")) + 
                theme_bw() + 
                scale_color_brewer(palette='Dark2',drop=FALSE) + 
                labs(x='Signal-to-noise',y='Correlation') + 
                theme(legend.position = 'top',
                      legend.title = element_blank(),
                      text=element_text(size=18))
            p2
        })
        
        output$densitiesPlot <- renderPlot({
            s2n_to_plot <- twomixtures() %>%
                distinct(signal_to_noise) %>%
                arrange() %>%
                pull() %>%
                quantile(probs=c(.2,.4,.6,.8)) %>% 
                as.vector()
            var_epsilon0 <- input$var_epsilon0
            var_epsilon1 <- input$var_epsilon1
            noise <- var_epsilon0 + var_epsilon1
            dd <- twomixtures() %>%
                filter(repetition==1,
                       signal_to_noise %in% s2n_to_plot) %>%
                select(mu_y1, signal_to_noise) %>%
                distinct() %>% 
                mutate(mu_y0 = 0,
                       var_y0 = var_epsilon0 + input$var_y0,
                       var_y1 = var_epsilon1 + input$var_y1)
            mixture0 <- dd %>%
                select(mu_y0, var_y0, signal_to_noise) %>%
                plyr::ddply(., "signal_to_noise", function(df) {
                    data.frame( 
                        value = seq(-3, input$max_signal, length = 100),
                        normal_curve = dnorm(seq(-3, input$max_signal, length = 100), 
                                           df$mu_y0, 
                                           df$var_y0)
                    )
                }) %>%
                mutate(mixture = 'mixture0',
                       normal_curve = normal_curve * (1-(input$p/100)))
            mixture1 <- dd %>%
                select(mu_y1, var_y1, signal_to_noise) %>%
                plyr::ddply(., "signal_to_noise", function(df) {
                    data.frame( 
                        value = seq(-3, input$max_signal, length = 100),
                        normal_curve = dnorm(seq(-3, input$max_signal, length = 100), 
                                             df$mu_y1, 
                                             df$var_y1) 
                    )
                }) %>%
                mutate(mixture = 'mixture1',
                       normal_curve = normal_curve * (input$p/100))
            bind_rows(mixture0,mixture1) %>% 
                filter(normal_curve > 1e-6) %>% 
                mutate(mu_y1 = signal_to_noise,
                       signal_to_noise = paste('Signal-to-noise: ',round(signal_to_noise,digits=2),sep='')) %>%
                ggplot(aes(x=value, y = normal_curve, fill=mixture)) +
                geom_line(color='black') + 
                geom_area(alpha=.8,position = "identity") + 
                geom_hline(aes(yintercept=0),lwd=1,color='black') + 
                geom_vline(aes(xintercept=mu_y1/2),lty=2,color='black') +
                facet_wrap(vars(signal_to_noise), ncol=2) + 
                theme_bw() +
                ylab("") +
                scale_y_continuous(breaks = NULL) +
                xlab(latex2exp::TeX("$Y$")) + 
                scale_color_manual(values=c("black", "red")) + 
                scale_fill_manual(values=c("black", "red")) + 
                theme(legend.position = 'none',text=element_text(size=18)) 
        })
    }
)
