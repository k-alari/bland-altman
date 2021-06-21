#
#  This R Shiny applet was designed as a companion to the following publication: 
#   
#    https://doi.org/10.1080/1091367X.2020.1853130
#
#
# Paper Abstract: There are two schools of thought in statistical analysis, frequentist,
# and Bayesian. Though the two approaches produce similar estimations and predictions in 
# large-sample studies, their interpretations are different. Bland Altman analysis is a 
# statistical method that is widely used for comparing two methods of measurement. It 
# was originally proposed under a frequentist framework, and it has not been used under 
# a Bayesian framework despite the growing popularity of Bayesian analysis. It seems 
# that the mathematical and computational complexity narrows access to Bayesian Bland 
# Altman analysis. In this article, we provide a tutorial of Bayesian Bland Altman 
# analysis. One approach we suggest is to address the objective of Bland Altman 
# analysis via the posterior predictive distribution. We can estimate the probability 
# of an acceptable degree of disagreement (fixed a priori) for the difference between 
# two future measurements. To ease mathematical and computational complexity, an 
# interface applet is provided with a guideline.

library(shiny)
library(DT)
library(metRology)
library(rmarkdown)

# Define UI for application
shinyUI(fluidPage(

    # Application title
    titlePanel("Bayesian Bland Altman Analysis"),

    # Get numeric inputs
    sidebarLayout(
        sidebarPanel(
            
            numericInput( inputId = "delta",
                         label = "What is the acceptable degree of agreement between the two measurements, δ > 0?",
                         value = 0.1, min = 0), 

            # Receive raw data 
            textInput( inputId = "data.receive", 
                       label="Input your sample data points (differences between two measurements) separated by commas (i.e., 0.02, 0.04, 0.05, ....).",
                       value="0.04, 0.09, 0.05, 0.1, 0.07, 0.05, 0.08, 0.06, 0.09, 0.03"),
            
            numericInput( inputId = "n.samp",
                          label="Input the size of a posterior sample.", 
                          value=10000
                          ),
            
            selectInput(inputId="model.check",
                        label="Which statistic do you want to use to calculate the posterior predictive p-value (PPP)?",
                        choices =c("Proportion of differences within delta (default)" = "prop.agree",
                                   "Skewness" = "skewness")
                        ),
            
            selectInput( inputId="prior", 
                         label="How do you want to specify the prior?", 
                         choices = c("Vague prior (default)" = "Opt1", 
                                     "Normal-gamma prior with a0, b0, mu0, and lambda0" = "Opt2", 
                                     "Normal-gamma prior with hat(sigma), u(sigma), l(mu), and u(mu)" = "Opt3",
                                     "Independent normal and gamma priors with a0, b0, mu0, and lambda0" = "Opt4",
                                     "Independent normal and gamma priors with hat(sigma), u(sigma), l(mu), and u(mu)" = "Opt5",
                                     "Independent uniform (flat) priors with hat(sigma), u(sigma), l(mu), and u(mu)" = "Opt6")),
            
            # Only shows up if user chooses Opt2: Normal-Gamma prior with a0, b0, etc.
            conditionalPanel(condition = "input.prior=='Opt2'",
                             numericInput( inputId="a0",
                                           label="What is your value for a0?", 
                                           value = 2.71),
                             numericInput( inputId="b0",
                                           label="What is your value for b0?", 
                                           value=0.00093),
                             numericInput( inputId="mu0",
                                           label="What is your value for mu0?", 
                                           value=0),
                             numericInput( inputId="lambda0",
                                           label="What is your value for lambda0?", 
                                           value=0.086)
            ),
            
            # Only shows up if user chooses Opt3: Normal-Gamma prior with hat(sigma), u(sigma), etc
            conditionalPanel(condition = "input.prior=='Opt3'",
                             numericInput( inputId="sigma.hat",
                                           label="What is your best guess for sigma, hat(sigma)?", 
                                           value=0.05),
                             numericInput( inputId="u.sigma",
                                           label="What is your upper bound for sigma, u(sigma)?", 
                                           value=0.1),
                             numericInput( inputId="l.mu",
                                           label="What is your lower bound for mu, l(mu)?",
                                           value=-0.5),
                             numericInput( inputId="u.mu",
                                           label="What is your upper bound for mu, u(mu)?",
                                           value=0.5)
            ),
            
            # Only shows up if user chooses Opt4: Independent Normal-Gamma with a0, b0, etc.
            conditionalPanel(condition = "input.prior=='Opt4'",
                             numericInput( inputId="a0",
                                           label="What is your value for a0?", 
                                           value = 2.71),
                             numericInput( inputId="b0",
                                           label="What is your value for b0?", 
                                           value=0.00093),
                             numericInput( inputId="mu0",
                                           label="What is your value for mu0?", 
                                           value=0),
                             numericInput( inputId="lambda0",
                                           label="What is your value for lambda0?", 
                                           value=0.086)
            ),
            
            # Only shows up if user chooses Opt5: Independent Normal-Gamma Prior with hat(sigma), u(sigma), etc
            conditionalPanel(condition = "input.prior=='Opt5'",
                             numericInput( inputId="sigma.hat",
                                           label="What is your best guess for sigma, hat(sigma)?", 
                                           value=0.05),
                             numericInput( inputId="u.sigma",
                                           label="What is your upper bound for sigma, u(sigma)?", 
                                           value=0.1),
                             numericInput( inputId="l.mu",
                                           label="What is your lower bound for mu, l(mu)?",
                                           value=-0.5),
                             numericInput( inputId="u.mu",
                                           label="What is your upper bound for mu, u(mu)?",
                                           value=0.5)
            ),
            
            # Only shows up if user chooses Opt6: Independent Uniform Priors with hat(sigma), u(sigma), etc
            conditionalPanel(condition = "input.prior=='Opt6'",
                             numericInput( inputId="l.sigma",
                                           label="What is your lower bound for sigma, l(sigma)?", 
                                           value=0),
                             numericInput( inputId="u.sigma",
                                           label="What is your upper bound for sigma, u(sigma)?", 
                                           value=1),
                             numericInput( inputId="l.mu",
                                           label="What is your lower bound for mu, l(mu)?",
                                           value=-1),
                             numericInput( inputId="u.mu",
                                           label="What is your upper bound for mu, u(mu)?",
                                           value=1)
            ),
            
           
            submitButton( "Submit", icon("refresh") )
            
        ),
        
        
        
        ### OUTPUT

        # Show plots generated by the BA.Bayesian function and data table
        mainPanel(
            plotOutput("distPlot"),
            htmlOutput("headingtext"),
            htmlOutput("interpretation1"),
            htmlOutput("interpretation2"), 
            htmlOutput("interpretation3"),
            htmlOutput("interpretation4"),
            tabPanel("rslt", DT::dataTableOutput("table"))
        
        )
    )
))
