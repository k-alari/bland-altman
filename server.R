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

shinyServer(function(input, output) {

    output$distPlot <- renderPlot({
        
        # Input values - regardless of prior specification method
        delta=input$delta
        n.samp=input$n.samp
        model.check=input$model.check
        
        
        #Convert inputted data into numeric vector
        data.receive=input$data.receive
        data.receive = as.character(data.receive)  #Make sure it is a character
        temp = gsub( " ", "", data.receive )  #Remove spaces
        d = as.numeric( unlist( strsplit( data.receive, split="," ) ) )  #Split by comma and list as numeric vector
        n = length(d)
        dbar = mean(d)
        v = sum((d-dbar)^2) / n

        # Function to obtain a0, b0, mu0, and lambda0 - For Normal
        parameters = function(sigma.hat, u.sigma, l.mu, u.mu) {
            
            # Distribution on sigma
            f = function(a0, b0, sigma){
                ifelse( exp(-b0*(1/sigma^2)) == 0, 0, 
                        ((b0^(a0+0.5))/(gamma(a0 + 0.5))) * ((1/(sigma^2)) ^ (a0 + 0.5-1)) * exp(-b0*(1/sigma^2)) * abs(-2/(sigma^3))
                )
            }
            
            
            # Estimate integral using Riemann sums
            n = 1000
            delta.x = u.sigma / n
            epsilon = 0.000001
            
            #Starting points
            delta = 0
            gamma = 101
            
            #To get right endpoints
            right=rep(NA, n)
            right[1] = delta.x
            for(i in 2:n){
                right[i] = right[i-1]+delta.x
            }
            
            # Storage to evaluate f at right endpoints
            f.sigma = rep(NA, length(right))
            
            
            # To obtain values for a0 and b0 
            for (K in 1:n) {
                a0 = (delta[K] + gamma[K])/2
                b0 = (a0 + 1) * (sigma.hat ^ 2)
                
                f.sigma = f(a0=a0, b0=b0, sigma = right)
                
                # If probability is too big (> 0.95)
                if(sum(f.sigma*delta.x, na.rm = TRUE) > 0.95){
                    delta[K+1] = delta[K]
                    gamma[K+1] = a0
                }
                
                # If probability is too small (< 0.95)
                if(sum(f.sigma*delta.x, na.rm = TRUE) < 0.95){
                    delta[K+1] = a0
                    gamma[K+1] = gamma[K]
                }
                
                if(gamma[K] - delta[K] < epsilon) break
            }
            
            # To get value for mu0
            mu0 = (u.mu + l.mu)/2
            
            # Now to get lambda0
            n = 100000
            delta=0
            gamma=100
            
            for(K in 1:n){
                lambda0 = (delta[K] + gamma[K])/2
                
                p = pt.scaled(u.mu, 2*a0, mu0, sqrt(b0/(a0*lambda0)), mu0) - pt.scaled(l.mu, 2*a0, mu0, sqrt(b0/(a0*lambda0)), mu0)
                
                if(p > 0.95){
                    delta[K+1] = delta[K]
                    gamma[K+1] = lambda0 } 
                
                if(p < 0.95){
                    delta[K+1] = lambda0
                    gamma[K+1] = gamma[K]
                }
                
                if((gamma[K] - delta[K]) < epsilon) break
            }
            
            values = list( a0, b0, mu0, lambda0 )
            names(values) = c( "a0", "b0", "mu0", "lambda0" )
            values
        }
        
        # Function that performs analysis and produces graphs
        # Norma-Gamma Prior
        BA.Bayesian = function( d, delta, a0, b0, mu0, lambda0, n.samp=10000, model.check="prop.agree" ) {
            n = length(d)
            dbar = mean(d)
            v = sum( (d-dbar)^2 ) / n
            lambda1 = lambda0 + n
            mu1 = ( lambda0*mu0+n*dbar ) / lambda1
            a1 = a0 + n/2
            b1 = b0 + n/2 * ( v+lambda0*(dbar-mu0)^2/lambda1 )
            mu.samp = tau.samp = rep( NA, n.samp )
            mu.samp[1] = dbar
            tau.samp[1] = 1/v
            
            set.seed(123) ### set seed for consistent result
            for ( i in 2:n.samp ) {
                mu.samp[i] = rnorm( 1, mu1, 1/sqrt(lambda1*tau.samp[i-1]) )
                tau.samp[i] = rgamma( 1, a1+0.5, b1+0.5*lambda1*(mu.samp[i]-mu1)^2 ) }
            sigma.samp = 1/sqrt(tau.samp)
            theta1.samp = mu.samp-1.96*sigma.samp
            theta2.samp = mu.samp+1.96*sigma.samp
            d.samp = rnorm( n.samp, mu.samp, sigma.samp )
            
            
            ##### HIGHEST DENSITY INTERVAL
            
            p.grid = round( seq( 0.01, 0.99, 0.01 ), 2 )
            n.grid = 0.05 * length(p.grid)
            
            mu.Q = quantile( mu.samp, p.grid )
            sigma.Q = quantile( sigma.samp, p.grid )
            theta1.Q = quantile( theta1.samp, p.grid )
            theta2.Q = quantile( theta2.samp, p.grid )
            diff.Q = quantile( d.samp, p.grid )
            
            rslt.mu.Q = rslt.sigma.Q = rslt.theta1.Q = rslt.theta2.Q = rslt.diff.Q = matrix( NA, n.grid, 3 )
            for ( i in 1:n.grid ) {
                index1 = which( p.grid == p.grid[i] ); index2 = which( p.grid == p.grid[i] + 0.95 )
                rslt.mu.Q[i,] = c( mu.Q[index1], mu.Q[index2], mu.Q[index2] - mu.Q[index1] )
                rslt.sigma.Q[i,] = c( sigma.Q[index1], sigma.Q[index2], sigma.Q[index2] - sigma.Q[index1] )
                rslt.theta1.Q[i,] = c( theta1.Q[index1], theta1.Q[index2], theta1.Q[index2] - theta1.Q[index1] )
                rslt.theta2.Q[i,] = c( theta2.Q[index1], theta2.Q[index2], theta2.Q[index2] - theta2.Q[index1] )
                rslt.diff.Q[i,] = c( diff.Q[index1], diff.Q[index2], diff.Q[index2] - diff.Q[index1] )
            }
            mu.HDI = rslt.mu.Q[ which.min( rslt.mu.Q[,3] ), 1:2 ]
            sigma.HDI = rslt.sigma.Q[ which.min( rslt.sigma.Q[,3] ), 1:2 ]
            theta1.HDI = rslt.theta1.Q[ which.min( rslt.theta1.Q[,3] ), 1:2 ]
            theta2.HDI = rslt.theta2.Q[ which.min( rslt.theta2.Q[,3] ), 1:2 ]
            diff.HDI = rslt.diff.Q[ which.min( rslt.diff.Q[,3] ), 1:2 ]
            
            p = c( 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 ) 
            mu.rslt = c( round(mean(mu.samp), 3), round(quantile( mu.samp, prob=p ),3), paste("(", round(mu.HDI[1],3), ", ", round(mu.HDI[2], 3), ")", sep="" ) )
            sigma.rslt = c( round(mean(sigma.samp),3), round(quantile( sigma.samp, prob=p ),3), paste("(", round(sigma.HDI[1],3), ", ", round(sigma.HDI[2],3), ")", sep="" ) )
            theta1.rslt = c( round(mean(theta1.samp),3), round(quantile( theta1.samp, prob=p ),3), paste("(", round(theta1.HDI[1],3), ", ", round(theta1.HDI[2],3), ")", sep="" ) )
            theta2.rslt = c( round(mean(theta2.samp),3), round(quantile( theta2.samp, prob=p ),3), paste("(", round(theta2.HDI[1],3), ", ", round(theta2.HDI[2],3), ")", sep="" ) )
            d.rslt = c( round(mean(d.samp),3), round(quantile( d.samp, prob=p ), 3), paste("(", round(diff.HDI[1], 3), ", ", round(diff.HDI[2], 3), ")", sep="") )
            post = rbind( mu.rslt, sigma.rslt, theta1.rslt, theta2.rslt, d.rslt )
            rownames(post) = c( "mu", "sigma", "theta1", "theta2", "difference" )
            colnames(post) = c( "mean", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "95% HDI")
            post.h1 = mean( theta1.samp >= -delta & theta2.samp <= delta )
            post.pred.agree = mean( abs(d.samp) <= delta )
            
            ### added: model checking
            stat.new = stat.new2 = rep( NA, n.samp )
            stat.obs = mean( ( d - mean(d) ) ^ 3 ) / var(d) ^ (3/2)
            stat.obs2 = mean( abs(d) < delta )
            
            for ( i in 1:n.samp ) { 
                d.new = sample( d.samp, size=n, replace=TRUE )
                stat.new[i] = mean( ( d.new - mean(d.new) ) ^ 3 ) / var(d.new) ^ (3/2)
                stat.new2[i] = mean( abs(d.new) < delta )
            }
            if ( input$model.check == "skewness" ) ppp = mean( stat.new > stat.obs )
            if ( input$model.check == "prop.agree" ) ppp = mean( stat.new2 > stat.obs2 )
            
            par( mfrow=c(2,3) )
            hist( mu.samp, xlab=expression(mu), col="red", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ mu), cex=0.9 )
            axis(1)
            hist( sigma.samp, xlab=expression(sigma), col="orange", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ sigma), cex=0.9 )
            axis(1)
            plot( mu.samp, sigma.samp, xlab=expression(mu), ylab=expression(sigma), main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("of" ~ mu ~ "and" ~ sigma), cex=0.9 )
            axis(1); axis(2)
            hist( theta1.samp, xlab=expression(theta[1]), col="darkgreen", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[1]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(lower bound)" ), cex=0.9 )
            axis(1)
            hist( theta2.samp, xlab=expression(theta[2]), col="blue", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[2]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(upper bound)"), cex=0.9 )
            axis(1)
            hist( d.samp, xlab=expression(tilde(D)), col="purple", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Predictive"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("Distribution of" ~ tilde(D)), cex=0.9 )
            axis(1)
            out = list( post, post.h1, post.pred.agree, ppp )
            names(out) = c("post","post.h1","post.pred.agree","ppp" )
            out
        }
        
        
        ### Independent Normal-Gamma Prior: calculating (a0, b0, mu0, lambda0) given (sigma.hat, u.sigma, l.mu, u.mu)
        parameters2 = function( sigma.hat, u.sigma, l.mu, u.mu ) {
            
            # Distribution on sigma
            f = function(a0, b0, sigma){
                ifelse( exp(-b0*(1/sigma^2)) == 0, 0, 
                        ((b0^(a0+0.5))/(gamma(a0 + 0.5))) * ((1/(sigma^2)) ^ (a0 + 0.5-1)) * exp(-b0*(1/sigma^2)) * abs(-2/(sigma^3))
                )
            }
            
            
            # Estimate integral using Riemann sums
            n = 1000
            delta.x = u.sigma / n
            epsilon = 0.000001
            
            #Starting points
            delta = 0
            gamma = 101
            
            #To get right endpoints
            right=rep(NA, n)
            right[1] = delta.x
            for(i in 2:n){
                right[i] = right[i-1]+delta.x
            }
            
            # Storage to evaluate f at right endpoints
            f.sigma = rep(NA, length(right))
            
            
            # To obtain values for a0 and b0 
            for (K in 1:n) {
                a0 = (delta[K] + gamma[K])/2
                b0 = (a0 + 1) * (sigma.hat ^ 2)
                
                f.sigma = f(a0=a0, b0=b0, sigma = right)
                
                # If probability is too big (> 0.95)
                if(sum(f.sigma*delta.x, na.rm = TRUE) > 0.95){
                    delta[K+1] = delta[K]
                    gamma[K+1] = a0
                }
                
                # If probability is too small (< 0.95)
                if(sum(f.sigma*delta.x, na.rm = TRUE) < 0.95){
                    delta[K+1] = a0
                    gamma[K+1] = gamma[K]
                }
                
                if(gamma[K] - delta[K] < epsilon) break
            }
            
            # To get value for mu0
            mu0 = (u.mu + l.mu)/2
            
            # Now to get lambda0
            n = 100000
            delta=0
            gamma=1000
            
            for(K in 1:n){
                lambda0 = (delta[K] + gamma[K])/2
                
                p = pnorm( u.mu, mu0, sqrt(1/lambda0) ) ### this is the modification for independent NG
                
                if(p > 0.95){
                    delta[K+1] = delta[K]
                    gamma[K+1] = lambda0 } 
                
                if(p < 0.95){
                    delta[K+1] = lambda0
                    gamma[K+1] = gamma[K]
                }
                
                if((gamma[K] - delta[K]) < epsilon) break
            }
            
            values = list( a0, b0, mu0, lambda0 )
            names(values) = c( "a0", "b0", "mu0", "lambda0" )
            values
        }
        
        # Independent Normal Gamma Prior Posterior Analysis
        BA.Bayesian2 = function( d, delta, a0, b0, mu0, lambda0, n.samp=10000, model.check="prop.agree" ) {
            n = length(d)
            dbar = mean(d)
            v = sum( (d-dbar)^2 ) / n
            a1 = a0 + n/2
            se = sqrt(v/n)
            mu.grid = seq( dbar-5*se, dbar+5*se, 10*se/1000 )
            mu.samp = tau.samp = rep( NA, n.samp )
            mu.samp[1] = dbar
            b1 = b0 + 0.5 * ( lambda0 * ( mu.samp[1] - mu0 )^2 + n * v + n * ( dbar - mu.samp[1] )^2 )
            tau.samp[1] = rgamma( 1, a1, b1 )
            
            set.seed(123) ### set seed for consistent result
            for ( i in 2:n.samp ) {
                b1.grid = b0 + 0.5 * ( lambda0 * ( mu.grid - mu0 )^2 + n * v + n * ( dbar - mu.grid )^2 )
                f.grid = exp( -tau.samp[i-1] * b1.grid )
                mu.samp[i] = sample( mu.grid, prob=f.grid, size=1 )
                b1 = b0 + 0.5 * ( lambda0 * ( mu.samp[i] - mu0 )^2 + n * v + n * ( dbar - mu.samp[i] )^2 )
                tau.samp[i] = rgamma( 1, a1, b1 ) }
            
            sigma.samp = 1/sqrt(tau.samp)
            theta1.samp = mu.samp-1.96*sigma.samp
            theta2.samp = mu.samp+1.96*sigma.samp
            d.samp = rnorm( n.samp, mu.samp, sigma.samp )
            
            
            ##### HIGHEST DENSITY INTERVAL
            
            p.grid = round( seq( 0.01, 0.99, 0.01 ), 2 )
            n.grid = 0.05 * length(p.grid)
            
            mu.Q = quantile( mu.samp, p.grid )
            sigma.Q = quantile( sigma.samp, p.grid )
            theta1.Q = quantile( theta1.samp, p.grid )
            theta2.Q = quantile( theta2.samp, p.grid )
            diff.Q = quantile( d.samp, p.grid )
            
            rslt.mu.Q = rslt.sigma.Q = rslt.theta1.Q = rslt.theta2.Q = rslt.diff.Q = matrix( NA, n.grid, 3 )
            for ( i in 1:n.grid ) {
                index1 = which( p.grid == p.grid[i] ); index2 = which( p.grid == p.grid[i] + 0.95 )
                rslt.mu.Q[i,] = c( mu.Q[index1], mu.Q[index2], mu.Q[index2] - mu.Q[index1] )
                rslt.sigma.Q[i,] = c( sigma.Q[index1], sigma.Q[index2], sigma.Q[index2] - sigma.Q[index1] )
                rslt.theta1.Q[i,] = c( theta1.Q[index1], theta1.Q[index2], theta1.Q[index2] - theta1.Q[index1] )
                rslt.theta2.Q[i,] = c( theta2.Q[index1], theta2.Q[index2], theta2.Q[index2] - theta2.Q[index1] )
                rslt.diff.Q[i,] = c( diff.Q[index1], diff.Q[index2], diff.Q[index2] - diff.Q[index1] )
            }
            mu.HDI = rslt.mu.Q[ which.min( rslt.mu.Q[,3] ), 1:2 ]
            sigma.HDI = rslt.sigma.Q[ which.min( rslt.sigma.Q[,3] ), 1:2 ]
            theta1.HDI = rslt.theta1.Q[ which.min( rslt.theta1.Q[,3] ), 1:2 ]
            theta2.HDI = rslt.theta2.Q[ which.min( rslt.theta2.Q[,3] ), 1:2 ]
            diff.HDI = rslt.diff.Q[ which.min( rslt.diff.Q[,3] ), 1:2 ]
            
            p = c( 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 ) 
            mu.rslt = c( round(mean(mu.samp), 3), round(quantile( mu.samp, prob=p ),3), paste("(", round(mu.HDI[1],3), ", ", round(mu.HDI[2], 3), ")", sep="" ) )
            sigma.rslt = c( round(mean(sigma.samp),3), round(quantile( sigma.samp, prob=p ),3), paste("(", round(sigma.HDI[1],3), ", ", round(sigma.HDI[2],3), ")", sep="" ) )
            theta1.rslt = c( round(mean(theta1.samp),3), round(quantile( theta1.samp, prob=p ),3), paste("(", round(theta1.HDI[1],3), ", ", round(theta1.HDI[2],3), ")", sep="" ) )
            theta2.rslt = c( round(mean(theta2.samp),3), round(quantile( theta2.samp, prob=p ),3), paste("(", round(theta2.HDI[1],3), ", ", round(theta2.HDI[2],3), ")", sep="" ) )
            d.rslt = c( round(mean(d.samp),3), round(quantile( d.samp, prob=p ), 3), paste("(", round(diff.HDI[1], 3), ", ", round(diff.HDI[2], 3), ")", sep="") )
            post = rbind( mu.rslt, sigma.rslt, theta1.rslt, theta2.rslt, d.rslt )
            rownames(post) = c( "mu", "sigma", "theta1", "theta2", "difference" )
            colnames(post) = c( "mean", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "95% HDI")
            post.h1 = mean( theta1.samp >= -delta & theta2.samp <= delta )
            post.pred.agree = mean( abs(d.samp) <= delta )
            
            ### added: model checking
            stat.new = stat.new2 = rep( NA, n.samp )
            stat.obs = mean( ( d - mean(d) ) ^ 3 ) / var(d) ^ (3/2)
            stat.obs2 = mean( abs(d) < delta )
            
            for ( i in 1:n.samp ) { 
                d.new = sample( d.samp, size=n, replace=TRUE )
                stat.new[i] = mean( ( d.new - mean(d.new) ) ^ 3 ) / var(d.new) ^ (3/2)
                stat.new2[i] = mean( abs(d.new) < delta )
            }
            if ( input$model.check == "skewness" ) ppp = mean( stat.new > stat.obs )
            if ( input$model.check == "prop.agree" ) ppp = mean( stat.new2 > stat.obs2 )
            
            par( mfrow=c(2,3) )
            hist( mu.samp, xlab=expression(mu), col="red", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ mu), cex=0.9 )
            axis(1)
            hist( sigma.samp, xlab=expression(sigma), col="orange", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ sigma), cex=0.9 )
            axis(1)
            plot( mu.samp, sigma.samp, xlab=expression(mu), ylab=expression(sigma), main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("of" ~ mu ~ "and" ~ sigma), cex=0.9 )
            axis(1); axis(2)
            hist( theta1.samp, xlab=expression(theta[1]), col="darkgreen", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[1]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(lower bound)" ), cex=0.9 )
            axis(1)
            hist( theta2.samp, xlab=expression(theta[2]), col="blue", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[2]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(upper bound)"), cex=0.9 )
            axis(1)
            hist( d.samp, xlab=expression(tilde(D)), col="purple", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Predictive"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("Distribution of" ~ tilde(D)), cex=0.9 )
            axis(1)
            out = list( post, post.h1, post.pred.agree, ppp )
            names(out) = c("post","post.h1","post.pred.agree","ppp" )
            out}
        
        
        ### Independent Flat Prior: posterior analysis
        BA.Bayesian3 = function( d, delta, l.sigma, u.sigma, l.mu, u.mu, n.samp=10000, model.check="prop.agree" ) {
            n = length(d)
            dbar = mean(d)
            v = sum( (d-dbar)^2 ) / n
            
            if ( l.mu > u.mu ) stop( "invalid boundaries for mu" )
            if ( l.sigma > u.sigma ) stop( "invalid boundaries for sigma" )
            if ( l.mu > dbar | u.mu < dbar ) stop( "boundaries for mu do not cover sample mean; prior may deviate too much from data" )
            if ( l.sigma > sqrt(v) | u.sigma < sqrt(v) ) stop( "boundaries for sigma do not cover sample standard deviation; prior may deviate too much from data" )
            
            l.sigma = max( 0, l.sigma ) ### in case user inputs a negative value
            l.tau = 1 / u.sigma^2
            u.tau = 1 / l.sigma^2
            
            mu.samp = tau.samp = rep( NA, n.samp )
            mu.samp[1] = mu.new = dbar
            tau.samp[1] = tau.new = 1/v
            
            set.seed(123) ### set seed for consistent result
            for ( i in 2:n.samp ) {
                mu.new = l.mu - 1
                while( mu.new < l.mu | mu.new > u.mu ) mu.samp[i] = mu.new = rnorm( 1, dbar, 1/sqrt(n*tau.new) )
                tau.new = l.tau - 1
                while( tau.new < l.tau | tau.new > u.tau ) tau.samp[i] = tau.new = rgamma( 1, 0.5*n + 1, 0.5*n*v + 0.5*n*(dbar-mu.new)^2 ) 
            }
            
            sigma.samp = 1/sqrt(tau.samp)
            theta1.samp = mu.samp-1.96*sigma.samp
            theta2.samp = mu.samp+1.96*sigma.samp
            d.samp = rnorm( n.samp, mu.samp, sigma.samp )
            
            
            ##### HIGHEST DENSITY INTERVAL
            
            p.grid = round( seq( 0.01, 0.99, 0.01 ), 2 )
            n.grid = 0.05 * length(p.grid)
            
            mu.Q = quantile( mu.samp, p.grid )
            sigma.Q = quantile( sigma.samp, p.grid )
            theta1.Q = quantile( theta1.samp, p.grid )
            theta2.Q = quantile( theta2.samp, p.grid )
            diff.Q = quantile( d.samp, p.grid )
            
            rslt.mu.Q = rslt.sigma.Q = rslt.theta1.Q = rslt.theta2.Q = rslt.diff.Q = matrix( NA, n.grid, 3 )
            for ( i in 1:n.grid ) {
                index1 = which( p.grid == p.grid[i] ); index2 = which( p.grid == p.grid[i] + 0.95 )
                rslt.mu.Q[i,] = c( mu.Q[index1], mu.Q[index2], mu.Q[index2] - mu.Q[index1] )
                rslt.sigma.Q[i,] = c( sigma.Q[index1], sigma.Q[index2], sigma.Q[index2] - sigma.Q[index1] )
                rslt.theta1.Q[i,] = c( theta1.Q[index1], theta1.Q[index2], theta1.Q[index2] - theta1.Q[index1] )
                rslt.theta2.Q[i,] = c( theta2.Q[index1], theta2.Q[index2], theta2.Q[index2] - theta2.Q[index1] )
                rslt.diff.Q[i,] = c( diff.Q[index1], diff.Q[index2], diff.Q[index2] - diff.Q[index1] )
            }
            mu.HDI = rslt.mu.Q[ which.min( rslt.mu.Q[,3] ), 1:2 ]
            sigma.HDI = rslt.sigma.Q[ which.min( rslt.sigma.Q[,3] ), 1:2 ]
            theta1.HDI = rslt.theta1.Q[ which.min( rslt.theta1.Q[,3] ), 1:2 ]
            theta2.HDI = rslt.theta2.Q[ which.min( rslt.theta2.Q[,3] ), 1:2 ]
            diff.HDI = rslt.diff.Q[ which.min( rslt.diff.Q[,3] ), 1:2 ]
            
            p = c( 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 ) 
            mu.rslt = c( round(mean(mu.samp), 3), round(quantile( mu.samp, prob=p ),3), paste("(", round(mu.HDI[1],3), ", ", round(mu.HDI[2], 3), ")", sep="" ) )
            sigma.rslt = c( round(mean(sigma.samp),3), round(quantile( sigma.samp, prob=p ),3), paste("(", round(sigma.HDI[1],3), ", ", round(sigma.HDI[2],3), ")", sep="" ) )
            theta1.rslt = c( round(mean(theta1.samp),3), round(quantile( theta1.samp, prob=p ),3), paste("(", round(theta1.HDI[1],3), ", ", round(theta1.HDI[2],3), ")", sep="" ) )
            theta2.rslt = c( round(mean(theta2.samp),3), round(quantile( theta2.samp, prob=p ),3), paste("(", round(theta2.HDI[1],3), ", ", round(theta2.HDI[2],3), ")", sep="" ) )
            d.rslt = c( round(mean(d.samp),3), round(quantile( d.samp, prob=p ), 3), paste("(", round(diff.HDI[1], 3), ", ", round(diff.HDI[2], 3), ")", sep="") )
            post = rbind( mu.rslt, sigma.rslt, theta1.rslt, theta2.rslt, d.rslt )
            rownames(post) = c( "mu", "sigma", "theta1", "theta2", "difference" )
            colnames(post) = c( "mean", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%", "95% HDI")
            post.h1 = mean( theta1.samp >= -delta & theta2.samp <= delta )
            post.pred.agree = mean( abs(d.samp) <= delta )
            
            ### added: model checking
            stat.new = stat.new2 = rep( NA, n.samp )
            stat.obs = mean( ( d - mean(d) ) ^ 3 ) / var(d) ^ (3/2)
            stat.obs2 = mean( abs(d) < delta )
            
            for ( i in 1:n.samp ) { 
                d.new = sample( d.samp, size=n, replace=TRUE )
                stat.new[i] = mean( ( d.new - mean(d.new) ) ^ 3 ) / var(d.new) ^ (3/2)
                stat.new2[i] = mean( abs(d.new) < delta )
            }
            if ( input$model.check == "skewness" ) ppp = mean( stat.new > stat.obs )
            if ( input$model.check == "prop.agree" ) ppp = mean( stat.new2 > stat.obs2 )
            
            par( mfrow=c(2,3) )
            hist( mu.samp, xlab=expression(mu), col="red", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ mu), cex=0.9 )
            axis(1)
            hist( sigma.samp, xlab=expression(sigma), col="orange", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ sigma), cex=0.9 )
            axis(1)
            plot( mu.samp, sigma.samp, xlab=expression(mu), ylab=expression(sigma), main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("of" ~ mu ~ "and" ~ sigma), cex=0.9 )
            axis(1); axis(2)
            hist( theta1.samp, xlab=expression(theta[1]), col="darkgreen", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[1]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(lower bound)" ), cex=0.9 )
            axis(1)
            hist( theta2.samp, xlab=expression(theta[2]), col="blue", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Distribution of" ~ theta[2]), cex=0.9 )
            mtext( side=3, line=0.5, adj=0.1, text=bquote("(upper bound)"), cex=0.9 )
            axis(1)
            hist( d.samp, xlab=expression(tilde(D)), col="purple", main="")
            mtext( side=3, line=1.5, adj=0.1, text=bquote("Posterior Predictive"), cex=0.9 )
            mtext( side=3, line=0, adj=0.1, text=bquote("Distribution of" ~ tilde(D)), cex=0.9 )
            axis(1)
            out = list( post, post.h1, post.pred.agree, ppp )
            names(out) = c("post","post.h1","post.pred.agree","ppp" )
            out}
        
        
        
        # Set prior values based on inputs 
        # Default Vague Prior - Normal Gamma
        if(input$prior=="Opt1"){
            a0=0.5
            b0=1e-6
            mu0=0
            lambda0=1e-6 
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian(d=d, delta=input$delta, a0=a0, b0=b0, mu0=mu0, lambda0=lambda0, n.samp=n.samp )
        }
        
        # Manual input of prior values - Normal Gamma
        if(input$prior=="Opt2"){
            a0=input$a0
            b0=input$b0
            mu0=input$mu0
            lambda0=input$lambda0
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian(d=d, delta=input$delta, a0=a0, b0=b0, mu0=mu0, lambda0=lambda0, n.samp=n.samp )
        }
        
        # Calculates prior values based on our prior specification - Normal Gamma 
        if(input$prior=="Opt3"){
            # Takes in input values
            sigma.hat=input$sigma.hat
            u.sigma=input$u.sigma
            l.mu=input$l.mu
            u.mu=input$u.mu
            
            # Calculates prior values based on input
            temp = parameters(sigma.hat = sigma.hat, u.sigma = u.sigma, l.mu = l.mu, u.mu = u.mu)
            a0 = temp[[1]]
            b0 = temp[[2]]
            mu0 = temp[[3]]
            lambda0 = temp[[4]]
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian(d=d, delta=input$delta, a0=a0, b0=b0, mu0=mu0, lambda0=lambda0, n.samp=n.samp )
        }
        
        # Manual input of prior values - Independent Normal Gamma
        if(input$prior=="Opt4"){
            a0=input$a0
            b0=input$b0
            mu0=input$mu0
            lambda0=input$lambda0
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian2(d=d, delta=input$delta, a0=a0, b0=b0, mu0=mu0, lambda0=lambda0, n.samp=n.samp )
        }
        
        # Calculates prior values based on our prior specification - Independent Normal Gamma 
        if(input$prior=="Opt5"){
            # Takes in input values
            sigma.hat=input$sigma.hat
            u.sigma=input$u.sigma
            l.mu=input$l.mu
            u.mu=input$u.mu
            
            # Calculates prior values based on input
            temp = parameters2(sigma.hat = sigma.hat, u.sigma = u.sigma, l.mu = l.mu, u.mu = u.mu)
            a0 = temp[[1]]
            b0 = temp[[2]]
            mu0 = temp[[3]]
            lambda0 = temp[[4]]
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian2(d=d, delta=input$delta, a0=a0, b0=b0, mu0=mu0, lambda0=lambda0, n.samp=n.samp )
        }
        
        # Uses independent uniform prior - with inputs of sigma.hat, u.sigma, etc.
        if(input$prior=="Opt6"){
            # Takes in input values
            l.sigma=input$l.sigma
            u.sigma=input$u.sigma
            l.mu=input$l.mu
            u.mu=input$u.mu
            
            # Output from BA.Bayesian function (using appropriate prior values based on user input with Normal-Gamma Prior)
            out = BA.Bayesian3(d=d, delta=input$delta, l.sigma=l.sigma, u.sigma=u.sigma, u.mu=u.mu, l.mu=l.mu, n.samp=n.samp )
        }
        
        # Output from BA.Bayesian function (using appropriate prior values based on user input)
        rslt = out$post
        ci = out$post[1,]
        
        # Pull out 95% CI for mu
        l = ci[2]
        u = ci[8]
        
        # Pull out Posterior Probability of H1
        posth1 = out$post.h1
        
        # Pull out Posterior Probability Distribution
        postpred = out$post.pred.agree
        
        # Pull out PPP
        
        ppp = out$ppp
        
        # Outputs for text interpretations. Follows HTML formatting 
        output$headingtext <- renderText({
            paste("<B>Based on the prior and observed data...</B>")
        })
        
        output$interpretation1 <- renderText({
            paste("<B>95% Credible Interval for μ:</B> The true mean (μ) of differences between the two measurements is between", l, "and", u, "with a probability of 0.95.")
        })
        
        output$interpretation2 <- renderText({
            paste("<B>Posterior Probability of H1:</B> The true limits of agreement (θ1 = μ - 1.96σ and θ2 = μ + 1.96σ) are within delta =", input$delta, "units with a probability of", round(posth1, 3),".")
        })
        
        output$interpretation3 <- renderText({
            paste("<B>Posterior Probability Distribution:</B> Future differences (between the two measurement methods) will be within delta =", input$delta, "units with a probability of", round(postpred, 3), ".")
        })
        
        output$interpretation4 <- renderText({
            paste("<B>Posterior Predictive P-value:</B> The posterior predictive p-value (PPP) is ", round(ppp, 3), ". A very small (i.e. close to 0) or large (i.e. close to 1) PPP indicates violation of the normality assumption.")
        })
     
        # Outputs rslt table using the DT package for R Shiny
        output$table <- DT::renderDataTable({
            DT::datatable(rslt)
        })

    })
    

})
