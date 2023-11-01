# @knitr read_data_create_functions

rm(list = ls())
# setwd("~/GitHub/polls")

library(mvtnorm)
# library(knitr)

load("last_sim.RData")

logit <- function(x) log(x/(1-x))
inv_logit <- function(x) 1/(1 + exp(-x))


# Adding ME1 and ME2, NE1 NE2 to sim_forecast matrix and ev vector

ev <- ev_state[colnames(sim_forecast)]
ev["ME"] <- 1
ev["NE"] <- 2
ev <- c(ev, "DC" = 3, "ME1" = 1, "ME2" = 1, "NE1" = 1, "NE2" = 1, "NE3" = 1)

sim_forecast <- cbind(sim_forecast, # Taking the 2012 delta with ME and NE state vote + noise
                      "ME1" = sim_forecast[,"ME"] + rnorm(nrow(sim_forecast),  .03,.0075),
                      "ME2" = sim_forecast[,"ME"] + rnorm(nrow(sim_forecast), -.03,.0075),
                      "DC" = 0.80 + rnorm(nrow(sim_forecast), 0,.015), 
                      "NE1" = sim_forecast[,"NE"] + rnorm(nrow(sim_forecast),  .03, .0075),
                      "NE2" = sim_forecast[,"NE"] + rnorm(nrow(sim_forecast),  .08, .0075),
                      "NE3" = sim_forecast[,"NE"] + rnorm(nrow(sim_forecast), -.10, .0075))

ev <- ev[colnames(sim_forecast)]

Sigma <- cov(logit(sim_forecast))
mu <- colMeans(logit(sim_forecast))
names(mu) <- colnames(sim_forecast)


draw_samples <- function(clinton_states = NULL, trump_states = NULL, states = NULL, 
                         upper_clinton = NULL, lower_clinton = NULL, print_acceptance = FALSE, target_nsim = 1000){
    sim <- matrix(NA, nr = 1, nc = length(mu))
    n <- 0
    while(nrow(sim) < target_nsim){
        # randomly sample from the posterior distribution and reject when constraints are not met
        n <- n + 1
        proposals <- inv_logit(rmvnorm(1e5, mu, Sigma, method = "svd")) # "DC" is pretty much uncorrelated
        colnames(proposals) <- names(mu)
        if (!is.null(clinton_states)) proposals[which(proposals[,clinton_states] < .5)] <- NA
        if (!is.null(  trump_states)) proposals[which(proposals[,  trump_states] > .5)] <- NA
        if (!is.null(        states)){
            for (s in states){
                proposals[which(proposals[, s] > upper_clinton[s] | 
                                proposals[, s] < lower_clinton[s])] <- NA
            }
        }
        reject <- apply(proposals, 1, function(x) any(is.na(x)))
        sim <- rbind(sim, proposals[!reject,])
        if (nrow(sim) < target_nsim & nrow(sim)/(nrow(proposals)*n) < 1-99.99/100){
            stop(paste("rmvnorm() is working hard... but more than 99.99% of the samples are rejected; you should relax some contraints.", sep = ""))
        }
    }
    return(list("matrix" = sim[-1,], "acceptance_rate" = nrow(sim)/(nrow(proposals)*n)))
}

update_prob <- function(clinton_states = NULL, trump_states = NULL, clinton_scores_list = NULL, target_nsim = 1000, show_all_states = FALSE){
    states <- names(clinton_scores_list)
    lower_clinton <- sapply(clinton_scores_list, function(x) x[1]/100)
    upper_clinton <- sapply(clinton_scores_list, function(x) x[2]/100)
    sim <- draw_samples(clinton_states = clinton_states, trump_states = trump_states, states = states, 
                        upper_clinton = upper_clinton, lower_clinton = lower_clinton, 
                        target_nsim = target_nsim)
    ev_dist <- (sim[["matrix"]] > .5) %*% ev
    state_win <- colMeans(sim[["matrix"]] > .5)
    p <- mean(ev_dist >= 270)
    sd <- sqrt(p*(1-p)/length(ev_dist))
    if (show_all_states){
        cat("Pr(Clinton wins) by state, in %:\n")
        print(t(round(100*state_win)))
        cat("--------\n")
    }
    cat(paste("Pr(Clinton wins the electoral college) = ", round(100*p), "%\n[nsim = ", length(ev_dist), "; se = ", round(sd*100,1), "%]", sep = ""))
    if (show_all_states) cat("\n--------\n")
}

update_prob_edited <- function(clinton_states = NULL, trump_states = NULL, clinton_scores_list = NULL, target_nsim = 1000, show_all_states = FALSE){
    states <- names(clinton_scores_list)
    lower_clinton <- sapply(clinton_scores_list, function(x) x[1]/100)
    upper_clinton <- sapply(clinton_scores_list, function(x) x[2]/100)
    sim <- draw_samples(clinton_states = clinton_states, trump_states = trump_states, states = states, 
                        upper_clinton = upper_clinton, lower_clinton = lower_clinton, 
                        target_nsim = target_nsim)
    ev_dist <- (sim[["matrix"]] > .5) %*% ev
    state_win <- colMeans(sim[["matrix"]] > .5)
    p <- mean(ev_dist >= 270)
    sd <- sqrt(p*(1-p)/length(ev_dist))
    return(as.vector(state_win))
}

real_scores <-c("AK"=0,  "AL"=0, "AR"=0, "AZ"=0,  "CA"=1, "CO"=1,  "CT"=1,  "DE"=1, "FL"=0, "GA"=0,  "HI"=1, "IA"=0, "ID"=0,  "IL"=1, "IN"=0, "KS"=0, "KY"=0, "LA"=0,  "MA"=1,  "MD"=1,  "ME"=1, "MI"=0, "MN"=1, "MO"=0, "MS"=0, "MT"=0, "NC"=0, "ND"=0, "NE"=0, "NH"=1,  "NJ"=1, "NM"=1, "NV"=1,  "NY"=1, "OH"=0, "OK"=0, "OR"=1, "PA"=0,  "RI"=1, "SC"=0, "SD"=0, "TN"=0, "TX"=0, "UT"=0, "VA"=1, "VT"=1,  "WA"=1, "WI"=0, "WV"=0, "WY"=0, "ME1"=1, "ME2"=0,  "DC"=1, "NE1"=0, "NE2"=0, "NE3"=0)

brier <- function()
{
    win_vector <- update_prob_edited()
    print("Odds by state that Clinton will win:")
    print(win_vector)
    print("Brier score of predictions")
    print((win_vector-real_scores)^2)
}

brier_2 <- function(state1 = NULL, state2 = NULL)
{
    index1 <- which(names(real_scores) %in% state1)
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()

    win_vector <- update_prob_edited(state2)
    print(paste("The probability of Clinton winning", state1, "when Clinton wins", state2, "is", win_vector[index1]))
    if( real_scores[index1]==1 &&real_scores[index2]==1)
     {
        brier_score <- ((predictions_vector[index2]* win_vector[index1])-1)^2
        print(paste("The brier score is", brier_score))
    }
    else {
        brier_score <- ((predictions_vector[index2]* win_vector[index1])-0)^2
        print(paste("The brier score is", brier_score))
    }

    print(paste("The probability of Trump winning", state1, "when Clinton wins", state2, "is", (1-win_vector[index1])))
    if( real_scores[index1]==0 &&real_scores[index2]==1)
     {
        brier_score <- ((predictions_vector[index2]* (1-win_vector[index1]))-1)^2
        print(paste("The brier score is", brier_score))
    }
    else {
        brier_score <- ((predictions_vector[index2]* (1-win_vector[index1]))-0)^2
        print(paste("The brier score is", brier_score))
    }
    
    win_vector2 <- update_prob_edited(NULL, state2)
    print(paste("The probability of Clinton winning", state1, "when Trump wins", state2, "is",win_vector2[index1]))
    if( real_scores[index1]==1 &&real_scores[index2]==0)
        {
        brier_score <- (((1-predictions_vector[index2])* win_vector2[index1])-1)^2
        print(paste("The brier score is", brier_score))
    }
    else {
        brier_score <- (((1-predictions_vector[index2])* win_vector2[index1])-0)^2
        print(paste("The brier score is", brier_score))
    }
    print(paste("The probability of Trump winning", state1, "when Trump wins", state2, "is",(1-win_vector2[index1])))
    if( real_scores[index1]==0 &&real_scores[index2]==0)
        {
        brier_score <- (((1-predictions_vector[index2])* (1-win_vector2[index1]))-1)^2
        print(paste("The brier score is", brier_score))
    }
    else {
        brier_score <- (((1-predictions_vector[index2])* (1-win_vector2[index1]))-0)^2
        print(paste("The brier score is", brier_score))
    }   
    

}

#1-CC, 2-TC, 3-CT, 4-TT
brier_3 <- function(state1 = NULL, state2 = NULL, scenario = NULL)
{
    index1 <- which(names(real_scores) %in% state1)
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()
    win_vector <- update_prob_edited(state2)
    win_vector2 <- update_prob_edited(NULL, state2)

    if(scenario==1)
    {
        name1 = "Clinton"
        name2 = "Clinton"
    }
    else if(scenario==2)
    {
        name1 = "Trump"
        name2 = "Clinton"
    }
    else if(scenario==3)
    {
        name1 = "Clinton"
        name2 = "Trump"
    }
    else if(scenario==4)
    {
        name1 = "Trump"
        name2 = "Trump"
    }

    if(name1 == "Clinton")
    {
        print(paste("The probability", name1, "wins", state1, "is", predictions_vector[index1]))
    }
    else 
    {
        print(paste("The probability", name1, "wins", state1, "is", 1-predictions_vector[index1]))
    }
    if(name2 == "Clinton")
    {
        print(paste("The probability", name2, "wins", state2, "is", predictions_vector[index2]))
    }
    else 
    {
        print(paste("The probability", name2, "wins", state2, "is", 1-predictions_vector[index2]))
    }

    if(scenario==1)
    {
        print(paste("The probability of Clinton winning", state1, "when Clinton wins", state2, "is", win_vector[index1]))
        if( real_scores[index1]==1 &&real_scores[index2]==1)
        {
            brier_score <- (predictions_vector[index1] - 1)^2
            brier_score <- brier_score + (predictions_vector[index2] -1)^2
            brier_score <- brier_score + (win_vector[index1]-1)^2      
        }
        else {

            if(real_scores[index1]==1)
            {
                brier_score <- (predictions_vector[index1] - 1)^2
            }
            else
            {
                brier_score <- (predictions_vector[index1] - 0)^2
            }
            if(real_scores[index2]==1)
            {
                brier_score <- brier_score + (predictions_vector[index2] -1)^2
            }
            else
            {
                brier_score <- brier_score + (predictions_vector[index2] -0)^2
            }
            brier_score <- brier_score + (win_vector[index1]-0)^2
        }
    }
    if(scenario==2)
    {
        print(paste("The probability of Trump winning", state1, "when Clinton wins", state2, "is", (1-win_vector[index1])))
        if( real_scores[index1]==0 &&real_scores[index2]==1)
        {
            brier_score <- ((1-predictions_vector[index1]) - 1)^2
            brier_score <- brier_score + (predictions_vector[index2] -1)^2
            brier_score <- brier_score + ((1-win_vector[index1])-1)^2    
        }
        else {
            if(real_scores[index1]==1)
            {
                brier_score <- ((1-predictions_vector[index1]) - 1)^2
            }
            else
            {
                brier_score <- ((1-predictions_vector[index1]) - 0)^2
            }
            if(real_scores[index2]==1)
            {
                brier_score <- brier_score + (predictions_vector[index2] -1)^2
            }
            else
            {
                brier_score <- brier_score + (predictions_vector[index2] -0)^2
            }
            brier_score <- brier_score + ((1-win_vector[index1])-0)^2
        }
        }
    
    
    
    if(scenario==3)
    {
        
        print(paste("The probability of Clinton winning", state1, "when Trump wins", state2, "is",win_vector2[index1]))
        if( real_scores[index1]==1 &&real_scores[index2]==0)
            {
            brier_score <- (predictions_vector[index1] - 1)^2
            brier_score <- brier_score + ((1-predictions_vector[index2]) -1)^2
            brier_score <- brier_score + (win_vector2[index1]-1)^2  
        }
        else {
            if(real_scores[index1]==1)
            {
                brier_score <- (predictions_vector[index1] - 1)^2
            }
            else
            {
                brier_score <- (predictions_vector[index1] - 0)^2
            }
            if(real_scores[index2]==1)
            {
                brier_score <- brier_score + ((1-predictions_vector[index2]) -1)^2
            }
            else
            {
                brier_score <- brier_score + ((1-predictions_vector[index2]) -0)^2
            }
            brier_score <- brier_score + (win_vector2[index1]-0)^2
        }
        }
    
    
    if(scenario==4)
    {
        print(paste("The probability of Trump winning", state1, "when Trump wins", state2, "is",(1-win_vector2[index1])))
        if( real_scores[index1]==0 &&real_scores[index2]==0)
            {
            brier_score <- ((1-predictions_vector[index1]) - 1)^2
            brier_score <- brier_score + ((1-predictions_vector[index2]) -1)^2
            brier_score <- brier_score + ((1-win_vector2[index1])-1)^2  

        }
        else {
            if(real_scores[index1]==1)
            {
                brier_score <- ((1-predictions_vector[index1]) - 1)^2
            }
            else
            {
                brier_score <- ((1-predictions_vector[index1]) - 0)^2
            }
            if(real_scores[index2]==1)
            {
                brier_score <- brier_score + ((1-predictions_vector[index2]) -1)^2
            }
            else
            {
                brier_score <- brier_score + ((1-predictions_vector[index2]) -0)^2
            }
            brier_score <- brier_score + ((1-win_vector2[index1])-0)^2
        }
          
    
    }
    
    brier_score <- brier_score/3
    print(paste("The brier score of all these happening is", brier_score))
}



brier_2edited <- function(state2 = NULL)
{
    all_states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO",
                    "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK",
                    "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI",
                    "WV", "WY", "ME1", "ME2", "DC", "NE1", "NE2", "NE3")
    real_scores <-c("AK"=0,  "AL"=0, "AR"=0, "AZ"=0,  "CA"=1, "CO"=1,  "CT"=1,  "DE"=1, "FL"=0, "GA"=0,  "HI"=1, "IA"=0, "ID"=0,  "IL"=1, "IN"=0, "KS"=0, "KY"=0, "LA"=0,  "MA"=1,  "MD"=1,  "ME"=1, "MI"=0, "MN"=1, "MO"=0, "MS"=0, "MT"=0, "NC"=0, "ND"=0, "NE"=0, "NH"=1,  "NJ"=1, "NM"=1, "NV"=1,  "NY"=1, "OH"=0, "OK"=0, "OR"=1, "PA"=0,  "RI"=1, "SC"=0, "SD"=0, "TN"=0, "TX"=0, "UT"=0, "VA"=1, "VT"=1,  "WA"=1, "WI"=0, "WV"=0, "WY"=0, "ME1"=1, "ME2"=0,  "DC"=1, "NE1"=0, "NE2"=0, "NE3"=0)
    
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()
    brier_score <- 0
    win_vector <- update_prob_edited(state2)
    win_vector2 <- update_prob_edited(NULL, state2)
    for(i in 1:56)
    {

        if(i!=index2)
        {
            if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
            {
                
                if( real_scores[i]==1 &&real_scores[index2]==1)
                {
                    b <-((predictions_vector[index2]* win_vector[i])-1)^2
                    brier_score <- b +brier_score
                    cat(all_states[i], state2, "CC" ,predictions_vector[i], win_vector[i], predictions_vector[index2], real_scores[i], real_scores[index2], b, "\n")
                }
                else {
                    b<-((predictions_vector[index2]* win_vector[i])-0)^2
                    brier_score <- b +brier_score
                    cat(all_states[i], state2, "CC" ,predictions_vector[i], win_vector[i], predictions_vector[index2], real_scores[i], real_scores[index2], b, "\n")
                }

                if( real_scores[i]==0 &&real_scores[index2]==1)
                {
                    b<-((predictions_vector[index2]* (1-win_vector[i]))-1)^2
                    brier_score <- brier_score + b
                    cat(all_states[i], state2, "TC" ,(1-predictions_vector[i]), (1-win_vector[i]), predictions_vector[index2], real_scores[i], real_scores[index2], b, "\n")

                }
                else {
                    b<-((predictions_vector[index2]* (1-win_vector[i]))-0)^2
                    brier_score <- brier_score +b
                    cat(all_states[i], state2, "TC" ,(1-predictions_vector[i]), (1-win_vector[i]), predictions_vector[index2], real_scores[i], real_scores[index2], b, "\n")

                }
                
                
                if( real_scores[i]==1 &&real_scores[index2]==0)
                {
                    b<-(((1-predictions_vector[index2])* win_vector2[i])-1)^2
                    brier_score <- brier_score +b
                    cat(all_states[i], state2, "CT" ,predictions_vector[i], win_vector2[i], (1-predictions_vector[index2]), real_scores[i], real_scores[index2], b, "\n")
                }
                else {
                    b<-(((1-predictions_vector[index2])* win_vector2[i])-0)^2
                    brier_score <- brier_score +b
                    cat(all_states[i], state2, "CT" ,predictions_vector[i], win_vector2[i], (1-predictions_vector[index2]), real_scores[i], real_scores[index2], b, "\n")
                }
                if( real_scores[i]==0 &&real_scores[index2]==0)
                {
                    b<-(((1-predictions_vector[index2])* (1-win_vector2[i]))-1)^2
                    brier_score <- brier_score +b
                    cat(all_states[i], state2, "TT" ,(1-predictions_vector[i]), (1-win_vector2[i]), (1-predictions_vector[index2]), real_scores[i], real_scores[index2], b, "\n")
                }
                else {
                    b<-(((1-predictions_vector[index2])* (1-win_vector2[i]))-0)^2
                    brier_score <- brier_score +b
                    cat(all_states[i], state2, "TT" ,(1-predictions_vector[i]), (1-win_vector2[i]), (1-predictions_vector[index2]), real_scores[i], real_scores[index2], b, "\n")
                }
            }
        }
    }
    return(brier_score)
}

brier_3edited <- function(state2 = NULL, scenario= NULL)
{
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()
    win_vector <- update_prob_edited(state2)
    win_vector2 <- update_prob_edited(NULL, state2)
    brier_score <-0
    total <- 0

    all_states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO",
                    "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK",
                    "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI",
                    "WV", "WY", "ME1", "ME2", "DC", "NE1", "NE2", "NE3")

    for(i in 1:56)
    {
        
        if(i!=index2)
        {
            if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
            {
                cat(all_states[i], state2, "")
        if(scenario==1)
        {
                cat(predictions_vector[i], win_vector[i], predictions_vector[index2],"")
            if( real_scores[i]==1 &&real_scores[index2]==1)
            {
                b1 <- (predictions_vector[i] - 1)^2
                b2 <-   (predictions_vector[index2] -1)^2
                b3 <-   ((win_vector[i]*predictions_vector[index2])-1)^2      
            }
            else {

                if(real_scores[i]==1)
                {
                    b1 <- (predictions_vector[i] - 1)^2
                }
                else
                {
                    b1 <- (predictions_vector[i] - 0)^2
                }
                if(real_scores[index2]==1)
                {
                    b2 <-   (predictions_vector[index2] -1)^2
                }
                else
                {
                    b2 <-   (predictions_vector[index2] -0)^2
                }
                b3 <-  ((win_vector[i]*predictions_vector[index2])-0)^2
            }
        }
        if(scenario==2)
        {
            cat((1-predictions_vector[i]), (1-win_vector[i]), predictions_vector[index2],"")
            if( real_scores[i]==0 &&real_scores[index2]==1)
            {
                b1 <- ((1-predictions_vector[i]) - 1)^2
                b2 <-  (predictions_vector[index2] -1)^2
                b3 <- (((1-win_vector[i])*predictions_vector[index2])-1)^2    
            }
            else {
                if(real_scores[i]==1)
                {
                    b1 <- ((1-predictions_vector[i]) - 1)^2
                }
                else
                {
                    b1 <- ((1-predictions_vector[i]) - 0)^2
                }
                if(real_scores[index2]==1)
                {
                    b2 <-   (predictions_vector[index2] -1)^2
                }
                else
                {
                    b2 <-  (predictions_vector[index2] -0)^2
                }
                b3 <-  (((1-win_vector[i])*predictions_vector[index2])-0)^2
            }
            }
        
        
        
        if(scenario==3)
        {
            cat(predictions_vector[i], win_vector2[i], (1-predictions_vector[index2]),"")
            if( real_scores[i]==1 &&real_scores[index2]==0)
                {
                b1 <- (predictions_vector[i] - 1)^2
                b2 <- ((1-predictions_vector[index2]) -1)^2
                b3 <- ((win_vector2[i]*(1-predictions_vector[index2]))-1)^2  
            }
            else {
                if(real_scores[i]==1)
                {
                    b1 <- (predictions_vector[i] - 1)^2
                }
                else
                {
                    b1 <- (predictions_vector[i] - 0)^2
                }
                if(real_scores[index2]==1)
                {
                    b2 <-  ((1-predictions_vector[index2]) -1)^2
                }
                else
                {
                    b2 <- ((1-predictions_vector[index2]) -0)^2
                }
                b3 <-  ((win_vector2[i]*(1-predictions_vector[index2]))-0)^2
            }
            }
        
        
        if(scenario==4)
        {
            cat((1-predictions_vector[i]), (1-win_vector2[i]), (1-predictions_vector[index2]),"")
            if( real_scores[i]==0 &&real_scores[index2]==0)
                {
                b1 <- ((1-predictions_vector[i]) - 1)^2
                b2 <- ((1-predictions_vector[index2]) -1)^2
                b3 <-  (((1-win_vector2[i])*(1-predictions_vector[index2]))-1)^2  

            }
            else {
                if(real_scores[i]==1)
                {
                    b1 <- ((1-predictions_vector[i]) - 1)^2
                }
                else
                {
                    b1 <- ((1-predictions_vector[i]) - 0)^2
                }
                if(real_scores[index2]==1)
                {
                    b2 <- ((1-predictions_vector[index2]) -1)^2
                }
                else
                {
                    b2 <-  ((1-predictions_vector[index2]) -0)^2
                }
                b3 <- (((1-win_vector2[i])*(1-predictions_vector[index2]))-0)^2
            }
            
        
        }
    
    brier_score <- b1 +b2+b3
    brier_score <- brier_score/3 
    total <- brier_score +total
    cat(real_scores[i], real_scores[index2], b1, b2,b3, brier_score, "\n")
    }}}
    return(total)
}
# average <- function()
# {
#     win_vector <- update_prob_edited()
#     brier_vector <- (win_vector-real_scores)^2
#     total <- 0
#     for (i in 1:56) {
#         total <- total + brier_vector[i]
#     }
#     print(paste("The average state wide brier score is", (total/56)))
  
#     all_states <-c("AK",  "AL", "AR", "AZ",  "CA", "CO",  "CT", "DC", "DE", "FL", "GA",  "HI", "IA", "ID",  "IL", "IN", "KS", "KY", "LA", "ME1", "ME2", "MA",  "MD",  "ME", "MI", "MN", "MO", "MS", "MT", "NE1", "NE2", "NE3", "NC", "ND", "NE", "NH",  "NJ", "NM", "NV",  "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT",  "WA", "WI", "WV", "WY")
#     total_2 <- 0
#     for (i in 1:56)
#     {
#         for (j in 1:56)
#         {
#             if(i !=j)
#             {
#                 if(j != 1 &&j !=2 &&j !=3 &&j !=5&&j !=7&&j !=8&&j !=9&&j !=12&&j !=14&&j !=15&&j !=16&&j !=17&&j !=18&&j !=19&&j !=20&&j !=22&&j !=23&&j !=29&&j !=30&&j !=31&&j !=32&&j !=33&&j !=34&&j !=35&&j !=37&&j !=40&&j !=42&&j !=43&&j !=45&&j !=47&&j !=48&&j !=50&&j !=52&&j !=53&&j !=55&&j !=56) 
#                 {
#                     if(i != 1 &&i !=2 &&i !=3 &&i !=5&&i !=7&&i !=8&&i !=9&&i !=12&&i !=14&&i !=15&&i !=16&&i !=17&&i !=18&&i !=19&&i !=20&&i !=22&&i !=23&&i !=29&&i !=30&&i !=31&&i !=32&&i !=33&&i !=34&&i !=35&&i !=37&&i !=40&&i !=42&&i !=43&&i !=45&&i !=47&&i !=48&&i !=50&&i !=52&&i !=53&&i !=55&&i !=56) 
#                     {   
#                         total_2 <-  total_2 + brier_2edited(state1 = all_states[i], state2= all_states[j])
#                     }
#                 }
#             }
#         }     
#     }
#     print(paste("The average brier score with the second way of testing it is", (total_2/1520)))
#     #runs 380 times, each time making 4 scores
# }
# test <-function(){
#     all_states <-c("AK",  "AL", "AR", "AZ",  "CA", "CO",  "CT", "DC", "DE", "FL", "GA",  "HI", "IA", "ID",  "IL", "IN", "KS", "KY", "LA", "ME1", "ME2", "MA",  "MD",  "ME", "MI", "MN", "MO", "MS", "MT", "NE1", "NE2", "NE3", "NC", "ND", "NE", "NH",  "NJ", "NM", "NV",  "NY", "OH", "OK", "OR", "PA",  "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT",  "WA", "WI", "WV", "WY")
#     total_3 <- 0

#     for (i in 1:56)
#     {
#         for (j in 1:56)
#         {
#             if(i !=j)
#             {
#                 if(j != 1 &&j !=2 &&j !=3 &&j !=5&&j !=7&&j !=8&&j !=9&&j !=12&&j !=14&&j !=15&&j !=16&&j !=17&&j !=18&&j !=19&&j !=20&&j !=22&&j !=23&&j !=29&&j !=30&&j !=31&&j !=32&&j !=33&&j !=34&&j !=35&&j !=37&&j !=40&&j !=42&&j !=43&&j !=45&&j !=47&&j !=48&&j !=50&&j !=52&&j !=53&&j !=55&&j !=56) 
#                 {
#                     if(i != 1 &&i !=2 &&i !=3 &&i !=5&&i !=7&&i !=8&&i !=9&&i !=12&&i !=14&&i !=15&&i !=16&&i !=17&&i !=18&&i !=19&&i !=20&&i !=22&&i !=23&&i !=29&&i !=30&&i !=31&&i !=32&&i !=33&&i !=34&&i !=35&&i !=37&&i !=40&&i !=42&&i !=43&&i !=45&&i !=47&&i !=48&&i !=50&&i !=52&&i !=53&&i !=55&&i !=56) 
#                     {
#                         print(i)
#                         total_3 <-  total_3 + brier_3edited(state1 = all_states[i], state2= all_states[j], scenario = 3)
#                     }
#                 }
#             }
#         }
#     }     
    
#     print(paste("The average brier score with the third way of testing it is(CT)", (total_3/380)))
#     #runs 3136 times

# }

a <-function()
{

    sink("./pkremp-2016-brier1.txt", append = T)
    all_states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO",
                    "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK",
                    "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI",
                    "WV", "WY", "ME1", "ME2", "DC", "NE1", "NE2", "NE3")
    real_scores <-c("AK"=0,  "AL"=0, "AR"=0, "AZ"=0,  "CA"=1, "CO"=1,  "CT"=1,  "DE"=1, "FL"=0, "GA"=0,  "HI"=1, "IA"=0, "ID"=0,  "IL"=1, "IN"=0, "KS"=0, "KY"=0, "LA"=0,  "MA"=1,  "MD"=1,  "ME"=1, "MI"=0, "MN"=1, "MO"=0, "MS"=0, "MT"=0, "NC"=0, "ND"=0, "NE"=0, "NH"=1,  "NJ"=1, "NM"=1, "NV"=1,  "NY"=1, "OH"=0, "OK"=0, "OR"=1, "PA"=0,  "RI"=1, "SC"=0, "SD"=0, "TN"=0, "TX"=0, "UT"=0, "VA"=1, "VT"=1,  "WA"=1, "WI"=0, "WV"=0, "WY"=0, "ME1"=1, "ME2"=0,  "DC"=1, "NE1"=0, "NE2"=0, "NE3"=0)

    win_vector <- update_prob_edited()
    brier_vector <- (win_vector-real_scores)^2

    total <- 0
    cat("State: Win Prob: Outcome: Brier Score: \n")
    for (i in 1:56) {
        #if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
        #{
            cat(all_states[i], win_vector[i], real_scores[i], brier_vector[i], "\n")
            total <- total + brier_vector[i]
        #}
    }
    cat("The average state wide brier score is", (total/56))
    sink()
}

b <-function()
{

    sink("./pkremp-2016-brier2.txt", append = T)
    all_states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO",
                    "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK",
                    "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI",
                    "WV", "WY", "ME1", "ME2", "DC", "NE1", "NE2", "NE3")

    total <- 0
    cat("State1: State2: Scenario: Win Prob(State1): Conditional Prob(State1): Win Prob(State2): Outcome1: Outcome2: Brier Score: \n")
    for (i in 1:56) {
        if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
            {
                total <- total + brier_2edited(state2 = all_states[i])
            }   
        }
    
    cat("The average brier score with the second way of testing it is", (total/1848))
    sink()
}

d <-function()
{

    sink("./pkremp-2016-brier3-TT.txt", append = T)
    all_states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA",
                    "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO",
                    "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK",
                    "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI",
                    "WV", "WY", "ME1", "ME2", "DC", "NE1", "NE2", "NE3")

    total <- 0
    cat("State1: State2: Win Prob(State1): Conditional Prob(State1): Win Prob(State2): Outcome1: Outcome2: Brier Score(State1): Brier Score(State2): Brier Score(Conditional): Brier Score(Average)\n")
    for (i in 1:56) {
        if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
            {
                total <-  total + brier_3edited(state2= all_states[i], scenario = 4)
            }   
        }
    
    cat("The average brier score with the third way of testing it is(TT)", (total/462))
    sink()
}

shares <-function()
{
    clinton_states = NULL
    trump_states = NULL 
    clinton_scores_list = NULL 
    target_nsim = 1000
    states <- names(clinton_scores_list)
    lower_clinton <- sapply(clinton_scores_list, function(x) x[1]/100)
    upper_clinton <- sapply(clinton_scores_list, function(x) x[2]/100)
    sim <- draw_samples(clinton_states = clinton_states, trump_states = trump_states, states = states, 
                        upper_clinton = upper_clinton, lower_clinton = lower_clinton, 
                        target_nsim = target_nsim)
    ev_dist <- (sim[["matrix"]] > .5) %*% ev
    state_win <- colMeans(sim[["matrix"]] > .5)
    shares_vector <- colMeans(sim[["matrix"]])
    real_shares <-c("AK"=37.6,  "AL"=34.7, "AR"=33.7, "AZ"=45.5,  "CA"=62.3, "CO"=48.2,  "CT"=54.7,  "DE"=53.4, "FL"=47.8, "GA"=45.9,  "HI"=62.2, "IA"=42.2, "ID"=27.5,  "IL"=56, "IN"=37.9, "KS"=36.3, "KY"=32.7, "LA"=38.4,  "MA"=61,  "MD"=61.3,  "ME"=48, "MI"=47.4, "MN"=46.9, "MO"=38.2, "MS"=40.1, "MT"=35.9, "NC"=46.8, "ND"=27.7, "NE"=34.4, "NH"=47.6,  "NJ"=55.5, "NM"=48.3, "NV"=47.9,  "NY"=59.5, "OH"=43.7, "OK"=28.9, "OR"=52, "PA"=47.9,  "RI"=55.5, "SC"=40.7, "SD"=31.7, "TN"=34.9, "TX"=43.5, "UT"=27.5, "VA"=50.2, "VT"=61.6,  "WA"=54.3, "WI"=47, "WV"=26.5, "WY"=22.5, "ME1"=0, "ME2"=0,  "DC"=92.8, "NE1"=0, "NE2"=0, "NE3"=0)
    real_shares = real_shares *.01
    total <- 0
    for(i in 1:56)
    {
        if(real_shares[i]!=0)
        {
            if(i!=0 && i!=1&& i!=2&& i!=3&& i!=5&& i!=7&& i!=8&& i!=11&& i!=13&& i!=14&& i!=15&& i!=16&& i!=17&& i!=18&& i!=19&& i!=20&& i!=26&& i!=28&& i!=29&& i!=31&& i!=34&& i!=36&& i!=37&& i!=39&& i!=41&& i!=42&& i!=44&& i!=46&& i!=47&& i!=49&& i!=50&& i!=51&& i!=53&& i!=54&& i!=56)
            {
                total <- total + ((shares_vector[i]-real_shares[i])^2)
            }
        }
    }
    print(shares_vector)
    print(real_shares)
    print(paste("The vote share brier score is ", total/20))

}

shares_3 <- function()
{
   clinton_states = NULL
    trump_states = NULL 
    clinton_scores_list = NULL 
    target_nsim = 1000
    states <- names(clinton_scores_list)
    lower_clinton <- sapply(clinton_scores_list, function(x) x[1]/100)
    upper_clinton <- sapply(clinton_scores_list, function(x) x[2]/100)
    sim <- draw_samples(clinton_states = clinton_states, trump_states = trump_states, states = states, 
                        upper_clinton = upper_clinton, lower_clinton = lower_clinton, 
                        target_nsim = target_nsim)
    ev_dist <- (sim[["matrix"]] > .5) %*% ev
    state_win <- colMeans(sim[["matrix"]] > .5)
    shares_vector <- colMeans(sim[["matrix"]])
    real_shares <-c("AK"=37.6,  "AL"=34.7, "AR"=33.7, "AZ"=45.5,  "CA"=62.3, "CO"=48.2,  "CT"=54.7,  "DE"=53.4, "FL"=47.8, "GA"=45.9,  "HI"=62.2, "IA"=42.2, "ID"=27.5,  "IL"=56, "IN"=37.9, "KS"=36.3, "KY"=32.7, "LA"=38.4,  "MA"=61,  "MD"=61.3,  "ME"=48, "MI"=47.4, "MN"=46.9, "MO"=38.2, "MS"=40.1, "MT"=35.9, "NC"=46.8, "ND"=27.7, "NE"=34.4, "NH"=47.6,  "NJ"=55.5, "NM"=48.3, "NV"=47.9,  "NY"=59.5, "OH"=43.7, "OK"=28.9, "OR"=52, "PA"=47.9,  "RI"=55.5, "SC"=40.7, "SD"=31.7, "TN"=34.9, "TX"=43.5, "UT"=27.5, "VA"=50.2, "VT"=61.6,  "WA"=54.3, "WI"=47, "WV"=26.5, "WY"=22.5, "ME1"=0, "ME2"=0,  "DC"=92.8, "NE1"=0, "NE2"=0, "NE3"=0)
    real_shares = real_shares *.01
    total <- 0

    for(i in 1:56)
    {
        for(j in 1:56)
        {
            if(real_shares[i]!=0)
            {
                if(real_shares[j]!=0)
                {
                    if(i!=j)
                    {
                        brier <- ((shares_vector[i]-real_shares[i])^2)
                        brier <- brier + ((shares_vector[j]-real_shares[j])^2)
                        brier <- brier + ((shares_vector[i]^2-real_shares[i]^2)^2)
                        brier <- brier + ((shares_vector[j]^2-real_shares[j]^2)^2)
                        brier <- brier + (((shares_vector[i]*shares_vector[j])-(real_shares[i]*real_shares[j]))^2)
                        brier <- brier/5
                        total <- total + brier
                    }
                }
            }
        }
    }
    print(paste("The vote share brier score third method is ",total/2550))   
}

