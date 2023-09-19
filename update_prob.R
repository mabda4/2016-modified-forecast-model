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



brier_2edited <- function(state1 = NULL, state2 = NULL)
{
    index1 <- which(names(real_scores) %in% state1)
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()

    win_vector <- update_prob_edited(state2)
    brier_score  <-0
    if( real_scores[index1]==1 &&real_scores[index2]==1)
     {
        brier_score <- brier_score +((predictions_vector[index2]* win_vector[index1])-1)^2
    }
    else {
        brier_score <- brier_score +((predictions_vector[index2]* win_vector[index1])-0)^2
    }

    if( real_scores[index1]==0 &&real_scores[index2]==1)
     {
        brier_score <- brier_score +((predictions_vector[index2]* (1-win_vector[index1]))-1)^2
    }
    else {
        brier_score <- brier_score +((predictions_vector[index2]* (1-win_vector[index1]))-0)^2
    }
    
    win_vector2 <- update_prob_edited(NULL, state2)
    if( real_scores[index1]==1 &&real_scores[index2]==0)
        {
        brier_score <- brier_score +(((1-predictions_vector[index2])* win_vector2[index1])-1)^2
    }
    else {
        brier_score <- brier_score +(((1-predictions_vector[index2])* win_vector2[index1])-0)^2
    }
    if( real_scores[index1]==0 &&real_scores[index2]==0)
        {
        brier_score <- brier_score +(((1-predictions_vector[index2])* (1-win_vector2[index1]))-1)^2
    }
    else {
        brier_score <- brier_score +(((1-predictions_vector[index2])* (1-win_vector2[index1]))-0)^2
    } 
    return(brier_score)
}

brier_3edited <- function(state1 = NULL, state2 = NULL, scenario= NULL)
{
    index1 <- which(names(real_scores) %in% state1)
    index2 <- which(names(real_scores) %in% state2)
    predictions_vector <- update_prob_edited()
    win_vector <- update_prob_edited(state2)
    win_vector2 <- update_prob_edited(NULL, state2)
    brier_score <-0


        if(scenario==1)
    {
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
    return(brier_score)
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

    sink("./FiveThirtyEight-2020-brier1.txt", append = T)
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

