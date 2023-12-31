rm(list = ls())
options(mc.cores = parallel::detectCores())

library(rstan)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(curl)
library(shinystan)
library(rmarkdown)


####################
# Useful functions #
####################

corr_matrix <- function(m){
    (diag(m)^-.5 * diag(nrow = nrow(m))) %*% m %*% (diag(m)^-.5 * diag(nrow = nrow(m))) 
}

cov_matrix <- function(n, sigma2, rho){
    m <- matrix(nrow = n, ncol = n)
    m[upper.tri(m)] <- rho
    m[lower.tri(m)] <- rho
    diag(m) <- 1
    (sigma2^.5 * diag(n))  %*% m %*% (sigma2^.5 * diag(n))
}

logit <- function(x) log(x/(1-x))
inv_logit <- function(x) 1/(1 + exp(-x))


start_date <- as.Date("2016-04-01") # Keeping all polls after April 1, 2016.

########################################################
# Downloading poll data from the HuffPost Pollster API #
########################################################

# Creating a vector of URL stubs, to fetch csv files for each state from HuffPost
# Note that FL and CA have different stubs.

state_name <- datasets::state.name
names(state_name) <- datasets::state.abb

# Note that URLS are formed differently for CA and FL
stubs <- state_name[-(which(names(state_name) %in% c("CA", "FL", "DC")))] %>% 
    tolower %>% sub(" ", "-", .) %>% 
    paste("2016", ., "president-trump-vs-clinton", sep = "-")

stubs <- c("2016-general-election-trump-vs-clinton",
           stubs, 
           "2016-california-presidential-general-election-trump-vs-clinton",
           "2016-florida-presidential-general-election-trump-vs-clinton")

names(stubs) <- c("--", 
                  names(state_name[-(which(names(state_name) %in% c("CA", "FL", "DC")))]), 
                  "CA", "FL")

stubs <- stubs[order(names(stubs))]

download_csv <- function(stub){
    url <- paste("http://elections.huffingtonpost.com/pollster/", stub, ".csv", sep = "")
    connection <- curl(url, "r") 
    df <- read.csv(connection, stringsAsFactors = FALSE)
    close(connection)
    message("Downloaded ", url)
    return(df)
}

# Download the data and put everything in a single df
all_polls <- map_df(stubs, download_csv, .id = "state")

colnames(all_polls) <- colnames(all_polls) %>% tolower

system("cp all_polls.csv all_polls-old.csv")
write.csv(all_polls, "all_polls.csv")
# all_polls <- read.csv("all_polls.csv", stringsAsFactors = FALSE, header = TRUE)

############################################################
# Cleaning up downloaded polling data / creating variables #
############################################################

df <- all_polls %>% 
    tbl_df %>%
    rename(pop = number.of.observations,
           vtype = population,
           method = mode) %>%
    mutate(begin = as.Date(start.date),
           end   = as.Date(end.date),
           t = end - (1 + as.numeric(end-begin)) %/% 2,
           entry_date = entry.date.time..et. %>% substr(1,10) %>% as.Date) %>%
    filter(t >= start_date & !is.na(t)
           & (vtype == "Likely Voters" | 
                  vtype == "Registered Voters" | 
                  vtype == "Adults") # This is to get rid of rows showing disaggregated polls by party ID.
           & pop > 1) %>%
    mutate(pollster = str_extract(pollster, pattern = "[A-z0-9 ]+") %>% sub("\\s+$", "", .),
           pollster = replace(pollster, pollster == "Fox News", "FOX"), # Fixing inconsistencies in pollster names
           pollster = replace(pollster, pollster == "WashPost", "Washington Post"),
           pollster = replace(pollster, pollster == "ABC News", "ABC"),
           undecided = ifelse(is.na(undecided), 0, undecided),
           other = ifelse(is.na(other), 0, other) + 
               ifelse(is.na(johnson), 0, johnson) + 
               ifelse(is.na(mcmullin), 0, mcmullin),
           sum = clinton + trump,
           week = floor_date(t - days(2), unit = "week") + days(2), # weeks start on Tuesdays.
           day_of_week = as.integer(t - week), # number of days since beginnning of the week
           # Trick to keep "likely voter" polls when multiple results are reported.
           polltype = as.integer(as.character(recode(vtype, "Likely Voters" = "0", 
                                                     "Registered Voters" = "1",
                                                     "Adults" = "2"))), 
           n_clinton = round(pop * clinton/100),
           # Only looking at trump or clinton voters, leaving 3rd party candidates and undecided voters out for now.
           n_respondents = round(pop*(clinton+trump)/100),
           p_clinton = clinton/(clinton+trump),
           # Numerical indices passed to Stan for states, days, weeks, pollsters
           index_s = as.numeric(as.factor(as.character(state))), 
           # Factors are alphabetically sorted: 1 = --, 2 = AL, 3 = AK, 4 = AZ...
           index_t = 1 + as.numeric(t) - min(as.numeric(t)),
           index_w = as.numeric(as.factor(week)),
           index_p = as.numeric(as.factor(as.character(pollster))))  %>%
    arrange(state, t, polltype, sum) %>% 
    distinct(state, t, pollster, .keep_all = TRUE) %>%
    select(state, t, begin, end, entry_date, pollster, polltype, method, 
           p_clinton, n_respondents, n_clinton, week, day_of_week, starts_with("index_")) 

print("New polls today:")
df %>% filter(entry_date >= Sys.Date()-1) %>% select(-polltype, -t, -method, -n_clinton)
cat(nrow(df), "rows in data frame\n")

##############################
# Removing overlapping polls #
##############################

print(nrow(df))

df$overlap_with_prev <- TRUE

while(any(df$overlap_with_prev == TRUE)){
    # Recursively drop polls if their dates overlap
    # Perhaps not the most efficient way, but it gets the job done, and the data frame is small.
    df <- df %>% 
        arrange(state, pollster, t) %>%
        group_by(state, pollster) %>% 
        mutate(end_prev_poll = lag(end), 
               overlap_with_prev = ifelse(!is.na(end_prev_poll), 
                                          end_prev_poll > begin,
                                          FALSE),
               overlap_with_next = ifelse(!is.na(lead(overlap_with_prev)), 
                                          lead(overlap_with_prev), 
                                          FALSE),
               latest_overlap = overlap_with_prev == TRUE & 
                   overlap_with_next == FALSE,
               # Drop all polls which overlap with next poll
               #                              if   next poll is the most recent (latest) poll 
               #                                   in the series of overlapping polls.
               drop_poll = (overlap_with_next == TRUE & 
                                (lead(latest_overlap) == TRUE | is.na(lead(latest_overlap))))
        ) %>%
        filter(drop_poll == FALSE) %>%
        ungroup
    # The while loop keeps going until no overlap is found.
    print(nrow(df))
}

##################
# Useful vectors #
##################

all_polled_states <- df$state %>% unique %>% sort
election_day <- as.Date("2016-11-08")
ndays <- max(df$t) - min(df$t)
all_t <- min(df$t) + days(0:(ndays))
all_t_until_election <- min(all_t) + days(0:(election_day - min(all_t)))
week_for_all_t <- floor_date(all_t - days(2), unit="week") + days(2)

all_weeks <- (floor_date(all_t - days(2), unit = "week") + days(2)) %>% unique # weeks start on Tuesdays
all_weeks_until_election <- (floor_date(all_t_until_election - days(2), unit = "week") + days(2)) %>% unique
all_pollsters <- levels(as.factor(df$pollster))



###################################################################
# Reading 2012 election data to (1) set priors on mu_b and alpha, #
#                               (2) get state_weights,            #
#                               (3) get state_names and EV        #
###################################################################

states2012 <- read.csv("2012.csv", 
                       header = TRUE, stringsAsFactors = FALSE) %>% 
    mutate(score = obama_count / (obama_count + romney_count),
           national_score = sum(obama_count)/sum(obama_count + romney_count),
           diff_score = score - national_score,
           share_national_vote = (total_count*(1+adult_pop_growth_2011_15))
           /sum(total_count*(1+adult_pop_growth_2011_15))) %>%
    arrange(state)
rownames(states2012) <- states2012$state

prior_diff_score <- states2012[all_polled_states[-1],]$diff_score
names(prior_diff_score) <- all_polled_states[-1]

state_weights <- c(0, states2012[all_polled_states[-1],]$share_national_vote / sum(states2012[all_polled_states[-1],]$share_national_vote))
names(state_weights) <- c("--", states2012[all_polled_states[-1],]$state)

all_states <- states2012$state
state_name <- states2012$state_name
names(state_name) <- states2012$state

# Electoral votes, by state:

ev_state <- states2012$ev
names(ev_state) <- states2012$state

#######################################################################
# Unique date indices for days in which national polls were conducted #
#######################################################################

# See the Stan model: trying to speed up computation of the likelihood by calculating
# the weighted mean of state latent vote intentions only for days national polls were conducted.

unique_ts <- unique(df$index_t[df$state == "--"])

natpolls_indices <- data.frame(unique_ts = unique_ts,
                               unique_ws = sapply(unique_ts, function(i) unique(df$index_w[df$index_t == i])))


df$index_t_unique <- sapply(1:nrow(df), 
                            function(i) ifelse(df$state[i] == "--", 
                                               which(natpolls_indices$unique_ts == df$index_t[i]), 
                                               0))

###################
# Creating priors #
###################

# Mean of the mu_b_prior
# 0.486 is the predicted Clinton share of the national vote according to the Time for Change model.
mu_b_prior <- logit(0.486 + c("--" = 0, prior_diff_score))

# The model uses national polls to complement state polls when estimating the national term mu_a.
# One problem until early September, was that voters in polled states were different from average voters :
# Several solid red states still hadn't been polled, the weighted average of state polls was slightly more pro-Clinton than national polls.

score_among_polled <- sum(states2012[all_polled_states[-1],]$obama_count)/
    sum(states2012[all_polled_states[-1],]$obama_count + 
            states2012[all_polled_states[-1],]$romney_count)
alpha_prior <- log(states2012$national_score[1]/score_among_polled)

sigma_mu_b_end <-cov_matrix(n = length(mu_b_prior) - 1, sigma2 = 1/20, rho = 0.5)
sigma_walk_b_forecast <- cov_matrix(length(mu_b_prior) - 1, 7*(0.015)^2, 0.75)
# sigma_poll_error <- cov_matrix(length(mu_b_prior) - 1, 0.08^2, .7) # About 0.08 = 2% sd on inv_logit scale; 0.04 = 1%
# sigma_poll_error <- cov_matrix(length(mu_b_prior) - 1, 0.05^2, .7) # = 1.25% sd = 1% mean absolute deviation (see: http://researchdmr.com/ProbabilityTotalError.pdf)

sigma_poll_error <- cov_matrix(length(mu_b_prior) - 1, 0.04^2, .7) # about 1% sd.

##################################################
# Passing the data to Stan and running the model #
##################################################

out <- stan("state and national polls.stan", 
            data = list(N = nrow(df),                  # Number of polls
                        S = max(df$index_s),           # Number of states
                        T = length(all_t_until_election),           # Number of days
                        W = length(all_weeks_until_election),           # Number of weeks
                        P = max(df$index_p),           # Number of pollsters
                        last_poll_T = length(all_t),
                        last_poll_W = length(all_weeks),
                        T_unique = max(df$index_t_unique), 
                        t_unique = df$index_t_unique,
                        unique_ts = natpolls_indices$unique_ts,
                        unique_ws = natpolls_indices$unique_ws,
                        s = df$index_s,
                        t = df$index_t,
                        w = df$index_w,
                        p = df$index_p,
                        n_clinton = df$n_clinton,
                        n_respondents = df$n_respondents,
                        state_weights = state_weights,
                        alpha_prior = alpha_prior,
                        mu_b_prior =  mu_b_prior,
                        sigma_mu_b_end = sigma_mu_b_end,
                        sigma_walk_b_forecast = sigma_walk_b_forecast,
                        sigma_poll_error = sigma_poll_error,
                        week = as.integer(as.factor(week_for_all_t)),
                        day_of_week = as.integer(all_t - week_for_all_t)),
            chains = 4, iter = 2000)

stan_summary <- capture.output(print(out, pars = c("alpha", "sigma_c", 
                                                   "sigma_u_state", "sigma_u_national",
                                                   "sigma_walk_a_past", "sigma_walk_b_past",
                                                   paste("mu_b[33,", as.character(2:max(df$index_s)),"]", 
                                                         sep =""))))
stan_summary

######################
# Extracting results #
######################

p <- rstan::extract(out, pars = "predicted_score")[[1]]
alpha <- rstan::extract(out, pars = "alpha")[[1]]
mu_a <- rstan::extract(out, pars = "mu_a")[[1]]
mu_b <- rstan::extract(out, pars = "mu_b")[[1]]
mu_c_standardized <- rstan::extract(out, pars = "mu_c")[[1]]
sigma_c <- rstan::extract(out, pars = "sigma_c")[[1]]
mu_c <- as.vector(sigma_c)*mu_c_standardized
sigma_walk_b_past <- rstan::extract(out, pars = "sigma_walk_b_past")[[1]]
sigma_walk_a_past <- rstan::extract(out, pars = "sigma_walk_a_past")[[1]]
sigma_u_state <- rstan::extract(out, pars = "sigma_u_state")[[1]]
sigma_u_national <- rstan::extract(out, pars = "sigma_u_national")[[1]]
poll_error <- rstan::extract(out, pars = "poll_error")[[1]]

dates <- sort(c(all_t[all_t <= all_weeks[length(all_weeks)]], 
                unique(setdiff(all_weeks_until_election + days(3), all_weeks + days(3)))))
dates <- c(dates[-length(dates)], election_day)
# dates = all dates until last Tuesday; followed by:       # for mu_a daily + interpolated mu_b weekly components
#         all Fridays until election day; followed by:     # because mu_b weekly components are centered on Fridays 
# (weeks start on Tuesdays; the midpoint between 2 successive Tuesdays is Friday).
#         election day.
dimnames(p) <- list(1:nrow(p), as.character(dates), all_polled_states)
dimnames(mu_a) <- list(1:nrow(mu_a), as.character(all_t))
dimnames(mu_b) <- list(1:dim(mu_b)[1],
                       as.character(all_weeks_until_election + days(3)),
                       all_polled_states)
dimnames(mu_c) <- list(1:nrow(mu_c), all_pollsters)

# hist(apply(p[,as.character(election_day),-1], 1, function(vec) cor(vec, inv_logit(mu_b_prior[-1]))))

pred <- data.frame(t = rep(dates, length(all_polled_states)),
                   state = as.character(rep(all_polled_states, each = length(dates))),
                   p =    apply(p, c(2,3), median) %>% as.vector,
                   p_sd = apply(p, c(2,3), sd) %>% as.vector,
                   high = apply(p, c(2,3), function(x) quantile(x, .95)) %>% as.vector,
                   low =  apply(p, c(2,3), function(x) quantile(x, .05))  %>% as.vector,
                   clinton_win = apply(p, c(2,3), function(x) mean(x > .5))  %>% as.vector)

# graphs.R will need indexes for last poll day. But the value is not computed by Stan
# in the generated quantities block in Stan.
# Solution: recovering values for pred[pred$t == all_t[length(all_t)],] by interpolation

now <- all_t[length(all_t)]
before <- dates[max(which(dates <= now))]
after  <- dates[min(which(dates >= now))]
days_before <- days(now-before)@day
days_after <-days(after-now)@day
w_before <- days_after/(days_before + days_after)
w_after <-  days_before/(days_before + days_after)
w <- data.frame(w = c(w_before, w_after), t = c(before, after))

interpolated_rows <- pred %>% filter(t == before | t == after) %>% 
    left_join(w, by = "t") %>% 
    group_by(state) %>%
    summarize(p = weighted.mean(p, w), 
              p_sd = weighted.mean(p_sd, w), 
              high = weighted.mean(high, w), 
              low = weighted.mean(low, w), 
              clinton_win = weighted.mean(clinton_win, w)) %>% 
    mutate(t = now) 

pred <- bind_rows(pred, interpolated_rows)


# cov_logit_p <- lapply(as.character(all_t), function(t) cov(logit(p[,t,-1])))
# names(cov_logit_p) <- as.character(all_t)

# median_logit_forecast <- logit(pred$p[pred$t == election_day & pred$state != "--"])
# names(median_logit_forecast) <- pred$state[pred$t == election_day & pred$state != "--"]
# cov_logit_forecast <- cov(logit(p[,as.character(election_day),-1]))
# cov_forecast <- cov(p[,as.character(election_day),-1])

sim_forecast <- p[,as.character(election_day),-1]

#################################################
# Predicted electoral votes for each simulation #
#################################################

sim_win <- sim_forecast > 0.5

all_non_polled_states <- setdiff(all_states, all_polled_states[-1])
non_polled_win        <- states2012[all_non_polled_states,]$score > .5
names(non_polled_win) <- all_non_polled_states

non_polled_win_matrix <- rep(non_polled_win, nrow(sim_win)) %>%
    matrix(nr = nrow(sim_win), byrow = TRUE,
           dimnames = list(1:nrow(sim_win), all_non_polled_states))

sim_win_all_states <- cbind(sim_win, non_polled_win_matrix)

result_ev_all_states <- sim_win_all_states %*% ev_state[colnames(sim_win_all_states)]

# P(win electoral college)
mean(result_ev_all_states >= 270)

# P(win national vote)
mean(p[, as.character(election_day), "--"] > .5)

##################################################################################
# Sorting predictions for every state, ready to be passed to ggplot for dotchart #
##################################################################################

pr_polled_states <- pred %>% filter(t == election_day & state != "--") %>%
    arrange(-clinton_win) %>%
    transmute(state, p=100*clinton_win, polled = TRUE)
all_non_polled_states <- setdiff(all_states, all_polled_states[-1])
non_polled_win <- states2012[all_non_polled_states,]$score > .5
names(non_polled_win) <- all_non_polled_states
pr_other_states <- data.frame(state = names(non_polled_win), p = 100*non_polled_win, polled = FALSE)

pr_all_states <- rbind(pr_other_states, pr_polled_states)

pr_all_states$position <- ifelse(pr_all_states$p == 100 & !pr_all_states$polled,
                                 0, ifelse(pr_all_states$p == 0 & !pr_all_states$polled,
                                           2,1))
pr_all_states <- pr_all_states %>% arrange(-p, position)

pr_all_states$state_name <- state_name[as.character(pr_all_states$state)]
pr_all_states$order_states <- nrow(pr_all_states):1
pr_all_states$state_name <- factor(pr_all_states$state_name, levels = pr_all_states$state_name[pr_all_states$order_states])

#############################################################################
# Take a subset of 100 simulations of predicted scores to display in ggplot #
#############################################################################

p_subset <- p[sample(1:nrow(p), 100, replace = FALSE),,]

######################
# Info on last polls #
######################

# Recording date/time of last model update and recording time zone.
time_lastrun <- Sys.time()
attributes(time_lastrun)$tzone <- Sys.timezone()

last_polls <- df %>% 
    arrange(desc(entry_date)) %>% 
    filter(entry_date >= Sys.Date() -1) %>%
    mutate(p_clinton = round(p_clinton*100, 1), 
           p_trump = 100-p_clinton, 
           N=n_respondents) %>% 
    select(entry_date, pollster, state, p_clinton, p_trump, N)

save(list = c("sim_forecast", "ev_state"), file = "last_sim.RData")

rm(out)
rm(p)
save.image("out.RData")
render("report.Rmd")