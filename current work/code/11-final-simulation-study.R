### SIMULATION STUDY ###

# This script contains the simulation study actually reported in the final 
# paper.

# final-sim-results-20190712-12.Rdata: truefunc = 3*x^2, run time 5302 seconds
# final-sim-results-20190712-11.Rdata: truefunc = 10 * (.5*x^2 - 0.5*x + 1.125 + .05*sin(4*pi*x)) run time 5063 seconds

### SETUP ###
global_path <- "/home/alex1/phd/projects/case-crossover/"

PARALLEL_EXECUTION <- TRUE
RW2BINS <- 50
POLYNOMIAL_DEGREE <- 1

source(stringr::str_c(global_path,"biometrics-paper/code/01-functions-for-case-crossover.R"))
library(parallel)
options(mc.cores = parallel::detectCores())

randomseed <- 36547658
set.seed(randomseed)
N <- 10000 # Number of subjects
control_days <- rep(3,N) # Number of control days for each subject

# Simulate the design matrix. First simulate a list, one element for each subject.
# Then after we compute the linear predictor, we will randomly assign one day to be the case day,
# with probabilities defined by the model.
Xlist <- map(control_days,~cbind(runif(.x + 1,0,1)))

names(Xlist) <- 1:N

# The true exposure risk function
# truefunc <- function(x) exp(sin(5*pi*x))
# truefunc <- function(x) exp(x)
# truefunc <- function(x) 3*x^2
truefunc <- function(x) 10 * (.5*x^2 - 0.5*x + 1.125 + .05*sin(4*pi*x))
# truefunc <- function(x) 10 * x

tibble(x = c(0,1)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = truefunc)

# The linear predictor
true_eta_list <- map(Xlist,~data.frame("X" = .x,"eta" = truefunc(.x) + rnorm(length(.x),sd = exp(-.5*12))))

# xx <- Xlist %>% reduce(c)
# yy <- true_eta_list %>% reduce(c)
# plot(yy ~ xx)


# Response
response_list <- true_eta_list %>%
  map(~mutate(.data = .x,
              prob = exp(eta - matrixStats::logSumExp(eta)),
              response = as.numeric(rmultinom(1,1,prob))))

# Now generate the full dataset
sim1data <- map2(response_list,names(response_list),~mutate(.data = .x,subject = as.numeric(.y))) %>%
  reduce(rbind) %>%
  rename(exposure = X,case = response) %>%
  as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "equal")) %>%
  arrange(subject,case)

cor(sim1data$case,sim1data$eta)
sim1data %>%
  group_by(case) %>%
  summarize(prob = mean(prob))


# Take out subjects whose binned covariates are equal on case and control day- they don't contribute
# to likelihood.
equalsubjects <- sim1data %>%
  group_by(subject) %>%
  summarize(numdays = n(),numvals = n_distinct(exposure_binned)) %>%
  filter(numvals < numdays)

# Standardize and remove subjects
# orig_xmean <- mean(sim1data$exposure_binned)
# orig_xsd <- sd(sim1data$exposure_binned)
orig_xmean <- 0
orig_xsd <- 1
sim1data <- sim1data %>%
  anti_join(equalsubjects,by = "subject") %>%
  mutate_at(c("exposure","exposure_binned"),~(.x - orig_xmean)/orig_xsd)

if (nrow(equalsubjects) > 0) {
  control_days <- control_days[-equalsubjects$subject]
}





# Create the model data
# Build in one linear constraint.
cat("Creating model data...\n")
Alist <- list()
Alist$exposure <- create_alist_element(sim1data$exposure_binned)

model_data <- list(
  A = Alist,
  M = Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1,
  n = length(unique(sim1data$subject)),
  X = sparse.model.matrix(case ~ -1 + poly(exposure,degree = POLYNOMIAL_DEGREE,raw = TRUE),data = sim1data)
)
model_data$p <- ncol(model_data$X)
model_data$control_days <- control_days
names(model_data$control_days) <- sim1data$subject %>% unique() %>% sort()
model_data$Nd <- sum(model_data$control_days)
model_data$Ne <- model_data$Nd + model_data$n
model_data$Wd <- model_data$M + model_data$p + model_data$Nd
model_data$Wdf <- model_data$M + model_data$p + model_data$Ne
# Priors. Prior u is P(sigma > u) = alpha. It's the prior median if alpha = .5
# So set u = log(1.25), 50% chance that a unit increase in exposure yields a 
# 25% increase in risk. 
model_data$theta_logprior <- function(theta,prior_alpha = .75,prior_u = log(20)) {
  # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
  # Implement the PC prior directly.
  # P(sigma > u) = alpha.
  # See inla.doc("pc.prec")
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}

# log(precision) for prior on beta
model_data$beta_logprec <- log(.05)


# Differenced matrices...
model_data$diffmat <- create_diff_matrix(model_data$control_days)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$control_days)
model_data$A$exposure$Ad <- model_data$diffmat %*% model_data$A$exposure$A
model_data$Xd <- model_data$diffmat %*% model_data$X
# Random effect model specification data
model_data$modelspec <- model_data$A %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)

model_data$vectorofcolumnstoremove <- round(RW2BINS/2)

cat("Finished creating model data!\n")



# Fit the model. Pick a grid for sigma based on the prior.
tibble(x = c(-15,15)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = function(x) exp(model_data$theta_logprior(x)))

# thetagrid <- as.list(seq(-10,-4,by = .1)) # This is the log(precision)
thetagrid <- as.list(seq(-6,2,by = .1)) # This is the log(precision)

control1 <- list(
  prec = 1e-06,
  stop.trust.radius = 1e-03,
  report.freq = 1,
  report.level = 4,
  start.trust.radius = 100,
  contract.threshold = .25,
  contract.factor = .1,
  expand.factor = 5,
  preconditioner = 1,
  trust.iter = 2000000,
  cg.tol = 1e-06,
  maxit = 1000
)

tm <- proc.time()
sim1opt <- optimize_all_thetas_parallel(
  theta = thetagrid,
  # theta = -6,
  model_data = model_data,
  startingvals = rep(0,model_data$Wd),
  optcontrol = control1,
  doparallel = PARALLEL_EXECUTION
)
rt <- proc.time() - tm

saveresults <- list(
  sim1opt = sim1opt,
  model_data = model_data,
  simmetadata = list(
    N = N,
    RW2BINS = RW2BINS,
    control_days = control_days,
    truefunc = truefunc,
    sim1data = sim1data,
    simtime = rt,
    randomseed = randomseed
  )
)

save(saveresults,file = paste0(global_path,"final-sim-results-20190712-12.Rdata"))
