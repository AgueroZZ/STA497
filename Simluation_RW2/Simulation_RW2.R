#Step1: Create the simulated Dataset:
PARALLEL_EXECUTION <- TRUE
RW2BINS <- 100
POLYNOMIAL_DEGREE <- 1
library(parallel)
options(mc.cores = parallel::detectCores())

randomseed <- 12346
set.seed(randomseed)
N <- 280 # Number of subjects

# Simulate the design matrix:
Xlist <- cbind(seq(from = 1, to = 28, length.out = N))
truefunc <- function(x) 1.2 * (.5*x^0.35 - 0.2*log(x+3) + 0.02 + 1.9*sin(0.9*x))


tibble(x = c(1,20)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = truefunc)
# Simulate the true eta:
true_eta_list <- data.frame("X" = Xlist,"eta" = truefunc(Xlist) + rnorm(length(Xlist),sd = exp(-.5*12)))

# Simulate the observed times based on multinomial likelihood:
data <- true_eta_list
data$ID <- 1:nrow(data)
data$times <- rep(0,nrow(data))
for (i in 1:nrow(data)) {
  eta <- data$eta
  prob <- rep(0,nrow(data))
  prob[data$times==0] <- exp(eta[data$times==0] - matrixStats::logSumExp(eta[data$times==0]))
  response = as.numeric(rmultinom(1,1,prob))
  data$times[which(response %in% 1)] <- max(data$times) + rexp(1,3)
}


# Making the data into the form that required by our algorithm:
data <- rename(data,exposure = X)
data <- data %>% as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "equal"))
data$entry <- rep(0,length(data$times))
data$censoring <- rep(1,length(data$times))
data <- arrange_data(data)
data$ID <- 1:N
Alist <- list()
Alist$exposure <- create_alist_element(data$exposure_binned)
model_data <- list(
  A = Alist,
  M = Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1,
  n = length(unique(data$ID)),
  X = sparse.model.matrix(eta ~ -1 + poly(exposure,degree = POLYNOMIAL_DEGREE,raw = TRUE),data = data)
)
model_data$theta_logprior <- function(theta,prior_alpha = .85,prior_u = log(20)) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
model_data$beta_logprec <- log(.5)
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Xd <- model_data$diffmat %*% model_data$X
model_data$p <- ncol(model_data$X)
model_data$Nd <- model_data$n - 1
model_data$Ne <- model_data$n
model_data$Wd <- model_data$M + model_data$p + model_data$Nd
model_data$Wdf <- model_data$M + model_data$p + model_data$Ne
model_data$times <- data$times
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
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
thetagrid <- as.list(seq(-3,2,by = 0.1)) # This is the log(precision)
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
  model_data = model_data,
  startingvals = rep(0,model_data$Wd),
  optcontrol = control1,
  doparallel = PARALLEL_EXECUTION
)
rt <- proc.time() - tm

# current time:
# 227.555   6.055  95.096

# New time:
#  35.586   2.922  15.738



create_single_lincomb <- function(u,idx,degree = POLYNOMIAL_DEGREE) {
  # u: value of the covariate
  # idx: index of random effect U to which u corresponds
  # degree: the degree of polynomial used in the model.
  Uvec <- rep(0,model_data$M)
  if (idx > 0) Uvec[idx] <- 1
  betavec <- u^(1:degree)
  ll <-   c(
    rep(0,model_data$Nd),
    Uvec,
    betavec
  )
  as(ll,"sparseMatrix")
}




uu <- sort(unique(model_data$A$exposure$u))
ii <- c(1:(model_data$vectorofcolumnstoremove-1),0,model_data$vectorofcolumnstoremove : (RW2BINS-1))
lincomb <- map2(uu,ii,~create_single_lincomb(.x,.y,POLYNOMIAL_DEGREE)) %>% reduce(cbind)


node1 <- model_data$vectorofcolumnstoremove - 1

constrA <- Matrix(
  c(
    # Sum to zero
    # rep(0,model_data$Nd),
    # rep(1,model_data$M),
    # rep(0,model_data$p),
    # node1 zero
    rep(0,model_data$Nd),
    c(rep(0,node1 - 1),1,rep(0,model_data$M - node1)),
    rep(0,model_data$p)
    # # node2 zero
    # rep(0,model_data$Nd),
    # c(rep(0,node2 - 1),1,rep(0,model_data$M - node2)),
    # rep(0,model_data$p)
    # # node3 zero
    # rep(0,model_data$Nd),
    # c(rep(0,node3 - 1),1,rep(0,model_data$M - node3)),
    # rep(0,model_data$p)
  ),
  ncol = 1,
  sparse = TRUE
)


tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd-model_data$M-model_data$p + 1):(model_data$Wd),
  model_results <- sim1opt,
  model_data <- model_data,
  lincomb = lincomb,
  constrA = constrA
)
rt2 <- proc.time() - tm2

# current time:
# 71.169   1.793  75.716
# New time:
# 10.083   0.441  10.339 

margmeans <- margmeans_and_vars$mean[1:(model_data$M)]
margbetas <- margmeans_and_vars$mean[(model_data$M+1):(model_data$M+model_data$p)]
margvars <- margmeans_and_vars$variance[1:(model_data$M)]
margbetavars <- margmeans_and_vars$variance[(model_data$M+1):(model_data$M+model_data$p)]
marglincombvars <- margmeans_and_vars$lincombvars

vv <- model_data$vectorofcolumnstoremove
margmeanall <- c(
  margmeans[1:(vv-1)],
  0,
  margmeans[vv:length(margmeans)]
)

margsd <- sqrt(marglincombvars)

fit_poly <- function(x){
  xx <- poly(x,degree = POLYNOMIAL_DEGREE,raw = T)
  as.numeric(xx %*% cbind(margbetas))
}

true_eta_list$u <- model_data$A$exposure$u
graph_true <- true_eta_list
graph_true <- distinct(graph_true,u,.keep_all = T)
graph_true <- arrange(graph_true,u)

PLOT_TEXT_SIZE = 8
simplot <- tibble(
  x = sort(unique(model_data$A$exposure$u)),
  mypoly = fit_poly(x) - fit_poly(x[vv]),
  mymean = mypoly + margmeanall,
  mymeanlower = mymean - 2*margsd,
  mymeanupper = mymean + 2*margsd
) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  geom_ribbon(aes(ymin = mymeanlower,ymax = mymeanupper),fill = "orange",alpha = .1) +
  geom_line(aes(y = mypoly),colour = 'red',linetype = 'dashed') +
  geom_line(aes(y = mymean),colour = 'orange',linetype = 'solid') + 
  geom_line(aes(y = truefunc(x) - truefunc(x[vv])),colour = 'blue',linetype = 'solid') + 
  labs(title = "Linear Term vs Liner Term plus Smooth Term curve",
       subtitle = "Red = Linear Term; Orange = Linear Term + Smooth Term; Blue = True function",
       x = "Covariate", y = "eta") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

simplot







