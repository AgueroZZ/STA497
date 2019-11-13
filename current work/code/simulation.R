#Use sim:
set.seed(520)
N = 200
library(coxed)
simdata <- sim.survdata(N=N, T=80, xvars=1, censor=.3, num.data.frames = 1)

#True Beta
simdata$betas


#Coxph
model <- coxph(Surv(y, failed) ~ X, data=simdata$data)
model$coefficients


#INLA
formula <- inla.surv(y,failed) ~ X
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),data = simdata$data, family = "coxph",
                   control.hazard = list(model="rw2",n.intervals = 20))
Inlaresult$summary.fixed


#Proposed Method:
data <- data_frame(x = simdata$data$X, ID = 1:length(simdata$data$X), censoring = ifelse(simdata$data$failed,1,0), times = simdata$data$y, entry = rep(0,length(simdata$data$y)))
true_etas <- simdata$xb
model_data <- list(n=N,A=NULL,M=0,p=1,Ne=length(true_etas),Nd = N-1,X=matrix(data = data$x,ncol = 1), Wd = N, Wdf = N)

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
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Xd <- model_data$diffmat %*% model_data$X
tibble(x = c(-15,15)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = function(x) exp(model_data$theta_logprior(x)))
model_data$p <- ncol(model_data$X)
model_data$Ne <- model_data$Nd + model_data$n
model_data$Wd <- model_data$M + model_data$p + model_data$Nd
model_data$Wdf <- model_data$M + model_data$p + model_data$Ne
model_data$times <- simdata$data$y
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID

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



global_path <- "~/Desktop/STA497/current work/code/"













system.time(hessian_log_likelihood(W=rnorm(N),model_data))

thetagrid <- list(0) # This is the log(precision)


PARALLEL_EXECUTION <- TRUE



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






compute_marginal_means(N,sim1opt,model_data)



