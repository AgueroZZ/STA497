#Use sim:
#Need to install package "coxed" first!
set.seed(10086)
N = 200
library(coxed)
simdata <- sim.survdata(N=N, T=200, xvars=3, censor=.3, num.data.frames = 1)


#True Betas:
simdata$betas
# simdata$betas
#[1,] -0.1546995
#[2,]  0.1045637
#[3,] -0.1277426


#Coxph Estimates:
model <- coxph(Surv(y, failed) ~ X1+X2+X3, data=simdata$data,ties = "breslow")
model$coefficients
#       coef    exp(coef)  se(coef)  z     p
#  X1 -0.28287   0.75362  0.18062 -1.566 0.117
#  X2  0.22606   1.25365  0.16029  1.410 0.158
#  X3  0.02978   1.03023  0.16186  0.184 0.854


#INLA Estimates:
formula <- inla.surv(y,failed) ~ X1+X2+X3
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),data = simdata$data, family = "coxph",
                   control.hazard = list(model="rw1",n.intervals = 20),control.inla = list(strategy = 'gaussian'),
                   control.fixed(prec = 0.05))
Inlaresult$summary.fixed
#                 mean        sd  0.025quant    0.5quant  0.975quant        mode          kld
# (Intercept) -5.44862936 0.2550881 -5.98322061 -5.43729910 -4.97832627 -5.41532815 2.588519e-04
# X1          -0.31576677 0.1820669 -0.67322705 -0.31577118  0.04139071 -0.31576460 5.036596e-09
# X2           0.25294302 0.1592903 -0.05979761  0.25293852  0.56542298  0.25294295 1.382366e-10
# X3           0.04123724 0.1596495 -0.27220552  0.04123151  0.35442691  0.04123351 7.617497e-09


#Now, let's do it in our proposed method:
#Manually setting up model_data that our algorithm will use 
data <- data_frame(X1=simdata$data$X1,X2=simdata$data$X2,X3=simdata$data$X3, ID = 1:N, censoring = ifelse(simdata$data$failed,1,0), times = simdata$data$y, entry = rep(0,length(simdata$data$y)))
true_etas <- simdata$xb
model_data <- list(n=N,A=NULL,M=0,p=1,Nd = N-1,X=as.matrix(simdata$xdata))
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
#
#
#Let's check how long our "hessian_log_likelihood" function is gonna work for this sample size:
system.time(hessian_log_likelihood(W=rnorm(N),model_data))
#
#No theta in this model:
thetagrid <- list(0) # This is the log(precision)
#
#
#
#OPTIMAZTION:
PARALLEL_EXECUTION <- TRUE

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


#Estimates using our algorithm:
compute_marginal_means(c(N,N+1,N+2),sim1opt,model_data)
#[1] -0.28244408  0.22586592  0.02978612




#Variance of Beta:
compute_marginal_variances(c(N,N+1,N+2),sim1opt,model_data)
#[1,] 0.0325692 0.00000000 0.00000000
#[2,] 0.0000000 0.02566175 0.00000000
#[3,] 0.0000000 0.00000000 0.02616222





#Comparison between ours, INLA and CoxPH:
#Our algorithm:
OUR_MEAN <- c(-0.28244408,0.22586592,0.02978612)
OUR_SD <- c(0.1804694,0.1601929,0.1617474)

#INLA's algorithm:
INLA_MEAN <- c(-0.31576677,0.25294302,0.04123724)
INLA_SD <- c(0.1820669,0.1592903,0.1596495)

#Coxph's algorithm: 
COXPH_MEAN <- c(-0.28287174,0.22606094,0.02978013)
COXPH_SD <- c(0.18062,0.16029,0.16186)



#Compare with STAN:
library(rstanarm)

# Data Manipulation for Poisson Trick:
data_STAN <- data
cut_one <-  max(data_STAN$times) + 1
cut_20 <- c(quantile(data_STAN$times, probs = c(1:19)/20), cut_one) %>% round()
data_STAN_20times <- survival::survSplit(formula = Surv(times, censoring) ~., data = data_STAN, cut = cut_20) %>%
  mutate(interval = factor(entry),
         interval_length = times - entry) %>%
  as_data_frame
data_STAN_20times

# Poisson Trick: Using CoxPH

glm(formula = censoring ~ X1+X2+X3 + offset(log(interval_length)),
    family  = poisson(link = "log"),
    data    = data_STAN_20times) %>% summary

#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -7.246544   0.089413 -81.046   <2e-16 ***
#X1          -0.261169   0.184245  -1.418    0.156    
#X2           0.210611   0.161317   1.306    0.192    
#X3           0.008258   0.166106   0.050    0.960 


# Poisson Trick: Using STAN
fit_STAN <- rstanarm::stan_glm(formula = censoring ~ X1+X2+X3 + offset(log(interval_length)),
                                       data = data_STAN_20times,
                                       family = poisson(link = "log"),
                                       prior_intercept = normal(location = 0, scale = 4.472136),
                                       prior = normal(location = 0, scale = 4.472136), chains = 100, iter = 10000)


fit_STAN$coefficients
#(Intercept)          X1          X2          X3 
#-7.25792890    -0.26103220  0.20946623  0.00739524 

fit_STAN$ses
#(Intercept)          X1          X2          X3 
# 0.08962793  0.18435508  0.16112624  0.16623627 



#Final Comparison:
par(mfrow=c(1,1))

#For X1:
hist(as.data.frame(fit_STAN)$X1,freq = F, main = "Comparison of Algorithms for coefficent of X1", xlab = "X1",
     breaks = seq(-1.8,1,by=0.01)
                  ,col="aliceblue",border = "aliceblue")
lines(dnorm(seq(-1.8,1,by=0.001),mean = OUR_MEAN[1], sd = OUR_SD[1])~seq(-1.8,1,by=0.001),type="l",col="Red")
lines(dnorm(seq(-1.8,1,by=0.001),mean = INLA_MEAN[1], sd =INLA_SD[1])~seq(-1.8,1,by=0.001),type="l",col="Brown")
lines(dnorm(seq(-1.8,1,by=0.001),mean = COXPH_MEAN[1], sd =COXPH_SD[1])~seq(-1.8,1,by=0.001),type="l",col="Blue")
legend(-1.8, 2.2, legend=c("Proposed", "INLA", "Coxph(Assume Normal)"),
       col=c("Red", "Brown","Blue"), lty=1, cex=0.8)

#For X2:
hist(as.data.frame(fit_STAN)$X2,freq = F, main = "Comparison of Algorithms for coefficent of X2", xlab = "X2",
     breaks = seq(-0.6,1,by=0.01)
     ,col="aliceblue",border = "aliceblue")
lines(dnorm(seq(-0.6,1,by=0.001),mean = OUR_MEAN[2], sd = OUR_SD[2])~seq(-0.6,1,by=0.001),type="l",col="Red")
lines(dnorm(seq(-0.6,1,by=0.001),mean = INLA_MEAN[2], sd =INLA_SD[2])~seq(-0.6,1,by=0.001),type="l",col="Brown")
lines(dnorm(seq(-0.6,1,by=0.001),mean = COXPH_MEAN[2], sd =COXPH_SD[2])~seq(-0.6,1,by=0.001),type="l",col="Blue")
legend(-0.6, 2.5, legend=c("Proposed", "INLA", "Coxph(Assume Normal)"),
       col=c("Red", "Brown","Blue"), lty=1, cex=0.8)



#For X3:
hist(as.data.frame(fit_STAN)$X3,freq = F, main = "Comparison of Algorithms for coefficent of X3", xlab = "X3",
     breaks = seq(-0.8,0.8,by=0.01)
     ,col="aliceblue",border = "aliceblue")
lines(dnorm(seq(-0.8,0.8,by=0.001),mean = OUR_MEAN[3], sd = OUR_SD[3])~seq(-0.8,0.8,by=0.001),type="l",col="Red")
lines(dnorm(seq(-0.8,0.8,by=0.001),mean = INLA_MEAN[3], sd =INLA_SD[3])~seq(-0.8,0.8,by=0.001),type="l",col="Brown")
lines(dnorm(seq(-0.8,0.8,by=0.001),mean = COXPH_MEAN[3], sd =COXPH_SD[3])~seq(-0.8,0.8,by=0.001),type="l",col="Blue")
legend(-0.8, 2.5, legend=c("Proposed", "INLA", "Coxph(Assume Normal)"),
       col=c("Red", "Brown","Blue"), lty=1, cex=0.8)
