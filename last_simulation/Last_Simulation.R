source("~/last_simulation/00-load-packages.R")
source("~/last_simulation/1.function for PH Model.R")
cores <- detectCores()


set.seed(123)
tdom <- seq(0, 1000, by=0.001)
haz <- rep(0, length(tdom))
cut <- 50
for (i in 1:cut) {
  low <- as.numeric(quantile(tdom,(i-1)/cut))
  high <- as.numeric(quantile(tdom,(i)/cut))
  if(i %% 2 == 1){
    haz[tdom<=high & tdom > low] <- 1/800
  }
  if(i %% 2 == 0){
    haz[tdom<=high & tdom > low] <- 1/200
  }
}

plot(tdom, haz, type='l', xlab='Time domain', ylab='Hazard')
cumhaz <- cumsum(haz*0.001)
Surv <- exp(-cumhaz)
u <- runif(800)
failtimes <- tdom[colSums(outer(Surv, u, `>`))]
hist(failtimes,breaks = 100)




# generate 1000 random samples:
N = 1000
RW2BINS = 30
POLYNOMIAL_DEGREE = 1
PARALLEL_EXECUTION = T

u <- runif(1000)
x <- seq(from = -10, to = 10, length.out = 1000)
eta <- 20/(1+exp(-x))-10
truefunc <- function(x) 20/(1+exp(-x))-10
tibble(x = c(-10,10)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = truefunc)


failtimes <- c()
for (i in 1:1000) {
  hazz <- haz * exp(eta[i])
  cumhaz <- cumsum(hazz*0.001)
  Surv <- exp(-cumhaz)
  failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
}



data <- data_frame(x=x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes==1000,yes = 0, no=1))
data <- rename(data,exposure = x)
data <- data %>% as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "equal"))
data <- arrange_data(data)
data$ID <- 1:length(u)
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
thetagrid <- as.list(seq(-4,2,by = 0.1)) # This is the log(precision)

# Random effect model specification data
model_data$modelspec <- model_data$A %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)
model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
cat("Finished creating model data!\n")

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

save(sim1opt,file = "~/result/Optim.Rdata")






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


save(margmeans_and_vars,file = "~/result/PosteriorMeanVar.Rdata")
















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

ggsave(filename = "~/result/Proposed.pdf",
       plot = simplot,
       width = 3,
       height = 3.5)











sim1optlogpost <- add_log_posterior_values(sim1opt,model_data)
normconsttheta <- normalize_log_posterior(sim1optlogpost$theta_logposterior,sim1optlogpost$theta)
normconstsigma <- normalize_log_posterior(sim1optlogpost$sigma_logposterior,sim1optlogpost$sigma)
sim1optlogpost$theta_logposterior <- sim1optlogpost$theta_logposterior - normconsttheta
sim1optlogpost$theta_post <- exp(sim1optlogpost$theta_logposterior)
sim1optlogpost$sigma_logposterior <- sim1optlogpost$sigma_logposterior - normconstsigma
sim1optlogpost$sigma_post <- exp(sim1optlogpost$sigma_logposterior)


sum(exp(sim1optlogpost$theta_logposterior) * c(0,diff(sim1optlogpost$theta)))
sum(rev(exp(sim1optlogpost$sigma_logposterior)) * c(0,diff(rev(sim1optlogpost$sigma))))

priorfunc <- function(x) exp(model_data$theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(model_data$theta_logprior(-2*log(x)))

thetapostplot <- plot(sim1optlogpost$theta_post~sim1optlogpost$theta,type="l")
ggsave(filename = "~/result/PosterTheta.pdf",
       plot = thetapostplot,
       width = 3,
       height = 3.5)


sigmapostplot <- plot(sim1optlogpost$sigma_post~sim1optlogpost$sigma,type="l")
ggsave(filename = "~/result/PosterSigma.pdf",
       plot = sigmapostplot,
       width = 3,
       height = 3.5)












#Comparison with INLA:
formula <- inla.surv(times,censoring) ~ -1+exposure + f(exposure_binned,model = 'rw2',constr = T)
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),data = data, family = "coxph",
                   control.hazard = list(model="rw2",n.intervals = 20),
                   num.threads = 4)
fhat <- Inlaresult$summary.random$exposure_binned$mean
f.ub <- Inlaresult$summary.random$exposure_binned$`0.975quant`
f.lb <- Inlaresult$summary.random$exposure_binned$`0.025quant`
plotINLA <- data.frame(exposure = Inlaresult$summary.random$exposure_binned$ID)
fit_poly2 <- function(x){
  xx <- poly(x,degree = POLYNOMIAL_DEGREE,raw = T)
  as.numeric(xx %*% cbind(as.numeric(Inlaresult$summary.fixed[1])))
}
mypoly = fit_poly2(plotINLA$exposure) - fit_poly2(plotINLA$exposure[vv])
meanhere <- fhat-fhat[vv] + mypoly

ggplot(plotINLA, aes(x = exposure)) + 
  geom_line(aes(y = meanhere)) + 
  geom_line(aes(y = truefunc(exposure) - truefunc(exposure[vv])),colour = 'blue',linetype = 'solid') +
  theme_bw(base_size = 20)
ggsave(filename = "~/result/INLAFit.pdf",
       plot = sigmapostplot,
       width = 3,
       height = 3.5)



