source("~/STA497/last_simulation/00-load-packages.R")
source("~/STA497/last_simulation/1.function_for_PH_Model_2.R")
cores <- detectCores()


set.seed(100)
tdom <- seq(0, 2000, by=0.001)
haz <- rep(0, length(tdom))
cut <- 50
for (i in 1:cut) {
  low <- as.numeric(quantile(tdom,(i-1)/cut))
  high <- as.numeric(quantile(tdom,(i)/cut))
  if(i %% 2 == 1){
    haz[tdom<=high & tdom > low] <- sample(size=1,x=c(1/1200,1/400),prob=c(0.5,0.5))
  }
  if(i %% 2 == 0){
    haz[tdom<=high & tdom > low] <- sample(size=1,x=c(1/1000,1/600),prob=c(0.5,0.5))
  }
}




# generate 1200 random samples:
N = 1200
RW2BINS = 100
POLYNOMIAL_DEGREE = 1
PARALLEL_EXECUTION = F

u <- runif(N)
x <- seq(from = -20, to = 20, length.out = N)
eta <- 2/(1+exp(-0.2*x)) + rnorm(length(x),sd = exp(-.5*12))
truefunc <- function(x) (2/(1+exp(-0.2*x)))
tibble(x = c(-20,20)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = truefunc)


failtimes <- c()
for (i in 1:N) {
  hazz <- haz * exp(eta[i])
  cumhaz <- cumsum(hazz*0.001)
  Surv <- exp(-cumhaz)
  failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
}



data <- data_frame(x=x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes>=2000,yes = 0, no=1))
for (i in 1:length(data$censoring)) {
  if(data$censoring[i]==1) data$censoring[i] <- rbinom(n=1,size=1,p=0.9)
}

data <- rename(data,exposure = x)
data <- data %>% as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
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
model_data$theta_logprior <- function(theta,prior_alpha = .5,prior_u = log(10)) {
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
thetagrid <- as.list(seq(4,14,by = 0.5)) # This is the log(precision)

# Random effect model specification data
model_data$modelspec <- model_data$A %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)
model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
cat("Finished creating model data!\n")


control1 <- list(
  prec = 1e-08,
  stop.trust.radius = 1e-10,
  report.freq = 10,
  report.level = 4,
  start.trust.radius = 100,
  contract.threshold = .25,
  contract.factor = .5,
  expand.factor = 3,
  trust.iter = 2000,
  maxit = 3000,
  preconditioner = 0
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






sim1optlogpost <- add_log_posterior_values(sim1opt,model_data)
normconsttheta <- normalize_log_posterior(sim1optlogpost$theta_logposterior,sim1optlogpost$theta)
normconstsigma <- normalize_log_posterior(sim1optlogpost$sigma_logposterior,sim1optlogpost$sigma)
sim1optlogpost$theta_logposterior <- sim1optlogpost$theta_logposterior - normconsttheta
sim1optlogpost$theta_post <- exp(sim1optlogpost$theta_logposterior)
sim1optlogpost$sigma_logposterior <- sim1optlogpost$sigma_logposterior - normconstsigma
sim1optlogpost$sigma_post <- exp(sim1optlogpost$sigma_logposterior)









tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd-model_data$M-model_data$p + 1):(model_data$Wd),
  model_results <- sim1optlogpost,
  model_data <- model_data,
  lincomb = lincomb,
  constrA = constrA
)
rt2 <- proc.time() - tm2

print(rt2)



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

ggsave(filename = "~/STA497/result/Proposed.pdf",
       plot = simplot,
       width = 3,
       height = 3.5)












print("The followings should be ones")
sum(exp(sim1optlogpost$theta_logposterior) * c(0,diff(sim1optlogpost$theta)))
sum(rev(exp(sim1optlogpost$sigma_logposterior)) * c(0,diff(rev(sim1optlogpost$sigma))))

priorfunc <- function(x) exp(model_data$theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(model_data$theta_logprior(-2*log(x)))

thetapostplot <- ggplot(mapping = aes(sim1optlogpost$theta,sim1optlogpost$theta_post)) + geom_line()
ggsave(filename = "~/STA497/result/PosterTheta.pdf",plot = thetapostplot)

sigmapostplot <- ggplot(mapping = aes(sim1optlogpost$sigma,sim1optlogpost$sigma_post)) + geom_line()
ggsave(filename = "~/STA497/result/PosterSigma.pdf", plot = sigmapostplot)












#Comparison with INLA:
formula <- inla.surv(times,censoring) ~ -1+exposure  + f(exposure_binned,model = 'rw2',constr = T)
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),control.inla = list(strategy = 'gaussian',int.strategy = 'grid', correct = FALSE),data = data, family = "coxph",
                   control.hazard = list(model="rw2",n.intervals = 20),
                   num.threads = 4)
fhat <- Inlaresult$summary.random$exposure_binned$mean
f.ub <- Inlaresult$summary.random$exposure_binned$`0.975quant`
f.lb <- Inlaresult$summary.random$exposure_binned$`0.025quant`
plotINLA <- data.frame(exposure = Inlaresult$summary.random$exposure_binned$ID)
fit_poly2 <- function(x){
  xx <- poly(x,degree = POLYNOMIAL_DEGREE,raw = T)
  as.numeric(xx %*% cbind(Inlaresult$summary.fixed[,1]))
}
mypoly = fit_poly2(plotINLA$exposure) - fit_poly2(plotINLA$exposure[vv])
meanhere <- fhat-fhat[vv] + mypoly

inlaplot <- ggplot(plotINLA, aes(x = exposure)) + 
   geom_line(aes(y = meanhere)) + 
    geom_line(aes(y = truefunc(exposure) - truefunc(exposure[vv])),colour = 'blue',linetype = 'solid') +
    theme_bw(base_size = 20)

new_compare <- simplot + geom_line(aes(y = meanhere)) +labs(title = "Linear Term vs Liner Term plus Smooth Term curve",
                                                            subtitle = "Red = Linear Term; Orange = Linear Term + Smooth Term; Blue = True function ; Black = INLA",
                                                            x = "Covariate", y = "eta")

ggsave(filename = "~/STA497/result/INLAplot.pdf", plot = new_compare)

