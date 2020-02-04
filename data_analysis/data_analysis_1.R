source("~/STA497/last_simulation/00-load-packages.R")
source("~/STA497/last_simulation/SingleTheta_new.R")
set.seed(123)

cores <- 8
RW2BINS <- 50
POLYNOMIAL_DEGREE <- 1
PARALLEL_EXECUTION <- TRUE
N <- nrow(Leuk)

data <- sample_n(Leuk, N, replace = FALSE, weight = NULL, .env = NULL)


data <- data_frame(tpi=data$tpi,times = data$time, entry = rep(0,nrow(data)),censoring = data$cens)

data <- rename(data, exposure = tpi)
data <- data %>% as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
data <- arrange_data(data)
data$ID <- 1:nrow(data)
Alist <- list()
Alist$exposure <- create_alist_element(data$exposure_binned)
model_data <- list(
  A = Alist,
  M = Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1,
  n = length(unique(data$ID))
)

model_data$theta_logprior <- function(theta,prior_alpha = .5,prior_u = log(10)) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}

model_data$p <- 0
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Nd <- model_data$n - 1
model_data$Ne <- model_data$n
model_data$Wd <- model_data$M  + model_data$Nd
model_data$Wdf <- model_data$M  + model_data$Ne
model_data$times <- data$times
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$A$exposure$Ad <- model_data$diffmat %*% model_data$A$exposure$A

model_data$thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 30)
mvQuad::rescale(model_data$thetagrid,domain = c(-1,16))
thetalist <- split(mvQuad::getNodes(model_data$thetagrid),rep(1:nrow(mvQuad::getNodes(model_data$thetagrid))))

# Random effect model specification data
model_data$modelspec <- model_data$A %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)
model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
cat("Finished creating model data!\n")


#Optimization:
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
  theta = model_data$thetagrid,
  model_data = model_data,
  optcontrol = control1,
  doparallel = PARALLEL_EXECUTION
)
rt <- proc.time() - tm




#Hyperparameter:Posterior
optresults_withlogpostpre <- add_log_posterior_values(sim1opt,model_data = model_data)
optresults_withlogpost <- normalize_optresults_logpost(optresults_withlogpostpre) 

margpost1 <- marginal_hyperparameter_posterior(1,optresults_withlogpost)
priorfunc <- function(x) exp(model_data$theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(model_data$theta_logprior(-2*log(x)))


thetapostplot1 <- margpost1$margpost %>%
  mutate(theta_post = exp(thetalogmargpost)) %>%
  ggplot(aes(x = theta)) +
  theme_classic() +
  geom_line(aes(y = theta_post),colour = "red",size = 1) +
  geom_line(aes(y = priorfunc(theta)),colour = "black",size = 0.5) +
  labs(x = "",y = "") +
  theme(text = element_text(size = 8))

sigmapostplot1 <- margpost1$margpost %>%
  mutate(sigma_post = exp(sigmalogmargpost)) %>%
  ggplot(aes(x = sigma)) +
  theme_classic() +
  geom_line(aes(y = sigma_post),colour = "red",size = 1) +
  geom_line(aes(y = priorfuncsigma(sigma)),colour = "black",size = 0.5) +
  labs(x = "",y = "") +
  theme(text = element_text(size = 8))


#Ploting:
ggsave(filename = "~/STA497/Leuk_PosterTheta.pdf",plot = thetapostplot1)
ggsave(filename = "~/STA497/Leuk_PosterSigma.pdf", plot = sigmapostplot1)




#Final Comparison:
#Our Latent Field Computation:
tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd-model_data$M + 1):(model_data$Wd),
  model_results <- optresults_withlogpostpre,
  model_data <- model_data,
  lincomb = NULL,
  constrA = NULL
)
rt2 <- proc.time() - tm2
print(rt2)

margmeans_rw2 <- margmeans_and_vars$mean[1:(RW2BINS-1)]
margmeans_iid <- margmeans_and_vars$mean[RW2BINS:model_data$M]
margvars_rw2 <- margmeans_and_vars$variance[1:(RW2BINS-1)]
margsd <- sqrt(margvars_rw2)
vv <- model_data$vectorofcolumnstoremove

margmeanall <- c(
  margmeans_rw2[1:(vv-1)],
  0,
  margmeans_rw2[vv:length(margmeans_rw2)]
)
margsd <- c(
  margsd[1:(vv-1)],
  0,
  margsd[vv:length(margsd)]
)

#INLA:
cnsA1 <- matrix(rep(0,RW2BINS),nrow = 1)
cnsA1[model_data$vectorofcolumnstoremove] <- 1
conse <- matrix(0, nrow = 1, ncol = 1)


formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse))
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),control.inla = list(strategy = 'gaussian',int.strategy = 'grid', correct = FALSE),data = data, family = "coxph")
fhat <- Inlaresult$summary.random$exposure_binned$mean
fhat[model_data$vectorofcolumnstoremove] = 0
meanhere <- fhat
plotINLA <- data.frame(exposure = Inlaresult$summary.random$exposure_binned$ID)



#GAM:
data1 <- data
x <- sort(unique(model_data$A$exposure$u))
b <- gam(times~s(exposure,pc = plotINLA$exposure[vv]),
         family=cox.ph(),data=data1,weights=censoring)
y <- as.numeric(predict.gam(b,data_frame(exposure = x)))-as.numeric(predict.gam(b,data_frame(exposure=plotINLA$exposure[vv])))


#Plot:
PLOT_TEXT_SIZE = 8
simplot <- tibble(
  x = sort(unique(model_data$A$exposure$u)),
  mymean = margmeanall,
  mymeanlower = mymean - 2*margsd,
  mymeanupper = mymean + 2*margsd
) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  geom_ribbon(aes(ymin = mymeanlower,ymax = mymeanupper),fill = "orange",alpha = .1) +
  geom_line(aes(y = mymean),colour = 'orange',linetype = 'solid') + 
  labs(title = "Comparison of Performance",
       subtitle = "Orange = Proposed;",
       x = "tpi", y = "eta") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))


new_compare <- simplot + geom_line(aes(y = meanhere),colour = "blue") +labs(title = "Comparison of Performance",
                                                                            subtitle = "Orange = Proposed; Blue = INLA ; Green = GAM",
                                                                            x = "Covariate", y = "eta")
new_compare <- new_compare + geom_line(aes(y = y),colour = "green")

ggsave(filename = "~/STA497/leuk_FinalPlot.pdf", plot = new_compare)

