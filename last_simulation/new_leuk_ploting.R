source("~/STA497/last_simulation/00-load-packages.R")
source("~/STA497/last_simulation/SingleTheta_new.R")

set.seed(123)

options(mc.cores = 10L)
cores <- 10
RW2BINS <- 50
POLYNOMIAL_DEGREE <- 1
PARALLEL_EXECUTION <- TRUE
N <- nrow(Leuk)

data <- sample_n(Leuk, N, replace = FALSE, weight = NULL, .env = NULL)


data <- data_frame(tpi=data$tpi,times = data$time, entry = rep(0,nrow(data)),censoring = data$cens, wbc = data$wbc, age = data$age, sex = data$sex)

data <- data %>% as_tibble() %>%
  mutate(tpi_binned = bin_covariate(tpi,bins = RW2BINS,type = "quantile"))
data <- arrange_data(data)
data$ID <- 1:nrow(data)
Alist <- list()
Alist$tpi <- create_alist_element(data$tpi_binned)
model_data <- list(
  A = Alist,
  M = Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1,
  n = length(unique(data$ID)),
  X = as.matrix(data_frame(sex = data$sex, age = data$age, wbc = data$wbc))
)



model_data$theta_logprior <- function(theta,prior_alpha = .05,prior_u = 0.3) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
model_data$beta_logprec <- log(.05)


model_data$p <- 3
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Nd <- model_data$n - 1
model_data$Ne <- model_data$n
model_data$Wd <- model_data$M  + model_data$Nd + model_data$p
model_data$Wdf <- model_data$M  + model_data$Ne + model_data$p
model_data$times <- data$times
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$A$tpi$Ad <- model_data$diffmat %*% model_data$A$tpi$A
model_data$Xd <- model_data$diffmat %*% model_data$X

model_data$thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 50)
mvQuad::rescale(model_data$thetagrid,domain = c(-2,18))
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
  geom_line(aes(y = theta_post),colour = "black",size = 1) +
  geom_line(aes(y = priorfunc(theta)),colour = "black",linetype = "dashed",size = 0.5) +
  # coord_cartesian(xlim = c(0,20)) +
  labs(y = "Density",x = "") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

sigmapostplot1 <- margpost1$margpost %>%
  mutate(sigma_post = exp(sigmalogmargpost)) %>%
  ggplot(aes(x = sigma)) +
  theme_classic() +
  geom_line(aes(y = sigma_post),colour = "black",size = 1) +
  geom_line(aes(y = priorfuncsigma(sigma)),colour = "black",linetype = "dashed",size = 0.5) +
  labs(x = "",y = "Density") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))



#Ploting:
ggsave(filename = "~/STA497/Leuk_PosterTheta.pdf",plot = thetapostplot1)
ggsave(filename = "~/STA497/Leuk_PosterSigma.pdf", plot = sigmapostplot1)



#Final Comparison:
#Our Latent Field Computation:
tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Nd + 1):(model_data$Wd),
  model_results <- optresults_withlogpostpre,
  model_data <- model_data,
  lincomb = NULL,
  constrA = NULL
)
rt2 <- proc.time() - tm2
print(rt2)

margmeans_rw2 <- margmeans_and_vars$mean[1:(RW2BINS-1)]
margmeans_linear <- margmeans_and_vars$mean[RW2BINS:length(margmeans_and_vars$mean)]
margsd_linear <- sqrt(margmeans_and_vars$variance[RW2BINS:length(margmeans_and_vars$mean)])
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
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(0.3, 0.05)))

formula <- inla.surv(times,censoring) ~ sex + age + wbc + f(tpi_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse), hyper = prior.prec)
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),control.inla = list(strategy = 'gaussian',int.strategy = 'grid', correct = FALSE), control.fixed = list(prec = 0.05), data = data, family = "coxph")
fhat <- Inlaresult$summary.random$tpi_binned$mean
fhat[model_data$vectorofcolumnstoremove] = 0
meanhere <- fhat
plotINLA <- data.frame(tpi = Inlaresult$summary.random$tpi_binned$ID)
plot(fhat,type = 'l')



#GAM:
x <- sort(unique(model_data$A$tpi$u))
b <- gam(times~sex + age + wbc + s(tpi,pc = plotINLA$tpi[vv]),
         family=cox.ph(),data=data,weights=censoring)
y <- as.numeric(predict.gam(b,data_frame(tpi = x,sex=rep(1,length(x)),age = rep(16,length(x)), wbc = rep(35,length(x)))))-as.numeric(predict.gam(b,data_frame(tpi=plotINLA$tpi[vv],sex = 1, age = 16, wbc = 35)))



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
  geom_ribbon(aes(ymin = mymeanlower,ymax = mymeanupper),fill = "lightgrey",alpha = .5) +
  geom_line(aes(y = mymeanupper),colour = "black",linetype = "blank") +
  geom_line(aes(y = mymeanlower),colour = "black",linetype = "blank") +
  geom_line(aes(y = mymean),colour = 'black',linetype = 'solid') + 
  geom_line(aes(y = truefunc(x) - truefunc(x[vv])),colour = 'black',linetype = 'dotdash') + 
  theme(text = element_text(size = PLOT_TEXT_SIZE))


new_compare <- simplot + geom_line(aes(y = meanhere),colour = "black",linetype = "dashed") 


ggsave(filename = "~/STA497/leuk_FinalPlot.pdf", plot = new_compare)

brinla::bri.basehaz.plot(Inlaresult)