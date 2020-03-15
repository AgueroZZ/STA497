set.seed(123)

cores <- 4
RW2BINS <- 50
PARALLEL_EXECUTION <- TRUE

data <- kidney

data <- data_frame(times = data$time, entry = rep(0,nrow(data)),censoring = data$status, age = data$age, sex = data$sex, disease = data$disease, id = data$id)
data$disease <- relevel(data$disease,ref = "Other")
data <- arrange_data(data)
data$ID <- 1:nrow(data)

xdata <- data.frame(age = data$age,sex = data$sex,GN = ifelse(data$disease == "GN",1,0),AN = ifelse(data$disease == "AN",1,0),PKD = ifelse(data$disease == "PKD",1,0))


Blist <- list()
Blist$id <- create_blist_element(data$id)


iidnum <- Blist %>% map("B") %>% map(ncol) %>% reduce(sum)





model_data <- list(
  X=as.matrix(xdata),
  B = Blist,
  M = iidnum,
  n = length(unique(data$ID)),
  Nd = length(unique(data$ID))-1
)


model_data$beta_logprec <- log(.05)
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Xd <- model_data$diffmat %*% model_data$X
model_data$B$id$Bd <- model_data$diffmat %*% model_data$B$id$B


model_data$p <- ncol(model_data$X)
model_data$Ne <- model_data$n
model_data$Wd <- model_data$M + model_data$p + model_data$Nd
model_data$Wdf <- model_data$M + model_data$p + model_data$Ne


model_data$times <- data$times
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID


model_data$theta_logprior <- function(theta,prior_alpha = c(.5),prior_u = c(2)) {
  # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
  # Implement the PC prior directly.
  # P(sigma > u) = alpha.
  # See inla.doc("pc.prec")
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
model_data$thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 50)
mvQuad::rescale(model_data$thetagrid,domain = c(-1,8))
thetalist <- split(mvQuad::getNodes(model_data$thetagrid),rep(1:nrow(mvQuad::getNodes(model_data$thetagrid))))


# Random effect model specification data
model_data$modelspec <- model_data$B %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)
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



# Plots of marginal posteriors for hyperparameter
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
  geom_line(aes(y = sigma_post),colour = "black",linetype = "solid",size = 1) +
  labs(x = "",y = "") +
  geom_line(aes(y = priorfuncsigma(sigma)),colour = 'black',linetype = 'dashed',size = 0.5) + 
  theme(text = element_text(size = 8))


#Ploting:
ggsave(filename = "~/Kidney_PosterTheta.pdf",plot = thetapostplot1)
ggsave(filename = "~/Kidney_PosterSigma.pdf", plot = sigmapostplot1)




# Posterior for Latent Field:
tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd-model_data$M - model_data$p + 1):(model_data$Wd),
  model_results <- optresults_withlogpostpre,
  model_data <- model_data,
  lincomb = NULL,
  constrA = NULL
)
rt2 <- proc.time() - tm2
print(rt2)

margmeans_iid <- margmeans_and_vars$mean[1:model_data$M]
margmeans_beta <- margmeans_and_vars$mean[(model_data$M+1):length(margmeans_and_vars$mean)]


frailty <- margmeans_iid
u <- sort(unique(model_data$B$id$u))

table_for_fra_ours <- data_frame(id = u, frailty = frailty)



# Compare with INLA:
data$GN <- xdata$GN
data$AN <- xdata$AN
data$PKD <- xdata$PKD

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(2, 0.5)))


formula <- inla.surv(times,censoring) ~ age + sex + GN + AN + PKD + f(id,model = 'iid',hyper = prior.prec)
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),control.inla = list(strategy = 'gaussian',int.strategy = 'grid'), control.fixed = list(prec = 0.05), data = data, family = "coxph")
Inlaresult$summary.fixed[,1][-1]
Inlaresult$summary.random$id[,1:2]$mean



# Compare with Coxph:
coxphfit <- coxph(Surv(times, censoring) ~ age + sex + GN + AN + PKD + frailty(id,dist = "gauss", sparse = F), data=data,ties = "breslow")
coxphfit$coefficients


table_in_total <- data_frame(id = u, Coxph = as.numeric(coxphfit$coefficients[6:length(coxphfit$coefficients)]),Ours = frailty, INLA = Inlaresult$summary.random$id[,1:2]$mean)
table_for_fixed <- data_frame(fixed_Mean = c("age","sex","GN","AN","PKD"),Ours = margmeans_beta, Coxph = as.numeric(coxphfit$coefficients[1:5]), INLA = Inlaresult$summary.fixed[,1][-1])
table_for_fixedSD <- data_frame(fixed_SD = c("age","sex","GN","AN","PKD"),Ours = sqrt(margmeans_and_vars$variance[(model_data$M+1):length(margmeans_and_vars$mean)]), Coxph = sqrt(diag(coxphfit$var)[1:5]), INLA = Inlaresult$summary.fixed[,2][-1])

library(gridExtra)
pdf("RandomEffects.pdf", height=11, width=8.5)
grid.table(table_in_total)
dev.off()




pdf("FixedEffectsMean.pdf", height=11, width=8.5)
grid.table(table_for_fixed)
dev.off()


pdf("FixedEffectsSD.pdf", height=11, width=8.5)
grid.table(table_for_fixedSD)
dev.off()

