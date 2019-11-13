### ANALYZE OUTPUT OF FINAL SIMULATION STUDIES ###

global_path <- "/home/alex1/phd/projects/case-crossover/"

PLOT_TEXT_SIZE = 12

# final-sim-results-20190712-12.Rdata: truefunc = x^2, run time 5302 seconds
# final-sim-results-20190712-11.Rdata: truefunc = 10 * (.5*x^2 - 0.5*x + .05*sin(4*pi*x)) run time 675 seconds

whichmodel <- "final-sim-results-20190712-12.Rdata"
plotfilename <- "final-plot-truefunc1.Rdata"
thetaplotfilename <- "final-plot-truefunc1-thetapost.Rdata"
RW2BINS <- 50
POLYNOMIAL_DEGREE <- 1

orig_xmean <- 0
orig_xsd <- 1

# truefunc <- function(x) 3*x^2
truefunc <- function(x) 10 * (.5*x^2 - 0.5*x + 1.125 + .05*sin(4*pi*x))

source(stringr::str_c(global_path,"biometrics-paper/code/01-functions-for-case-crossover.R"))
library(parallel)
options(mc.cores = parallel::detectCores())



e1 <- new.env()
load(file = paste0(global_path,whichmodel),envir = e1)
model_data <- e1$saveresults$model_data
sim1opt <- e1$saveresults$sim1opt
rm(e1)

# Diagostics
sim1opt$optimizer %>%
  map("status") %>%
  reduce(c) %>%
  tibble(status = .) %>%
  group_by(status) %>%
  summarize(cnt = n())
# # Plot the posterior for theta.
# sim1optlogpost <- add_log_posterior_values(sim1opt,model_data)
# # # Normalize for plotting.
# normconst <- normalize_log_posterior(sim1optlogpost$theta_logposterior,sim1optlogpost$theta)
# sim1optlogpost$theta_logposterior <- sim1optlogpost$theta_logposterior - normconst
# sum(exp(sim1optlogpost$theta_logposterior) * c(0,diff(sim1optlogpost$theta))) # Should be 1

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

thetapostplot <- sim1optlogpost %>%
  mutate(theta_post = exp(theta_logposterior),
         prior = exp(model_data$theta_logprior(theta))) %>%
  ggplot(aes(x = theta)) +
  theme_classic() +
  geom_line(aes(y = theta_post),colour = "darkorange") +
  stat_function(fun = priorfunc,colour = "red",linetype = "dashed") +
  geom_point(aes(y = theta_post),colour = "black",fill = "orange",pch = 21,alpha = .4) +
  scale_x_continuous(breaks = -8:-4) +
  coord_cartesian(xlim = c(-8,-4)) +
  labs(title = "Posterior distribution for log(precision)",
       subtitle = "Dashed line: prior",
       x = "Theta = log(precision)",
       y = "Density") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

thetapostplot

sigmapostplot <- sim1optlogpost %>%
  mutate(prior = priorfuncsigma(sigma)) %>%
  ggplot(aes(x = sigma)) +
  theme_classic() +
  geom_line(aes(y = sigma_post),colour = "black") +
  stat_function(fun = priorfuncsigma,colour = "black",linetype = "dashed") +
  # geom_point(aes(y = sigma_post),colour = "black",fill = "orange",pch = 21,alpha = .4) +
  # scale_x_continuous(breaks = seq(0,20,by=2)) +
  scale_x_continuous(breaks = seq(10,50,by=10)) +
  # coord_cartesian(xlim = c(0,20)) +
  coord_cartesian(xlim = c(8,50)) +
  labs(y = "Density",x = "") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

sigmapostplot

# # Smooth on the log scale then transform
# sigmalogpostsmoothed <- loess(sigma_logposterior ~ sigma,data = sim1optlogpost)
# sim1optlogpost$sigmalogpostsmoothed <- sigmalogpostsmoothed$fitted
# 
# logpostsigmaplotsmoothed <- sim1optlogpost %>%
#   mutate(prior = priorfuncsigma(sigma)) %>%
#   ggplot(aes(x = sigma)) +
#   theme_classic() +
#   coord_cartesian(xlim = c(0.03,.2)) +
#   geom_line(aes(y = exp(sigmalogpostsmoothed)),colour = "darkorange") +
#   stat_function(fun = priorfuncsigma,colour = "red",linetype = "dashed") +
#   labs(title = "Posterior for sigma (orange) and prior (red)",
#        subtitle = paste0(CAUSE_OF_DEATH," deaths, ",REGION," region(s), ",RW2_BINS," bins"),
#        x = "Standard Deviation (sigma)",
#        y = "Density") +
#   scale_x_continuous(breaks = seq(0,1,by=.02)) +
#   theme(text = element_text(size = PLOT_TEXT_SIZE))
# 
# logpostsigmaplotsmoothed


# Constraints, or "reparametrizations"
# Set the nodes on either side of the centre node equal to zero
node1 <- model_data$vectorofcolumnstoremove - 1
# node2 <- model_data$vectorofcolumnstoremove
# node3 <- model_data$vectorofcolumnstoremove + 1
#
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
# constrA <- NULL
# Linear combinations. Need one per unique value of the covariate
# Function to create linear combination vector for a single value of the covariate
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
# Create a matrix of linear combinations corresponding to the elements of u
# The middle element, where U = 0, should still have a variance, due to beta
uu <- sort(unique(model_data$A$exposure$u))
ii <- c(
  1:(model_data$vectorofcolumnstoremove - 1),
  0,
  (model_data$vectorofcolumnstoremove:(RW2BINS-1))
)


lincomb <- purrr::map2(uu,ii,
                       ~create_single_lincomb(.x,.y,POLYNOMIAL_DEGREE)) %>%
  purrr::reduce(cbind)

# lincomb <- NULL


# NEW mean/variance function
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd - model_data$M - model_data$p + 1):(model_data$Wd), # The first M are U; the last p are beta
  model_results = sim1optlogpost,
  model_data = model_data,
  constrA = constrA,
  lincomb = lincomb
)

margmeans <- margmeans_and_vars$mean[1:(model_data$M)]
margbetas <- margmeans_and_vars$mean[(model_data$M + 1):(model_data$M + model_data$p)]
margvars <- margmeans_and_vars$variance[1:(model_data$M)]
margbetavars <- margmeans_and_vars$variance[(model_data$M + 1):(model_data$M + model_data$p)]

marglincombvars <- margmeans_and_vars$lincombvars

# Add back in the one that was made to be zero
vv <- model_data$vectorofcolumnstoremove
margmeanall <- c(
  margmeans[1:(vv-1)],
  0,
  margmeans[vv:length(margmeans)]
)
margvarall <- c(
  margvars[1:(vv-1)],
  0,
  margvars[vv:length(margvars)]
)

margsd <- sqrt(marglincombvars)



fit_poly <- function(x) {
  xx <- poly(x,degree = POLYNOMIAL_DEGREE,raw = TRUE)

  as.numeric(xx %*% cbind(margbetas))
}


# # Plot the results...
plotfunc <- function(x) exp(truefunc(x) - truefunc(0))
simplot <- tibble(
  x = sort(unique(model_data$A$exposure$u)),
  mypoly = fit_poly(x) - fit_poly(x[vv]) + truefunc(x[vv]) - truefunc(0),
  mymean = mypoly + margmeanall,
  mymeanlower = mymean - 2*margsd,
  mymeanupper = mymean + 2*margsd
) %>%
  mutate_at(c("mypoly","mymean","mymeanlower","mymeanupper"),exp) %>%
  mutate(x = x * orig_xsd + orig_xmean) %>%
  ggplot(aes(x = x)) +
  theme_classic() +
  geom_ribbon(aes(ymin = mymeanlower,ymax = mymeanupper),fill = "lightgrey",alpha = .2) +
  geom_line(aes(y = mymean),colour = 'black') +
  geom_line(aes(y = mymeanupper),colour = "black",linetype = "dashed") +
  geom_line(aes(y = mymeanlower),colour = "black",linetype = "dashed") +
  stat_function(fun = plotfunc,linetype = "dotdash") +
  labs(x = "", y = "exp(eta)") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

simplot

# Save workspace for future changes
# save.image(file = paste0(global_path,"biometrics-paper/data-analysis/sim-workspace-1.Rdata"))

load(file = paste0(global_path,"biometrics-paper/data-analysis/sim-workspace-2.Rdata"))

# Save plots. Save them as R objects then load them together in a separate
# script, combine in one cowplot, and save to disk.
# saveplt <- cowplot::plot_grid(
#                 sigmapostplot + theme(axis.title.x = element_blank()),
#                 simplot + theme(axis.title.x = element_blank()),
#                 
#                 # sigmapostplot + theme(plot.title = element_blank(),plot.subtitle = element_blank()),
#                 # simplot + theme(plot.title = element_blank(),plot.subtitle = element_blank()),
#                 nrow = 1
#               )
# save(saveplt,
#   file = paste0(global_path,"biometrics-paper/data-analysis/figures/",plotfilename)
# )
# # Also save the theta posterior:
# save(thetapostplot,file = paste0(global_path,"biometrics-paper/data-analysis/figures/",thetaplotfilename))

# UPDATE: save plots separately as pdfs. Combine in latex.
# No longer need the separate file where we combine them.

# Sigma post.
ggsave(filename = paste0(global_path,"biometrics-paper/data-analysis/figures/sim-2-sigmapostplot.pdf"),
       plot = sigmapostplot,
       width = 3,
       height = 3.5)

# Simulation results
ggsave(filename = paste0(global_path,"biometrics-paper/data-analysis/figures/sim-2-resultsplot.pdf"),
       plot = simplot,
       width = 3,
       height = 3.5)



