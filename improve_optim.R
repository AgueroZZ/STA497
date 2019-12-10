optimize_latentfield_trustoptim <- function(theta,model_data,startingvals=NULL,random_start_sd = .1,Q = NULL,report_freq = 1,report_level = 2,trcontrol = NULL) {
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  if (length(startingvals) != model_data$Wd) stop(stringr::str_c("Length of starting values: ",
                                                                 length(startingvals),
                                                                 " but length of latent field: ",
                                                                 model_data$Wd))
  # Set up the functions.
  optfunc <- function(W) -1 * log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) -1 * grad_log_posterior_W(W,theta,model_data,Q)
  
  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  # If desired, estimate hessian using sparse FD
  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  if (is.null(trcontrol)) {
    trcontrol <- list(
      prec = 1e-06,
      report.freq = 1,
      report.level = 4,
      start.trust.radius = 100,
      contract.threshold = .25,
      contract.factor = .25,
      expand.factor = 5,
      trust.iter = 2000000,
      cg.tol = 1e-06,
      maxit = 1000,
      preconditioner = 0,
      stop.trust.radius = 1e-08
    )
  }
  opt <- trust.optim(
    x = startingvals,
    fn = optfunc,
    gr = optfuncgrad,
    method = "SR1",
    control = trcontrol
  )
  opt$hessian <- optfunchess(opt$solution)
  list(
    optimizer = list(
      name = "trust_optim",
      hessian = opt$hessian,
      trust_radius = opt$trust.radius,
      norm_grad = sqrt(sum(opt$gradient^2)),
      scaled_norm_grad = sqrt(sum(opt$gradient^2)) / sqrt(model_data$Wd),
      status = opt$status,
      nnz = opt$nnz
    ),
    theta = theta,
    starting = startingvals,
    solution = opt$solution,
    function_value = opt$fval,
    iterations = opt$iterations
  )
}
