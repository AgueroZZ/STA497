### Functions to implement INLA for case-crossover model ###
# Contents:



### SET GLOBAL VARIABLES ###
# global_path <- "/home/alex1/phd/projects/case-crossover/"



### LOAD PACKAGES AND HELPER FUNCTIONS ###
source(stringr::str_c(global_path,"biometrics-paper/code/00-load-packages.R"))
source(stringr::str_c(global_path,"diag-of-inv.R"))

### MODEL SETUP ###
# Function to return items needed to fit model
model_setup <- function(formula,data,models) {
  # This function takes a formula and data source and returns a list
  # containing the fixed and random effects design matrices plus 
  # other necessary quantities
  # Arguments:
  #   formula: a formula compatible with lme4::lmer containing fixed and random effects.
  #            The response is a 0/1 indicator indiciating the case day, and the formula
  #            must contain a variable "id" also present in data which groups days from
  #            the same subject together.
  #   data: a data frame, ideally a tibble but definitely inheriting from class
  #         data.frame, containing all variables specified in formula. Data should be
  #         (and will be internally) SORTED by id and date; if the actual date is not
  #         available, a new column containing the sort order will be added.
  #   models: a character vector of the same length as the number of random effects in
  #           formula, containing the random effect structure for each (e.g. iid, rw1, rw2)
  #
  # Returns: a list with the following items:
  #   A: random effects design matrix
  #   X: fixed effects design matrix
  #   M, n, p, Ne, Nd : # of random effects params (ncol(A)); # of subjects (length(unique(data$id)));
  #               # of fixed effects params (ncol(X)); length of linear predictor (sum of all)
  #               subjects' control days + case days); length of differenced linear predictor.
  #   control_days: a named vector containing each subject's id mapped to # of control days
  #   Wd, Wdf: dimensions of the "constrained" latent field and full latent field.
  
  # TODO: implement this function, lol.
}

# Functions to create the differencing matrix
# First create a differencing matrix of supplied dimension
create_single_diff_matrix <- function(d) {
  # d: ROW dimension. Result will be d x (d+1). d is the number of control days in our application
  cbind(Diagonal(d,-1),Matrix(1,d,1))
}
# Now create a function to create the whole big differencing matrix, given a vector of
# number of control days
create_diff_matrix <- function(control_days) {
  purrr::map(control_days,create_single_diff_matrix) %>% bdiag()
}

# dd <- create_diff_matrix(model_data$control_days) # dim = 6,835; size = 220 kB

# Now the inverse of the diff mat transpose crossproduct (dtcp)
# Function to create a single diffmat inverse
create_single_dtcp_inverse <- function(d) {
  Diagonal(d,1) - Matrix(1/(d+1),d,d)
}
# Function to create the whole block diagonal thing
# Will be slow because of bdiag(), but only have to do it once
create_full_dtcp_matrix <- function(control_days) {
  purrr::map(control_days,create_single_dtcp_inverse) %>% bdiag()
}

# Function to take a covariate and return the appropriate element of "Alist"
create_alist_element <- function(u,constraint = NULL) {
  # u: covariate. NOT sorted and MAY contain ties/repeated values, in general.
  # constraint: vector containing values of u for which random effect U should be
  # constrained to be zero. 
  lu <- length(u)
  A <- Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "rw2"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }
  
  list(u = u,A = A,model = model,constrzero = constrzero)
}

# Bin the covariates into 100 bins
# This gets the dimension of the RW part of the latent field down 
# from 3919 + 3919 + 3922 = 11,760 to 100 x 3 = 300 (!)
bin_covariate <- function(u,bins = 100,type = "quantile",custombins = NULL,decimals = 5) {
  if (min(u) < 0) {
    lowerend <- min(u)*1.01
  } else {
    lowerend <- min(u) * .99
  }
  if (type == "quantile") {
    # bin into quantiles. quantile gives the upper end of the range.
    bininfo <- tibble(
      upper = quantile(u,probs = (1:bins)/bins),
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  } else if (type == "equal") {
    # Bin into equally spaced bins
    bininfo <- tibble(
      upper = seq(min(u),max(u),length.out = bins),
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  } else if (type == "custom") {
    if (is.null(custombins)) stop("type == custom but no custom bins provided")
    # The provided bins are treated as upper end points
    # Make sure to provide a big range that covers all the data.
    lowerend <- min(custombins)
    bininfo <- tibble(
      upper = custombins,
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  }
  
  # Now if u[i] is between lower and upper, out[i] = midpoint
  out <- numeric(length(u))
  for (i in 1:length(u)) {
    tmp <- bininfo %>%
      filter(lower < u[i],upper >= u[i]) %>%
      pull(midpoint)
    # if (length(tmp) == 0) cat(i,"\n")
    out[i] <- tmp
  }
  round(out,decimals)
}


### LOG-LIKELIHOOD ###
# This section contains all quantities which directly relate to the log-likelihood

# Function to implement the case-crossover log-likelihood as a function of the
# constrained latent field. So the input is W, where the first Wd elements are Deltas.
# The weird part about this likelihood is that it depends on the data only through
# metadata (the indices). So first create a prep function:
prep_data_for_log_lik <- function(W,model_data) {
  # Take the parameters and model data; return a list of parameters grouped by subject.
  # This is "deltasplit" in the below functions.
  # Get delta
  delta <- W[1:model_data$Nd]
  
  # Split delta into a list containing the values for each control day for each subject
  split(delta,as.numeric(rep(names(model_data$control_days),model_data$control_days)))
}


log_likelihood <- function(W,model_data) {
  # Compute the log-likelihood for data y and parameters W.
  # Arguments:
  #   W: the latent field, a vector of dimension Wd. The first Nd elements are the
  #      "deltas", and are what's actually used to compute the likelihood.
  #   model_data: a list containing the output of model_setup; A and X don't need to be
  #               included (for efficiency). The dimensions of everything are what's needed,
  #               as well as control_days, for grouping when calculating the log-likelihood
  #
  # Returns: a numeric value, the log-likelihood.
  
  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)
  
  # Helper function to compute the log(1 + sum(exp(-delta))) for each person.
  # Use logsumexp to avoid underflow in the individual delta terms inside the sum
  # Then exponentiate the result- if the WHOLE SUM is so close to -Inf that exponentiating
  # it gives nearly zero, it's not a problem, since 1 + x basically equals 1 in that case
  # anyways. Underflow is only a problem when adding the exp(-delta_it) terms.
  compute_loglik_term <- function(deltavec) {
    # deltavec is a vector of deltas for each person
    # Note we are adding both minus signs in here.
    -log(1 + exp(matrixStats::logSumExp(-deltavec)))
  }
  
  # Now apply that function to the whole list and sum the result to get the answer
  purrr::map(deltasplit,compute_loglik_term) %>% reduce(sum)
}

# W <- rnorm(model_data$Wd)
# log_likelihood(W,model_data)
# Timing for Nd = 6,835:
# > microbenchmark::microbenchmark(log_likelihood(W,model_data))
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# log_likelihood(W, model_data) 63.51671 65.55222 68.79937 69.22114 71.00985 78.80851   100

# Function to implement the gradient of the case-crossover log-likelihood
grad_log_likelihood <- function(W,model_data) {
  # Arguments: see log_likelihood
  # Returns: numeric vector of the same length as W representing the gradient.
  #          It will have zeroes on the end of it (see appendix of paper).
  
  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)
  # Helper to compute the gradient term
  compute_gradient_term <- function(deltavec) {
    # Compute the denominator
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    # Return exp(-deltavec)/denom
    exp(-deltavec)/denom
  }
  # The gradient is the concatenation of all these (vector) terms,
  # plus zeroes on the end to make it as long as W
  gradient_front <- purrr::map(deltasplit,compute_gradient_term) %>% reduce(c)
  gradient_back <- rep(0,length(W) - length(gradient_front))
  c(gradient_front,gradient_back)
}

# grad_log_likelihood(W,model_data)
# microbenchmark::microbenchmark(grad_log_likelihood(W,model_data))
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# grad_log_likelihood(W, model_data) 125.3096 130.0573 144.8127 132.9739 136.6671 521.6215   100

# Function to implement the (sparse) NEGATED hessian of the case-crossover log-likelihood
# This is MINUS the second derivative. Because that's what's used in the paper
# So remember when optimizing: provide the NEGATIVE of log_likelihood and gradient_...,
# but provide the POSITIVE of this function.

# Two functions are written. hessian_log_likelihood_x() provides the ELEMENTS of the hessian only
# This is very fast (134 milliseconds on average for the air pollution example)
# hessian_log_likelihood_structure() brute-force computes the hessian using bdiag,
# which is an order of magnitude slower (around a second on average). But the structure
# never changes, so we only need to compute this once.
# hessian_log_likelihood() computes the hessian, but you can supply it the structure and it
# will basically just wrap hessian_log_likelihood_x(). Much faster.

hessian_log_likelihood_structure <- function(W,model_data,triplet = FALSE) {
  # Arguments: see log_likelihood
  # Returns: a list with i and p for the sparse hessian. NOTE: currently does it
  # for dgCMatrix- i.e. not symmetric. I couldn't figure out how to do it for symmetric.
  # Still fast and memory efficient.
  # Returns the structure as (i,p), in column-compressed format. To get (i,j) triplet
  # format do triplet = TRUE
  
  # Split the parameter vector
  deltasplit <- prep_data_for_log_lik(W,model_data)
  
  # Helper function to create the blocks
  # Takes in a vector of deltasplit and returns a block
  # These blocks are dense.
  compute_hessian_block <- function(deltavec) {
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    expdeltavecscaled <- exp(-deltavec)/denom
    
    if (length(expdeltavecscaled) == 1) return(expdeltavecscaled * (1 - expdeltavecscaled))
    
    diag(expdeltavecscaled) - expdeltavecscaled %o% expdeltavecscaled
  }
  
  # Create a list of blocks
  blocklist <- purrr::map(deltasplit,compute_hessian_block)
  
  # Then create a big zero matrix of dimension equal to length(W) - sum(length(deltasplit))
  blocklist <- c(blocklist,Matrix(0,model_data$Wd - model_data$Nd,model_data$Wd - model_data$Nd))
  
  # Final result is a big block diagonal matrix consisting of the created blocks
  # concatenated with the zero block. The Matrix::bdiag function is REALLY slow for
  # the case of many small dense matrices, as we have here (see docs). 
  
  # UPDATE: changed to symmetric structure.
  structure <- bdiag(blocklist)
  # structure <- forceSymmetric(bdiag(blocklist))
  if (triplet) {
    structure <- as(structure,'dgTMatrix')
    return(list(i = structure@i,j = structure@j))
  } else {
    return(list(i = structure@i,p = structure@p))
  }
}

hessian_log_likelihood_x <- function(W,model_data) {
  deltasplit <- prep_data_for_log_lik(W,model_data)
  
  # Helper function to create the blocks
  # Takes in a vector of deltasplit and returns a block
  # These blocks are dense.
  compute_hessian_block <- function(deltavec) {
    denom <- 1 + exp(matrixStats::logSumExp(-deltavec))
    dveclen <- length(deltavec)
    expdeltavecscaled <- exp(-deltavec)/denom
    
    if (length(expdeltavecscaled) == 1) return(expdeltavecscaled * (1 - expdeltavecscaled))
    
    outmat <- diag(expdeltavecscaled) - expdeltavecscaled %o% expdeltavecscaled
    outmat[upper.tri(outmat,diag = TRUE)]
  }
  
  purrr::map(deltasplit,compute_hessian_block) %>% purrr::reduce(c)
  
}

# main_hess <- hessian_log_likelihood_structure(W,model_data)
# hess_x <- hessian_log_likelihood_x(W,model_data)
# 
# all(hess_x == main_hess@x) # TRUE!

hessian_log_likelihood <- function(W,model_data,structure=NULL,hessian_fd = FALSE) {
  # structure is a list containing elements i and p
  # If not provided, will call hessian_log_likelihood_structure()
  if (is.null(structure)) {
    if (hessian_fd) {
      # Need it in triplet format for FD
      structure <- hessian_log_likelihood_structure(W,model_data,triplet = TRUE)
    } else {
      # Else return in column-compressed format
      structure <- hessian_log_likelihood_structure(W,model_data)
    }
  }
  
  if (hessian_fd) {
    # Estimate the hessian by forward differences
    obj <- sparseHessianFD::sparseHessianFD(x = W, fn = log_likelihood, gr = grad_log_likelihood, rows = structure$i+1, cols = structure$j+1, index1 = TRUE,model_data = model_data)
    return(obj$hessian(W))
  } else{
    out <- new("dgCMatrix")
    out@i <- structure$i
    out@p <- structure$p
    out@Dim <- as.integer(rep(length(structure$p) - 1,2))
    out@x <- numeric(length(out@i))
    out <- as(out,'symmetricMatrix')
    out@x <- hessian_log_likelihood_x(W,model_data)
    return(as(out,'dgCMatrix'))
  }
}

# xx1 <- hessian_log_likelihood(W,model_data)
# tmpstruct <- hessian_log_likelihood_structure(W,model_data)
# tmpstruct <- list(i = tmpstruct@i,p = tmpstruct@p)
# xx2 <- hessian_log_likelihood(W,model_data,structure = tmpstruct)
# norm(xx1 - xx2)
# norm(xx2 - main_hess)


# hessian_log_likelihood(W,model_data)
# pryr::object_size(hessian_log_likelihood(W,model_data))
# microbenchmark::microbenchmark(hessian_log_likelihood(W,model_data))
# Unit: milliseconds
# expr                                       min       lq     mean   median       uq     max neval
# hessian_log_likelihood(W, model_data) 683.9824 826.7987 982.7159 1026.382 1095.087 1685.91   100

# DONE: bdiag is really slow, need to write something custom. Takes about a second (!)
# microbenchmark::microbenchmark(hessian_log_likelihood_x(W,model_data))
# Unit: milliseconds
# expr      min       lq     mean   median       uq
# hessian_log_likelihood_x(W, model_data) 119.6059 124.5707 131.2804 131.7892 136.0529
# max neval
# 154.7045   100
#
# Got it down by a factor of almost 10.



### LATENT FIELD ###

# Functions to implement the precision matrix for
# 1) Linear terms only
# 2) Random effect terms only
# 3) Both linear and random effect terms
# Case 3) is the most common but 1) and 2) are needed for the actual examples
# in the paper.

# Linear terms only
Q_matrix_linear <- function(theta,model_data,tau=exp(12),debug = FALSE) {
  # Compute the Q matrix for a model with only linear terms
  # Arguments:
  #   theta: scalar, value of LOG PRECISION log(1/variance) of prior distribution on beta
  #   model_data: usual model data object. Contains the differenced design matrix Xd
  #   tau: INLA's linear predictor noise term
  # Returns:
  #   A large, sparse matrix representing Q from the paper.
  
  # UPDATE: theta is obsolete. will pull this from model data
  if ("beta_logprec" %in% names(model_data)) theta <- model_data$beta_logprec
  
  # Construct the matrix
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(theta))
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% model_data$Xd),
    cbind(-tau*t(model_data$Xd) %*% model_data$lambdainv,Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd))
  )
}

# Q_matrix_linear(1,model_data)
# microbenchmark::microbenchmark(Q_matrix_linear(1,model_data))
# Unit: milliseconds
# expr     min       lq     mean  median       uq     max neval
# Q_matrix_linear(1, model_data) 6.09665 6.180281 6.529313 6.21383 6.266974 16.3048   100
# Surprisingly fast...

# RW2 Models: random effects design matrix and precision matrix
Q_matrix_rw2_one_component <- function(theta,model_data,covariate) {
  # theta: log(precision); log of 1/(rw2 smoothing variance)
  # covariate: character, name of covariate as appears in model_data$A
  # This function creates the component of the Q matrix corresponding to a single covariate.
  # model_data contains an element "A" which is a list; names(A) is the vector of names of
  # covariates to be modelled smoothly. A$covariate is itself a list containing Amat, the 
  # random effects design matrix, and u, the vector containing the UNIQUE values of the covariate. The columns 
  # of A$covariate$Amat should already be ordered according to order(unique(u)); u can thus be either
  # sorted or not. It will be sorted inside this function
  # u should ALREADY be unique and sorted. This just makes this explicit in the code.
  u <- sort(unique(model_data$A[[covariate]]$u))
  
  ul <- length(u)
  du <- diff(u)
  
  H <- Matrix::bandSparse(n = ul,
                          diagonals = list(
                            c(1/du[-(ul-1)],0),
                            c(0,-(1/du[-(ul-1)] + 1/du[-1]),0),
                            c(0,1/du[-1])
                          ),
                          k = c(-1,0,1))
  
  AA <- Matrix::Diagonal(n = ul,x = c(2/du[1],2/(du[-(ul-1)] + du[-1]),2/du[(ul-1)]))
  
  exp(theta) * forceSymmetric(crossprod(H,crossprod(AA,H)))
}

Q_matrix_rw2 <- function(theta,model_data,tau = exp(12)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta[.y],model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)
  
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Ad),
    cbind(-tau*t(Ad)%*% model_data$lambdainv,Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad)))
  )
}

# Compare sparisty structure (should be same...)
# cp1 <- as(crossprod(Ad),"dgTMatrix")
# cp2 <- as(crossprod(Ad,crossprod(model_data$lambdainv,Ad)),"dgTMatrix")
# compare_sparsity <- tibble(i1 = cp1@i,i2 = cp2@i,
#        j1 = cp1@j,j2 = cp2@j,
#        x1 = cp1@x,x2 = cp2@x) %>%
#   mutate(x3 = x1/2)
# 
# compare_sparsity %>% filter(i1 != i2)
# compare_sparsity %>% filter(j1 != j2)
# compare_sparsity %>% filter(x3 != x2) # All equal.

# A <- model_data$A$pm25$A
# D <- model_data$diffmat
# DDtinv <- model_data$lambdainv
# 
# D <- Matrix(c(-1,1,0,0,0,0,-1,1),nrow = 2,byrow = TRUE,sparse = TRUE)
# 
# D2 <- create_single_diff_matrix(2)
# D2 <- bdiag(D2,D2)
# 
# A2 <- Diagonal(6)

# microbenchmark::microbenchmark(Q_matrix_rw2(0,model_data))

# Q matrix for model containing both linear and rw2 terms
Q_matrix_both <- function(theta,model_data,tau = exp(12)) {
  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta[.y],model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)
  
  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv
  
  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }
  
  
  # Linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))
  
  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * crossprod(Ad,crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}

# Make a general function to compute the Q matrix for a model
# This can be written to just call one of the above... somehow. Maybe as an option inside model_data
Q_matrix <- function(theta,model_data,tau = exp(12),forcesymm = TRUE) {
  # theta: vector of hyperparameters. The structure of this will depend on the model,
  # as specified by model_data
  if (is.null(model_data$A)) {
    if (model_data$p == 0) stop("both X and A are null...")
    mat <- Q_matrix_linear(theta,model_data,tau)
  } 
  else {
    if (model_data$p > 0) {
      mat <- Q_matrix_both(theta,model_data,tau)
      # mat <- mat + (1/tau) * Diagonal(dim(mat)[1])   
    } else {
      # warning("random walk-only models are rank deficient for case crossover. adding fudge factor. consider adding a linear term.")
      mat <- Q_matrix_rw2(theta,model_data,tau)
      # mat <- mat + (1/tau) * Diagonal(dim(mat)[1])    
    }
  }
  
  if (forcesymm) return(forceSymmetric(mat))
  mat
}

# microbenchmark::microbenchmark(Q_matrix(.1,model_data))
# Unit: milliseconds
# expr      min       lq     mean   median      uq     max neval
# Q_matrix(0.1, model_data) 6.077121 6.137851 6.541133 6.178669 6.25778 12.9752   100

# (log) prior of W
logprior_W <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  -as.numeric((1/2)*crossprod(W,crossprod(Q,W)))
}

# (log) posterior for a given value of theta
log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  logprior_W(W,theta,model_data,Q) + log_likelihood(W,model_data)
}



# Gradient of log posterior
grad_log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  -as.numeric(crossprod(Q,W)) + grad_log_likelihood(W,model_data)
}

# Hessian of log posterior. Option to pass in structure of log likelihoood hessian.
# Option to pass in Q as well
hessian_log_posterior_W <- function(W,theta = NULL,Q = NULL,model_data,structure = NULL,hessian_fd = FALSE) {
  if (is.null(theta) & is.null(Q)) stop("One of Q or theta must be provided")
  if (is.null(Q)) Q <- Q_matrix(theta,model_data)
  -(Q + hessian_log_likelihood(W,model_data,structure,hessian_fd))
}

# Q <- Q_matrix(.1,model_data)
# microbenchmark::microbenchmark(
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data),
#   hessian_log_posterior_W(W,theta = .1,model_data = model_data)
# )
# Unit: milliseconds
# expr      min       lq     mean median       uq      max neval
# hessian_log_posterior_W(W, Q = Q, model_data = model_data) 863.1047 1033.178 1200.637 1198.622 1317.136 2038.524   100
# hessian_log_posterior_W(W, theta = 0.1, model_data = model_data) 881.3422 1022.531 1196.749 1239.993 1326.682 1895.039   100
# So providing Q or not makes basically no difference, it saves like 4 milliseconds on average.
# But providing the structure...
# tmpstructure <- hessian_log_likelihood_structure(W,model_data)
# tmpstructure <- list(i = tmpstructure@i,p = tmpstructure@p)
# microbenchmark::microbenchmark(
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data),
#   hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = tmpstructure)
# )
# Unit: milliseconds
# expr      min lq      mean    median        uq       max neval
# hessian_log_posterior_W(W, Q = Q, model_data = model_data) 865.1010 1064.3718 1213.0569 1223.1004 1344.2056 2047.0038   100
# hessian_log_posterior_W(W, Q = Q, model_data = model_data, structure = tmpstructure) 120.0503 148.8404  167.4932  166.9949  182.7649  280.0576   100
# Order of magnitude speedup, woohoo!





### HYPERPARAMETERS ###

# Function to compute the log posterior approximation for theta
# model_data will have a function theta_logprior() which takes vector theta
# and evaluates the log of the joint prior
log_posterior_theta <- function(theta,W,model_data,hessian_structure = NULL,Q = NULL) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(W,model_data)
  }
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure)
  
  term1 <- log_likelihood(W,model_data)
  dt <- determinant(Q,logarithm = TRUE) # For this, we DO need the determinant
  # term2_det <- (1/2) * as.numeric(dt$modulus * dt$sign)
  term2_det <- (1/2) * as.numeric(dt$modulus)
  term2 <- logprior_W(W,theta,model_data) # Doesn't contain the determinant
  term3 <- model_data$theta_logprior(theta)
  qcdet <- determinant(Q_p_C,logarithm = TRUE)
  # term4 <- -(1/2)*as.numeric(qcdet$modulus * qcdet$sign) # The gaussian approx evaluated at conditional mode
  term4 <- -(1/2) * as.numeric(qcdet$modulus)
  as.numeric(term1 + term2_det + term2 + term3 + term4)
}

# log_posterior_theta(.1,startingvals,model_data)

# Function to compute the log posterior for sigma, the standard deviation
# sigma = exp(-.5 * theta)
# jacobian is 2/sigma
log_posterior_sigma <- function(sigma,W,model_data,hessian_structure = NULL,Q = NULL) {
  log(2/sigma) + log_posterior_theta(theta = -2 * log(sigma),
                                  W = W,
                                  model_data = model_data,
                                  hessian_structure = hessian_structure,
                                  Q = Q)
}



### OPTIMIZATION ###

# Need functions to take in the parameter and model data and perform the optimization using
# trust.optim and ipoptr. Should return results in a common format.

optimize_latentfield_trustoptim <- function(theta,model_data,startingvals=NULL,random_start_sd = .1,hessian_structure = NULL,Q = NULL,report_freq = 1,report_level = 2,trcontrol = NULL,hessian_fd = FALSE) {
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
  # Get the hessian structure, if not provided
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(startingvals,model_data)
  }
  
  # If desired, estimate hessian using sparse FD
  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure,hessian_fd = hessian_fd)
  
  
  if (is.null(trcontrol)) {
    trcontrol <- list(
      prec = 1e-02,
      stop.trust.radius = sqrt(.Machine$double.eps),
      report.freq = report_freq,
      report.level = report_level,
      start.trust.radius = 1000,
      contract.threshold = .25,
      contract.factor = .5,
      expand.factor = 5,
      preconditioner = 1,
      trust.iter = 1e08,
      cg.tol = 1e-02
    )
  }
  
  opt <- trust.optim(
    x = startingvals,
    fn = optfunc,
    gr = optfuncgrad,
    hs = optfunchess,
    method = "Sparse",
    control = trcontrol
  )
  
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

# control1 <- list(
#   prec = 1e-02,
#   stop.trust.radius = sqrt(.Machine$double.eps),
#   report.freq = 1,
#   report.level = 4,
#   start.trust.radius = 100000,
#   contract.threshold = .25,
#   contract.factor = .5,
#   expand.factor = 5,
#   preconditioner = 1,
#   trust.iter = 1e08,
#   cg.tol = 1e-02,
#   maxit = 1
# )
# control2 <- list(
#   prec = 1e-02,
#   stop.trust.radius = sqrt(.Machine$double.eps),
#   report.freq = 1,
#   report.level = 4,
#   start.trust.radius = 500,
#   contract.threshold = .25,
#   contract.factor = .1,
#   expand.factor = 5,
#   preconditioner = 1,
#   trust.iter = 2000,
#   cg.tol = 1e-02,
#   maxit = 1000
# )


# Wbest <- optimize_latentfield_trustoptim(
#   theta = .001,
#   model_data = model_data,
#   random_start_sd = 0,
#   trcontrol = control1,
#   report_level = 4
# )
# 
# Wbest2 <- optimize_latentfield_trustoptim(
#   theta = 7,
#   model_data = model_data,
#   startingvals = Wbest$solution,
#   trcontrol = control2,
#   report_level = 4
# )

optimize_latentfield_ipopt <- function(theta,model_data,startingvals=NULL,hessian_structure = NULL,Q = NULL,ipcontrol = NULL,boxconstr = c(-Inf,Inf)) {
  optfunc <- function(W) -1 * log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) -1 * grad_log_posterior_W(W,theta,model_data,Q)
  
  
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd)
  if (length(startingvals) != model_data$Wd) stop(stringr::str_c("Length of starting values: ",
                                                                 length(startingvals),
                                                                 " but length of latent field: ",
                                                                 model_data$Wd))
  
  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  # Get the hessian structure if not provided
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(startingvals,model_data)
  }
  
  optfunchess_full <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure)
  ipopt_hess <- function(W,obj_factor,hessian_lambda) {
    obj_factor * forceSymmetric(optfunchess_full(W))@x
  }
  # Required by ipopt...
  optfunchess_struct <- function() {
    # Only need to call this one once...
    hess <- optfunchess_full(startingvals)
    hess <- as(forceSymmetric(hess),'dsTMatrix')
    forStruct = cbind(hess@i+1, hess@j+1)
    tapply(forStruct[,1], forStruct[,2], c)
  }
  ipopt_hess_structure <- optfunchess_struct()
  
  
  if (is.null(ipcontrol)) {
    ipcontrol <- list(
      max_iter = 100,
      print_level = 5
    )
  }
  
  opt <- ipoptr::ipoptr(
    x0 = startingvals,
    eval_f = optfunc,
    eval_grad_f = optfuncgrad,
    eval_h = ipopt_hess,
    eval_h_structure = ipopt_hess_structure,
    lb = rep(boxconstr[1],length(startingvals)),
    ub = rep(boxconstr[2],length(startingvals)),
    opts = ipcontrol
  )
  list(
    optimizer = list(
      name = "ipopt",
      status = opt$status,
      call = opt$call,
      message = opt$message
    ),
    theta = theta,
    starting = startingvals,
    solution = opt$solution,
    function_value = opt$objective,
    iterations = opt$iterations
  )
}
# 
# Wbest <- optimize_latentfield_ipopt(theta = 0,model_data = model_data,
#                            startingvals = rep(0,model_data$Wd),
#                            max_iter = 20)
# 
# Wbest2 <- optimize_latentfield_ipopt(theta = 0.5,model_data = model_data,
#                            startingvals = Wbest$solution,
#                            max_iter = 5)

# Benchmarking for model with RW2 components:
# > microbenchmark::microbenchmark(ipopt_hess(startingvals,1,1))
# Unit: milliseconds
# expr      min       lq     mean   median      uq      max neval
# ipopt_hess(startingvals, 1, 1) 136.0813 139.9142 151.7076 142.2197 144.721 383.3563   100
#
# > microbenchmark::microbenchmark(optfunc(startingvals))
# Unit: milliseconds
# expr      min       lq     mean   median      uq      max neval
# optfunc(startingvals) 66.26799 69.08459 71.67931 71.65213 73.4415 91.55854   100
#
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# optfuncgrad(startingvals) 119.5731 123.8929 133.3798 125.7219 128.1633 363.3237   100
#
# microbenchmark::microbenchmark(
#   grad_log_posterior_W(startingvals,theta,model_data,Q),
#   grad_log_posterior_W(startingvals,theta,model_data,Q = NULL)
# )
# Unit: milliseconds
#                                                   expr      min       lq     mean median       uq      max neval
# grad_log_posterior_W(startingvals, theta, model_data, Q) 121.4813 125.6133 134.6979 127.3386 128.9025 369.0929   100
# grad_log_posterior_W(startingvals, theta, model_data, Q = NULL) 157.9516 162.8538 171.9846 165.0233 167.6070 397.6115   100

# Implement newton's method manually
optimize_latentfield_newton <- function(theta,model_data,startingvals=NULL,random_start_sd = .1,hessian_structure = NULL,Q = NULL,maxiter = 100,eps = 1e-04,debug = FALSE) {
  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  # Get the hessian structure if not provided
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(startingvals,model_data)
  }
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  if (length(startingvals) != model_data$Wd) stop(stringr::str_c("Length of starting values: ",
                                                                 length(startingvals),
                                                                 " but length of latent field: ",
                                                                 model_data$Wd))
  
  optfunc <- function(W) log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) grad_log_posterior_W(W,theta,model_data,Q)
  
  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data,structure = hessian_structure)
  
  # Function to perform a single newton step
  newton_step <- function(W0) {
    # W0: current point.
    gradval <- optfuncgrad(W0)
    hessval <- forceSymmetric(optfunchess(W0))
    increment <- solve(hessval,gradval)
    list(
      normgrad = sqrt(sum(gradval^2)),
      normincr = sqrt(sum(increment^2)),
      solution = increment
    )
  }
  
  # Iterate
  k <- 1
  startingvals <- rnorm(model_data$Wd,sd = 0)
  opt <- newton_step(startingvals)
  while(k <= maxiter & opt$normincr > eps) {
    if (debug) {
      if (k %% 10 == 0 || k == 1) {
        cat("Iteration | fval | norm(grad) | norm(incr)\n")
      }
    }
    cat(k,optfunc(opt$solution),opt$normgrad,opt$normincr,"\n")
    opt <- newton_step(opt$solution)
    k <- k + 1
  }
  
  list(
    optimizer = list(
      name = "newton",
      normgrad = optfuncgrad(opt$solution),
      normincr = opt$normincr
    ),
    theta = theta,
    starting = startingvals,
    solution = as.numeric(opt$solution),
    function_value = optfunc(as.numeric(opt$solution)),
    iterations = k
  )
}

# Create a tibble from a list of list outputs from the opt functions
# These are all in a common format.
optlist_to_tibble <- function(lst) {
  # If lst is a single element of opt, wrap it in a list so the below works
  if ("optimizer" %in% names(lst)) lst <- list(lst)
  tibble(
    theta = purrr::map_dbl(lst,"theta"),
    optimizer = purrr::map(lst,"optimizer"),
    starting = purrr::map(lst,"starting"),
    solution = purrr::map(lst,"solution"),
    function_value = purrr::map_dbl(lst,"function_value"),
    iterations = purrr::map_dbl(lst,"iterations")
  )
}

# Now a function to do the full optimization.
# Starts by doing num_starts searches, in parallel, using ipopt, with random starts, for the first value of theta
# Then uses trustoptim for incremental values of theta
optimize_all_thetas <- function(theta,model_data,num_starts,startingvals=NULL,random_start_sd = .1,hessian_structure = NULL,starting_opt = "trust.optim",inner_opt = "newton",debug = FALSE,doitinparallel = FALSE,...) {
  # theta is a LIST of theta vectors.
  # num_starts is the number of times to run the initial optimization with random starting values
  # Then, the best value from all these runs is taken and used as the starting value for the next value
  # of theta, which is only run once. Then this solution is used for the NEXT theta, and the optimization
  # proceeds in this way in series.
  
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(startingvals,model_data)
  }
  
  
  # Do the first optimization
  if (starting_opt == "ipoptr") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 5
    do_initial_opt <- function(x) optimize_latentfield_ipopt(theta[[1]],
                                                             model_data,
                                                             startingvals,
                                                             hessian_structure,
                                                             print_level = howmuchtoprint)
  }
  else if (starting_opt == "trust.optim") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 1
    do_initial_opt <- function(x) optimize_latentfield_trustoptim(theta[[1]],
                                                                  model_data,
                                                                  startingvals,
                                                                  hessian_structure = hessian_structure,
                                                                  report_freq = howmuchtoprint)
  }
  else if (starting_opt == "newton") {
    warning("I wouldn't recommend using newton as the starting optimizer. Try trust.optim instead")
    do_initial_opt <- function(x) optimize_latentfield_newton(theta[[1]],
                                                              model_data,
                                                              startingvals,
                                                              hessian_structure = hessian_structure,
                                                              debug = debug)
  }
  else {
    stop(stringr::str_c("Unknown initial optimizer: ",optimizer))
  }
  cat("Performing initial optimization...\n")
  
  tm <- proc.time()
  if (num_starts > 1 & doitinparallel) {
    initial_opt <- mclapply(1:num_starts,do_initial_opt) %>% optlist_to_tibble()
  } else {
    initial_opt <- do_initial_opt(1) %>% optlist_to_tibble()
  }
  
  initial_opt_time <- proc.time() - tm
  initial_opt_time <- unname(initial_opt_time["elapsed"])
  
  # Form output table, by taking the row of initial_opt corresponding to the lowest function value
  final_output <- initial_opt %>%
    arrange(function_value) %>%
    slice(1)
  
  final_output <- list(
    list(theta = final_output$theta,
         optimizer = final_output$optimizer,
         starting = final_output$starting[[1]],
         solution = final_output$solution[[1]],
         function_value = final_output$function_value,
         iterations = final_output$iterations
    )
  )
  
  # Now do inner optimization
  # Use the answer to the previous theta as the starting value for the next theta
  K <- length(theta)
  if (inner_opt == "ipoptr") {
    do_inner_opt <- function(theta,startingvals) optimize_latentfield_ipopt(theta,
                                                                            model_data,
                                                                            startingvals,
                                                                            hessian_structure = hessian_structure,
                                                                            print_level = howmuchtoprint
    )
  }
  else if (inner_opt == "trust.optim") {
    do_inner_opt <- function(theta,startingvals) optimize_latentfield_trustoptim(theta,
                                                                                 model_data,
                                                                                 startingvals,
                                                                                 hessian_structure = hessian_structure,
                                                                                 report_freq = howmuchtoprint)
  }
  else if (inner_opt == "newton") {
    do_inner_opt <- function(theta,startingvals) optimize_latentfield_newton(theta,
                                                                             model_data,
                                                                             startingvals,
                                                                             hessian_structure = hessian_structure,
                                                                             maxiter = 1,
                                                                             debug = debug)
  }
  else {
    stop(stringr::str_c("Unknown inner optimizer: ",optimizer))
  }
  cat("Performing inner optimization...\n")
  tm <- proc.time()
  if (!debug) sink("~/tmp3645789078654")
  for (k in 2:K) {
    if (debug) {
      cat("k = ",k,", theta = ",theta[[k]],"\n")
    }
    final_output[[k]] <- do_inner_opt(
      theta = theta[[k]],
      startingvals = final_output[[k-1]]$solution
    )
  }
  if (!debug) {
    sink()
    system("rm ~/tmp3645789078654")
  }
  inner_opt_time <- proc.time() - tm
  inner_opt_time <- unname(inner_opt_time["elapsed"])
  
  cat("Time taken for optimization with ",num_starts," random starts: \nInitial: ",initial_opt_time," seconds.\nInner: ",inner_opt_time,"seconds.\n")
  final_output %>% optlist_to_tibble()
}

# thetagrid <- as.list(seq(0,1,by = .1))
# allthetaopt <- optimize_all_thetas(
#   theta = thetagrid,
#   model_data = model_data,
#   num_starts = 1,
#   random_start_sd = 0,
#   starting_opt = "trust.optim",
#   inner_opt = "trust.optim",
#   debug = TRUE
# )




# Now do a version that just does them all in parallel
optimize_all_thetas_parallel <- function(theta,model_data,startingvals = NULL,random_start_sd = 0,hessian_structure = NULL,optimizer = "trust.optim",optcontrol = NULL,debug = FALSE,doparallel = TRUE) {
  # Trick: to use 0 as a starting value, set random_start_sd = 0. rnorm(1,0,sd = 0) returns 0 wp1
  
  # Startingvals
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  
  if (is.null(hessian_structure)) {
    hessian_structure <- hessian_log_likelihood_structure(rep(0,model_data$Wd),model_data)
  }
  
  if (optimizer == "ipopt") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 5
    do_opt <- function(theta) {
      optimize_latentfield_ipopt(theta,
                                 model_data,
                                 startingvals,
                                 rnorm(model_data$Wd,sd = random_start_sd),
                                 hessian_structure = hessian_structure,
                                 print_level = howmuchtoprint)
    }
  }
  else if (optimizer == "trust.optim") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 1
    if (is.null(optcontrol)) {
      # Set default control parameters for trust.optim
      optcontrol <- list(
        prec = 1e-08,
        stop.trust.radius = sqrt(.Machine$double.eps),
        report.freq = 1,
        report.level = 1,
        start.trust.radius = 1000,
        contract.threshold = .25,
        contract.factor = .5,
        expand.factor = 5,
        preconditioner = 1,
        trust.iter = 200000,
        cg.tol = 1e-06,
        maxit = 1000
      )
    }
    do_opt <- function(theta) {
      cat("Optimizing with theta =",theta,"\n")
      optimize_latentfield_trustoptim(theta = theta,
                                      model_data = model_data,
                                      startingvals = startingvals,
                                      hessian_structure = hessian_structure,
                                      trcontrol = optcontrol)
    }
  }
  
  cat("Performing optimization, start time: ", format(Sys.time(),"%H:%M:%S"),"\n")
  
  tm <- proc.time()
  if (doparallel) {
    opt <- mclapply(theta,do_opt)
  } else {
    opt <- lapply(theta,do_opt)
  }
  opt_time <- proc.time() - tm
  opt_time <- unname(opt_time["elapsed"])
  cat("Time taken for optimization with ",length(theta)," values of theta:",opt_time,"seconds.\n")
  
  # Take the BEST mode for each theta
  opt %>% optlist_to_tibble() %>% group_by(theta) %>% arrange(function_value) %>% slice(1)
  
}

# all_theta_opt <- optimize_all_thetas(theta = as.list(seq(.01,1.01,by = .05)),model_data = model_data,num_starts = 20)
# For theta = as.list(seq(.01,1.01,by = .05)):
# ipopt took about 24 seconds per random start (in parallel with 2 cores)
# trust.optim took about 10 seconds per additional theta
# Need to get the hessian computation reduced by an order of magnitude
# Then with additional cores, this should run quite fast.

# After speeding up hessian computation:
# 9.7 seconds per random start (ipopt; 2 cores in parallel)
# 7.75 seconds per additional theta.

# Now for a bigger theta grid (because log post was strictly increasing over range used...)
# all_theta_opt_bigger <- optimize_all_thetas(theta = as.list(seq(0.01,5.01,by = .5)),model_data = model_data,num_starts = 20)
# Initial:  259.382  seconds.
# Inner:  1504.836 seconds.
# save(all_theta_opt_bigger,file = "/home/alex1/phd/projects/case-crossover/biometrics-paper/all_theta_opt_bigger.Rdata")
# load(file = "/home/alex1/phd/projects/case-crossover/biometrics-paper/all_theta_opt_bigger.Rdata")

# Try it with ipopt as the inner optimizer. Also use a more appropriate grid
# thetagrid <- as.list((1/sqrt(seq(0.01,32.01,by = .5))))
# thetagrid <- as.list(seq(-10,10,by = .5))
# 
# 
# all_theta_opt_bigger5 <- optimize_all_thetas_parallel(theta = thetagrid,
#                                             model_data = model_data,
#                                             num_starts = 10,
#                                             random_start_sd = .1,
#                                             debug = TRUE)


# 65 values of theta in parallel with random starts took 185 seconds (3 minutes)...
# 30 values of theta in parallel with 10 random starts EACH took 32 minutes on 16 cores.
# Trying many random starts resulted in MUCH better results- modes were stable
# It does appear as though there is a problem in the computation of logpi(theta|y) now.


### COMPUTATION OF FINAL RESULT AND SUMMARIES ###

# load(file = "/home/alex1/phd/projects/case-crossover/biometrics-paper/all_theta_opt_2.Rdata")


# Take in output of optimize_all_thetas(). Add the value of the log_posterior_theta().
# This will be used for posterior summary and density computation; but we also need to plot it directly.
add_log_posterior_values <- function(optresults,model_data) {
  optresults <- ungroup(optresults)
  hessian_structure <- hessian_log_likelihood_structure(rnorm(model_data$Wd),model_data)
  optresults$sigma <- exp(-.5 * optresults$theta)
  # Log posterior for theta
  logposttheta <- optresults %>%
    purrr::pmap(~log_posterior_theta(..1,unlist(..4),model_data,hessian_structure)) %>%
    purrr::reduce(c)
  optresults$theta_logposterior <- logposttheta
  # Log posterior for sigma
  # logpostsigma <- optresults %>%
  #   purrr::pmap(~log_posterior_sigma(..7,unlist(..4),model_data,hessian_structure)) %>%
  #   purrr::reduce(c)
  optresults <- optresults %>%
    mutate(sigma = exp(-.5 * theta),
           sigma_logposterior = log(2/sigma) + logposttheta)
  optresults
}

# Function to crudely integrate posterior given log posterior values and a grid of thetas
normalize_log_posterior <- function(pp,tt) {
  # pp: log posterior
  # tt: theta values at which pp is evaluated
  # function returns the LOG of the normalizing constant
  
  # Make sure tt is sorted
  df <- tibble(pp = pp,tt = tt) %>% arrange(tt)
  tt <- df$tt
  pp <- df$pp

  lp <- length(pp)
  matrixStats::logSumExp(c(
    matrixStats::logSumExp(pp[-1] + log(diff(tt)) + log(1/2)),
    matrixStats::logSumExp(pp[-lp] + log(diff(tt)) + log(1/2))
  ))
}



# all_theta_opt_withlogpost <- add_log_posterior_values(all_theta_opt_bigger5,model_data)
# # # Normalize and plot them
# normconst <- matrixStats::logSumExp(all_theta_opt_withlogpost$theta_logposterior)
# all_theta_opt_withlogpost %>%
#   mutate(logpostnormalized = theta_logposterior - normconst,
#          post_for_plot = exp(logpostnormalized)) %>%
#   filter(theta >= -5) %>%
#   ggplot(aes(x = theta,y = post_for_plot)) +
#   theme_light() +
#   geom_line(colour = "black",size = .2) +
#   geom_point(pch = 21,colour = "black",size = 2,fill = "red",alpha = .5) +
#   labs(title = "Posterior distribution for log(theta) (log-precision of beta)",
#        x = "log(theta)",
#        y = "Posterior density")

# Function to get the diagonal of the inverse of the precision matrix

# Function to compute value of the marginal density for W_i at supplied point
compute_log_density <- function(W,i,model_results,model_data,constrA = NULL,correction = -Inf) {
  # W: vector of values of W to compute log pi(W_i|y) at
  # i: index. which element of the latent field do you want the marginal for?
  # model_results: tibble containing output of optimization
  
  # Add log posterior values for theta if not present
  if (!("theta_logposterior" %in% names(model_results))) {
    model_results <- add_log_posterior_values(model_results,model_data)
  }
  # Normalize
  thetanormconst <- matrixStats::logSumExp(model_results$theta_logposterior)
  model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
  
  # model_results$theta <- round(model_results$theta,3)
  
  # Compute values of the Gaussian approximation for each theta, for each value
  # of W supplied. 
  
  # Marginal variances
  precision_matrices <- model_results %>% 
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  hessian_structure <- hessian_log_likelihood_structure(W = model_results$solution[[1]],model_data = model_data)
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data,structure = hessian_structure),
      theta = ..1)
    )
  QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]] + exp(correction) * Diagonal(dim(.x[["Q"]])[1]),theta = .x[["theta"]]))
  margvar <- QpC %>%
    purrr::map(~diagOfInv(x = .x[["QpC"]],constrA = constrA)) %>%
    purrr::map(~.x[i])
  
  # Marginal means:
  if (is.null(constrA)) {
    # No linear constraints- just grab the right element of the latent field.
    margmeans <- model_results %>%
      purrr::pmap(~list(mean = ..4[i],theta = ..1))
  } else {
    # Correct for linear constraints
    # List of uncorrected means, one for each theta
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))
    # Compute the Q^-1A^T term. Q here is really Q + C of course...
    
    WW <- purrr::map(QpC,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
    YY <- purrr::map2(WW,uncorrectedmean,~list(
      YY = solve(t(constrA) %*% .x[["WW"]],t(constrA) %*% .y[["mode"]]),
      theta = .y[["theta"]]
    )
    )
    # Now correct
    margmeans <- pmap(list(uncorrectedmean,WW,YY),
                      ~list(mean = ..1[["mode"]] - ..2[["WW"]] %*% ..3[["YY"]],theta = ..1[["theta"]])) %>%
      map(~list(mean = .x[["mean"]][i, ],theta = .x[["theta"]]))
  }
  
  # Now compute the log densities and return
  out <- purrr::map2(margmeans,margvar,~list(dens = dnorm(W,mean = .x[["mean"]],sd = sqrt(.y),log = TRUE),theta = .x[["theta"]])) %>%
    purrr::map(~.x[["dens"]] + filter(model_results,theta == as.numeric(.x[["theta"]])) %>% pull(theta_logposterior)) %>%
    purrr::reduce(rbind)
  
  if (nrow(model_results) == 1) return(out)
  
  apply(out,2,matrixStats::logSumExp)
  
}

# compute_log_density(c(-.01,-.001,0,.001,.01),i = 6836,model_results = all_theta_opt_withlogpost,model_data = model_data)
# 
# plotfunc <- function(x,i) {
#   logdens <- compute_log_density(x,i,
#                                  model_results = all_theta_opt_withlogpost,
#                                  model_data = model_data)
# 
#   exp(logdens)
# }
# 
# 
# betapm25plt <- tibble(x = c(-0.02,.1)) %>%
#   ggplot(aes(x = x)) +
#   theme_light() +
#   stat_function(fun = function(x) plotfunc(x,i = model_data$Wd-2),colour = "purple") +
#   labs(title = "pm25 beta",
#        subtitle = "Marginal posterior",
#        x = "pm25 beta",
#        y = "Density")
# 
# betano2plt <- tibble(x = c(-0.02,.1)) %>%
#   ggplot(aes(x = x)) +
#   theme_light() +
#   stat_function(fun = function(x) plotfunc(x,i = model_data$Wd-1),colour = "purple") +
#   labs(title = "no2 beta",
#        subtitle = "Marginal posterior",
#        x = "no2 beta",
#        y = "Density")
# 
# betao3plt <- tibble(x = c(-0.02,.1)) %>%
#   ggplot(aes(x = x)) +
#   theme_light() +
#   stat_function(fun = function(x) plotfunc(x,i = model_data$Wd),colour = "purple") +
#   labs(title = "o3 beta",
#        subtitle = "Marginal posterior",
#        x = "o3 beta",
#        y = "Density")
# 
# # For comparison, plot a normal with a really low variance
# reallylowvariancenormalplot <- tibble(x = c(-0.02,.1)) %>%
#   ggplot(aes(x = x)) + 
#   theme_light() + 
#   stat_function(fun = function(x) dnorm(x,mean = 0,sd = .002),colour = "blue") +
#   labs(title = "For comparsion, normal with sd = .002",
#        x = "",y = "")
# 
# cowplot::plot_grid(
#   betapm25plt,
#   betao3plt,
#   betao3plt,
#   reallylowvariancenormalplot,
#   nrow = 2
# )

# Compute marginal means and variances only

compute_marginal_means <- function(i,model_results,model_data,constrA = NULL) {
  # i: vector of indices. which elements of the latent field do you want the marginal for?
  # model_results: tibble containing output of optimization
  # Have to deal with the fact that the means can be negative. This function basically does
  # fancy averaging.
  # constrA is a matrix of linear constraints, AW = 0 where W is the latent field.
  
  # Add log posterior values for theta if not present
  if (!("theta_logposterior" %in% names(model_results))) {
    model_results <- add_log_posterior_values(model_results,model_data)
  }
  # Normalize
  thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,model_results$theta)
  model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
  
  # Marginal means:
  if (is.null(constrA)) {
    # No linear constraints- just grab the right element of the latent field.
    margmeans <- model_results %>%
      purrr::pmap(~..4[i]) %>%
      purrr::reduce(rbind)
  } else {
    # Correct for linear constraints
    # List of uncorrected means, one for each theta
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))
    # Compute the Q^-1A^T term. Q here is really Q + C of course...
    precision_matrices <- model_results %>% 
      purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                        theta = ..1)
      )
    hessian_structure <- hessian_log_likelihood_structure(W = model_results$solution[[1]],model_data = model_data)
    hessians <- model_results %>% 
      purrr::pmap(~list(
        C = hessian_log_likelihood(W = ..4,model_data = model_data,structure = hessian_structure),
        theta = ..1)
      )
    
    QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))
    
    WW <- purrr::map(QpC,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
    YY <- purrr::map2(WW,uncorrectedmean,~list(
      YY = solve(t(constrA) %*% .x[["WW"]],t(constrA) %*% .y[["mode"]]),
      theta = .y[["theta"]]
    )
    )
    # Now correct
    margmeans <- pmap(list(uncorrectedmean,WW,YY),
                      ~..1[["mode"]] - ..2[["WW"]] %*% ..3[["YY"]]) %>%
      purrr::map(~.x[i, ]) %>%
      purrr::reduce(rbind)
  }
  
  if (nrow(model_results) == 1) return(margmeans)
  
  postvals <- exp(model_results$theta_logposterior)
  sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
}

compute_marginal_variances <- function(i,model_results,model_data,constrA = NULL,correction = -Inf) {
  # W: vector of values of W to compute log pi(W_i|y) at
  # i: index. which element of the latent field do you want the marginal for?
  # model_results: tibble containing output of optimization
  
  # Add log posterior values for theta if not present
  if (!("theta_logposterior" %in% names(model_results))) {
    model_results <- add_log_posterior_values(model_results,model_data)
  }
  # Normalize
  thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,model_results$theta)
  model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
  
  
  precision_matrices <- model_results %>% 
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  hessian_structure <- hessian_log_likelihood_structure(W = model_results$solution[[1]],model_data = model_data)
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data,structure = hessian_structure),
      theta = ..1)
    )
  
  # Marginal variances: add the precision and the hessian and get diagOfInv
  margvar <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
    purrr::map(~diagOfInv(x = .x,constrA = constrA,i = i)) %>%
    purrr::reduce(rbind)
  
  if (nrow(model_results) == 1) return(margvar)
  
  postvals <- exp(model_results$theta_logposterior)
  sweep(margvar,1,postvals,"*") %>% apply(2,sum)
}


# beta_margmeans <- rginal_means(6836:6838,model_results,model_data)
# beta_margvars <- rginal_variances(6836:6838,model_results,model_data)
# 
# cowplot::plot_grid(
#   betapm25plt + stat_function(fun = function(x) dnorm(x,mean = beta_margmeans[1],sd = sqrt(beta_margvars[1])),colour = "orange"),
#   betano2plt + stat_function(fun = function(x) dnorm(x,mean = beta_margmeans[2],sd = sqrt(beta_margvars[2])),colour = "orange"),
#   betao3plt + stat_function(fun = function(x) dnorm(x,mean = beta_margmeans[3],sd = sqrt(beta_margvars[3])),colour = "orange"),
#   nrow = 2
# )


# The computation of marginal means and variances use a bunch of the same quantities, so do them together

compute_marginal_means_and_variances <- function(i,model_results,model_data,constrA = NULL,lincomb = NULL) {
  if (nrow(model_results) > 1) {
    # Add log posterior values for theta if not present
    if (!("theta_logposterior" %in% names(model_results))) {
      model_results <- add_log_posterior_values(model_results,model_data)
    }
    # Normalize
    thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,model_results$theta)
    model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
    # Get the integration weights
    dx1 <- diff(model_results$theta) # dx1[i] = x[i+1] - x[i]
    ld <- length(dx1)
    intweights <- c(
      dx1[1]/2,
      (dx1[1:(ld-1)] + dx1[2:ld])/2,
      dx1[ld]/2
    )
  }
  
  # Compute the precision matrices for each theta
  precision_matrices <- model_results %>% 
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  
  # Compute the hessians for each theta
  hessian_structure <- hessian_log_likelihood_structure(W = model_results$solution[[1]],model_data = model_data)
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data,structure = hessian_structure),
      theta = ..1)
    )
  # If linear combinations required, set up the relevant functions
  if (!is.null(lincomb)) {
    compute_var_one_lincomb <- function(a,Q) {
      # a <- cbind(a)
      ZZ <- solve(Q,a)
      as.numeric(crossprod(a,ZZ))
    }
    compute_var_all_lincombs <- function(A,Q) {
      # Coerce to list of sparse vectors
      AA <- list()
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[ ,j],"sparseVector")
      AA %>% purrr::map(~compute_var_one_lincomb(.x,Q)) %>% purrr::reduce(c)
    }
    compute_one_lincomb_correction <- function(a,WW,VV) {
      as.numeric(crossprod(crossprod(VV,a),crossprod(WW,a)))
    }
    compute_all_lincomb_correction <- function(A,WW,VV) {
      AA <- list()
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[ ,j],"sparseVector")
      AA %>% purrr::map(~compute_one_lincomb_correction(.x,WW,VV)) %>% purrr::reduce(c)
    }
  }
  
  # If no linear constraints, compute the marginal means and variances as normal
  if (is.null(constrA)) {
    margmeans <- model_results %>%
      purrr::pmap(~..4) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)
    
    # Marginal variances: add the precision and the hessian and get diagOfInv
    margvars <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
      purrr::map(~diagOfInv(x = .x,constrA = NULL,i = i)) %>%
      purrr::reduce(rbind)
    
    # If there are linear combinations, compute their variances separately from diagOfInv
    if (!is.null(lincomb)) {
      # lincomb is a column matrix. Change to list and map over the columns
      lincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
  } else {
    # If there are linear constraints, compute the corrected mean and variance
    # First compute the uncorrected mean
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))
    
    # Get the precision matrix of the GMRF- Q + C
    QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))
    
    # Compute the correction term
    WW <- purrr::map(QpC,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
    
    YY <- purrr::map2(WW,uncorrectedmean,~list(
      YY = solve(t(constrA) %*% .x[["WW"]],t(constrA) %*% .y[["mode"]]),
      theta = .y[["theta"]]
    ))
    correction_mean <- purrr::map2(WW,YY,~list(correction = .x[["WW"]] %*% .y[["YY"]],theta = .x[["theta"]]))
    
    # Now correct
    margmeans <- purrr::map2(uncorrectedmean,correction_mean,
                             ~.x[["mode"]] - .y[["correction"]]) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)
    
    # Now compute the variances
    # Add the corrected mean to the model_results
    model_results$corrected_mean <- vector(mode = "list",length = nrow(model_results))
    for (k in 1:nrow(model_results)) model_results$corrected_mean[[k]] <- margmeans[k, ]
    
    # Re-compute the hessians
    corrected_hessians <- model_results %>%
      purrr::pmap(~list(
        C = hessian_log_likelihood(W = ..12,model_data = model_data,structure = hessian_structure),
        theta = ..1)
      )
    # Get the corrected precision matrix of the GMRF- Q + C_correct
    QpC_corrected <- purrr::map2(precision_matrices,corrected_hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))
    
    # uncorrectedvariances <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = NULL,i = i))
    margvars <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = constrA,i = i)) %>%
      purrr::reduce(rbind)
    if (!is.matrix(margvars)) margvars <- matrix(margvars,nrow = 1)
    
    # If we require marginal variances for linear combinations, compute them separately from diagOfInv
    if (!is.null(lincomb)) {
      uncorrectedlincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]]))
      
      # Compute the corrections
      WW <- purrr::map(QpC_corrected,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
      VV <- purrr::map(WW,~list(VV = solve(t(.x[["WW"]]) %*% constrA,t(.x[["WW"]])),theta = .x[["theta"]])) %>% 
        purrr::map(~list(VV = t(.x[["VV"]]),theta = .x[["theta"]]))
      
      lincombvarcorrections <- purrr::map2(WW,VV,~list(lincombvarcorrection = compute_all_lincomb_correction(lincomb,.x[["WW"]],.y[["VV"]]),
                                                       theta = .x[["theta"]]))
      
      lincombvars <- purrr::map2(uncorrectedlincombvars,lincombvarcorrections,
                                 ~list(lincombvar = .x[["lincombvar"]] - .y[["lincombvarcorrection"]],
                                       theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
    
    
  }
  if (nrow(margmeans) == 1) {
    finalmeans <- as.numeric(margmeans)[i]
    finalvars <- as.numeric(margvars)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- as.numeric(lincombvars)
    
  } else {
    postvals <- exp(model_results$theta_logposterior + log(intweights))
    finalmeans <- sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
    finalvars <- sweep(margvars,1,postvals,"*") %>% apply(2,sum)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- sweep(lincombvars,1,postvals,"*") %>% apply(2,sum)
    finalmeans <- finalmeans[i]
  }
  
  list(mean = finalmeans,
       variance = finalvars,
       lincombvars = finallincombvars)
  
}

### Misc Functions ###
# From Pmisc:
precToSd <- function (densmat) 
{
  densmat[, "y"] = densmat[, "y"] * 2 * densmat[, "x"]^(3/2)
  densmat[, "x"] = 1/sqrt(densmat[, "x"])
  densmat
}



