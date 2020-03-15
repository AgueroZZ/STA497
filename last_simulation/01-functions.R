# Get the diagonal of Q^-1 via inverting sparse cholesky decomp
diagOfInv_original = function(x, verbose=FALSE,constrA = NULL,i = NULL) {
  # x: matrix to compute the diagonal of the inverse of
  # verbose: print progress?
  # constrA: optional n x k matrix of linear constraints to correct for. n = dim(x), k = number of constraints
  # i: optional vector of indices of x for which to return the diagonal. only used if correcting
  # for linear constraints. The whole matrix x is needed even if you only want a small subset i
  # of marginal variances. However, correcting for linear constraints needs only the subset i,
  # which if small can drastically decrease run time.
  
  if (is.null(i)) i <- 1:dim(x)[1]
  
  if(verbose) {
    cat("cholesky\n")
  }
  cholHere = Matrix::expand(Matrix::Cholesky(x, 
                                             LDL=FALSE, 
                                             super=FALSE))
  if(verbose) {
    cat("solve\n")
  }
  cholHere$Linv = Matrix::solve(cholHere$L)
  
  if(verbose) {
    cat("multiply\n")
  }
  # sum the columns of Linv^2
  cholHere$LinvDf = data.table::data.table(
    col = rep(1:nrow(cholHere$Linv), diff(cholHere$Linv@p)),
    x = cholHere$Linv@x^2
  )
  
  varDiag = cholHere$LinvDf[, .(sum = sum(x)), by = col]
  
  if(verbose) {
    cat("permute\n")
  }
  
  # do the permutation transform
  varDiagMat = Diagonal(nrow(varDiag), varDiag$sum)
  varDiagMatP = crossprod(cholHere$P, varDiagMat) %*% cholHere$P
  
  if (is.null(constrA)) {
    vars <- varDiagMatP@x
  } else {
    # Correct for linear constraints
    WW <- Matrix::solve(x,constrA)
    # VV <- Matrix::solve(t(constrA) %*% WW,WW)
    VV <- t(Matrix::solve(t(WW) %*% constrA,t(WW)))
    # Correct for the ones that are actually being returned
    correction <- rep(0,length(varDiagMatP@x))
    for (j in i) correction[j] <- VV[j, ] %*% WW[j, ]
    
    vars <- varDiagMatP@x - correction
  }
  vars[i]
}


diagOfInv <- function(x, verbose=FALSE,constrA = NULL,i = NULL) {
  if (is.null(i)){
    i <- 1:dim(x)[1]
  } 
  result <- diag(solve(x))
  if (is.null(constrA)) {
    vars <- result
  } else {
    # Correct for linear constraints
    WW <- Matrix::solve(x,constrA)
    # VV <- Matrix::solve(t(constrA) %*% WW,WW)
    VV <- t(Matrix::solve(t(WW) %*% constrA,t(WW)))
    # Correct for the ones that are actually being returned
    correction <- rep(0,length(result))
    for (j in i) correction[j] <- VV[j, ] %*% WW[j, ]
    vars <- result - correction
  }
  vars[i]
}




### MODEL SETUP ###
# Function to create the difference matrix D
create_diff_matrix <- function(n) {
  # n is the total # sample size.
  cbind(Matrix(1,n-1,1),Diagonal(n-1,-1))
}
# Function to create the crossproduct's inverse,(DD^T)^(-1)
create_full_dtcp_matrix <- function(n) {
  m <- Diagonal(n-1,1) - Matrix(1/n,n-1,n-1)
  m
}

# Use to order your data first!!! It orders your dataset based on observed times:
arrange_data <- function(data){
  newdata <- arrange(data,times)
  newdata
}

#Function that creates the adjust rank that will be used later for Breslow's adjustment:
Get_Adj <- function(model_data){
  new_data <- model_data
  new_data$adjRank <- rank(model_data$times,ties.method = "min")
  new_data
}

# Function to take a covariate and return the appropriate element of "Alist"(rw2)
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

# Function to take a covariate and return the appropriate element of "Blist"(iid)
create_blist_element <- function(u,constraint = NULL) {
  # u: covariate. NOT sorted and MAY contain ties/repeated values, in general.
  # constraint: vector containing values of u for which random effect U should be
  # constrained to be zero. 
  lu <- length(u)
  B <- Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "iid"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }
  list(u = u,B = B,model = model,constrzero = constrzero)
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
#Function that creates the k observed time
# Argument: Data that contains all times and all censoring indicators (0 means right censored)
# Return: A vector that contains all the observed times in the dataset
Get_Observed_time <- function(model_data){
  observed <- model_data$times[model_data$censoring == 1]
  observed
}

#Function that determines the k observed time's indexes 
# Argument: Model data that contains all times and all censoring indicators
# Return: A vector that contains all the indexes of observed times in the dataset
Get_Observed_index <- function(model_data){
  index <- which(model_data$ID %in% model_data$ID[model_data$censoring == 1])
  index
}

#Function that creates the i-th Risk Set.
# Argument: the i-th observed time, model data that contains all entry times and observed times
# Return: the indexes of observations that are in the i-th Risk Set.
Get_Risk_Set_i <- function(timei,model_data){
  set <- timei<=model_data$times & timei>=model_data$entry
  set <- which(set)
  set
}

#Function that creates all the Risk sets using Get_Risk_Set_i function.
Get_all_Risk_Set <- function(model_data){
  observed <- Get_Observed_time(model_data)
  lapply(observed,Get_Risk_Set_i,model_data)
}

#Function used to retrieve deltas which determine the log-likelihood.
#No Delta_1,1 in neither output nor input, as it is always 0, so # of delta is n-1
prep_data_for_log_lik <- function(W,model_data) {
  # Take the parameters and model data; return a list of parameters grouped by subject.
  # Get delta
  delta <- W[1:(model_data$n-1)]
  delta
}

# New method: That also allows the Breslow's method. Still fast but not require big storage for risk sets.Please use this one.
log_likelihood <- function(W,model_data) {
  # Arguments:
  #   W: the latent field, a vector of dimension Wd. The first Nd elements are the
  #      "deltas", and are what's actually used to compute the likelihood.
  #   model_data: a list containing the output of model_setup; A and X don't need to be
  #               included (for efficiency). The dimensions of everything are what's needed.
  delta <- prep_data_for_log_lik(W,model_data)
  delta <- c(0,delta)
  ob <- Get_Observed_index(model_data)
  adj <- Get_Adj(model_data)$adjRank
  compute_loglik_for_ith_obs <- function(i) {
    # deltavec is a vector of deltas for each Risk Set.
    deltavec <- delta[adj[i]:length(delta)]
    # Note we are adding both minus signs in here.
    -log(exp(matrixStats::logSumExp(delta[i]-deltavec)))
  }
  # Now apply that function to the whole list and sum the result to get the answer
  lapply(ob,compute_loglik_for_ith_obs) %>% reduce(sum)
}

# A faster function that compute the value of the gradient of log-likelihood at a specific value of W, but does not need risk sets.
# This function corrects for ties.
grad_log_likelihood <- function(W,model_data) {
  # Arguments: see log_likelihood
  # Returns: numeric vector of the same length as W representing the gradient.
  #          It will have zeroes on the end of it (see appendix of paper).
  # Prepare the deltas that we will be using for computation of gradient.
  delta <- prep_data_for_log_lik(W,model_data)
  delta <- c(0,delta)
  # Get the index of observed times and corresponding Risk sets.
  ob <- Get_Observed_index(model_data)
  adj <- Get_Adj(model_data)$adjRank
  # Helper to compute the ith gradient term.
  compute_i_gradient_term <- function(i) {
    grad <- rep(0,model_data$n)
    deltavec <- delta[adj[ob[i]]:length(delta)]
    deno <- exp(matrixStats::logSumExp(delta[ob[i]]-deltavec))
    grad[adj[ob[i]]:length(delta)] <- exp(delta[ob[i]]-deltavec)/(deno)
    grad[ob[i]] <- -(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec))-1)/(deno)
    grad[1] <- 0
    grad
  }
  # Go ahead to use compute_i_gradient_term to compute all the gradients, and sum to a vector
  final_grad <- lapply(1:length(ob),compute_i_gradient_term)
  final_grad <- Reduce('+',final_grad)
  final_grad <- final_grad[-1]
  # The gradient is the concatenation of all these (vector) terms,
  # plus zeroes on the end to make it as long as W
  gradient_back <- rep(0,length(W) - length(final_grad))
  c(final_grad,gradient_back)
}



#A Faster version of correct way to compute C matrix: allow Breslow Correction

hessian_log_likelihood <- function(W,model_data) {
  ob <- Get_Observed_index(model_data)
  delta <- prep_data_for_log_lik(W,model_data)
  adj <- Get_Adj(model_data)$adjRank
  #Temporally fill in Delta_11 into delta, to make the future computation easier.
  delta <- c(0,delta)
  #Define a function to compute one C_i matrix:
  GetCi <- function(i,model_data=model_data,delta=delta){
    if(i == adj[i]){
      deltavec <- delta[i:model_data$n]
      denom <- exp(matrixStats::logSumExp(deltavec[1]-deltavec))
      M <- matrix(0,nrow = length(deltavec), ncol = length(deltavec))
      for (j in 1:length(deltavec)) {
        for (k in j:length(deltavec)) {
          if(j==k & k==1){
            M[j,k] <- -(exp(matrixStats::logSumExp(deltavec[1]-deltavec))-1)/(denom^2)
          }
          else if(j==k & k!=1){
            deltavec1 <- deltavec[c(-1,-j)]
            if(length(deltavec1)==0){
              M[j,k] <- -(exp(deltavec[1]-deltavec[j]))/(denom^2)
            }
            if(length(deltavec1)!=0){
              M[j,k] <- -((exp(matrixStats::logSumExp(2*deltavec[1]-deltavec[j]-deltavec1)))+exp(deltavec[1]-deltavec[j]))/(denom^2)
            }
          }
          else if(j==1 & k!=1){
            M[j,k] <- (exp(deltavec[1]-deltavec[k]))/(denom^2)
          }
          else if(k==1 & j!=1){
            M[j,k] <- (exp(deltavec[1]-deltavec[j]))/(denom^2)
          }
          else {
            M[j,k]<- (exp(2*deltavec[1]- deltavec[j]-deltavec[k]))/(denom^2)
          }
        }
      }
      M <- bdiag(matrix(0,ncol = (model_data$n-ncol(M)),nrow = (model_data$n-nrow(M))),M)
      M
    }
    else{
      deltavec <- delta[adj[i]:model_data$n]
      dif <- i - adj[i]
      c_i <- 1+dif
      denom <- exp(matrixStats::logSumExp(deltavec[c_i]-deltavec))
      M <- matrix(0,nrow = length(deltavec), ncol = length(deltavec))
      for (j in 1:length(deltavec)) {
        for (k in j:length(deltavec)) {
          if(j==k & k== c_i){
            M[j,k] <- -(exp(matrixStats::logSumExp(deltavec[c_i]-deltavec))-1)/(denom^2)
          }
          else if(j==k & k!= c_i){
            deltavec1 <- deltavec[c(-c_i,-j)]
            if(length(deltavec1)==0){
              M[j,k] <- -(exp(deltavec[c_i]-deltavec[j]))/(denom^2)
            }
            if(length(deltavec1)!=0){
              M[j,k] <- -((exp(matrixStats::logSumExp(2*deltavec[c_i]-deltavec[j]-deltavec1)))+exp(deltavec[c_i]-deltavec[j]))/(denom^2)
            }
          }
          else if(j==1 & k!=c_i){
            M[j,k] <- (exp(deltavec[c_i]-deltavec[k]))/(denom^2)
          }
          else if(k==1 & j!=c_i){
            M[j,k] <- (exp(deltavec[c_i]-deltavec[j]))/(denom^2)
          }
          else {
            M[j,k]<- (exp(2*deltavec[c_i]- deltavec[j]-deltavec[k]))/(denom^2)
          }
        }
      }
      M <- bdiag(matrix(0,ncol = (model_data$n-ncol(M)),nrow = (model_data$n-nrow(M))),M)
      M
    }
  }
  R <- lapply(ob, GetCi, model_data=model_data,delta=delta)
  C_final <- Reduce('+',R)
  C_final <- bdiag(C_final[2:nrow(C_final),2:ncol(C_final)],diag(rep(0,model_data$Wd-model_data$Nd),nrow =model_data$Wd-model_data$Nd))
  C_final <- as(C_final,"dgCMatrix")
  return(-forceSymmetric(C_final))
}






### LATENT FIELD ###

# Functions to implement the precision matrix for
# 1) Linear terms only
# 2) IID Random effect only
# 3) RW2 Random effect terms only
# 4) Two of the above
# 5) All of the above


# Linear terms only
Q_matrix_linear <- function(theta,model_data,tau = exp(7),debug = FALSE) {
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

# iid random effect term
Q_matrix_iid_one_component <- function(theta,model_data,covariate) {
  # theta: log(precision); log of 1/(iid random effect variance)
  # covariate: character, name of covariate as appears in model_data$B
  # This function creates the component of the Q matrix corresponding to a single covariate.
  # model_data contains an element "B" which is a list; names(B) is the vector of names of
  # covariates to be modelled as iid random effect. B$covariate is itself a list containing Bmat, the 
  # random effects design matrix, and u, the vector containing the UNIQUE values of the covariate. The columns 
  # of B$covariate$Bmat should already be ordered according to order(unique(u)); u can thus be either
  # sorted or not. It will be sorted inside this function
  # u should ALREADY be unique and sorted. This just makes this explicit in the code.
  theta2 <- theta[length(theta)]
  u <- sort(unique(model_data$B[[covariate]]$u))
  ul <- length(u)
  AA <- Matrix::Diagonal(n = ul,1)
  exp(theta2) * AA
}


Q_matrix_iid <- function(theta,model_data,tau = exp(7)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$B)) stop("no iid components in model")
  
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)
  
  howmanyiid <- length(whichiid)
  
  
  Suinv <- purrr::map2(whichiid,1:howmanyiid,
                       ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)
  
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Bd),
    cbind(-tau*t(Bd)%*% model_data$lambdainv,Suinv + tau * crossprod(Bd,crossprod(model_data$lambdainv,Bd)))
  )
}

# RW2 Models: random effects design matrix and precision matrix
Q_matrix_rw2_one_component <- function(theta,model_data,covariate) {
  # theta1: log(precision); log of 1/(rw2 smoothing variance)
  # covariate: character, name of covariate as appears in model_data$A
  # This function creates the component of the Q matrix corresponding to a single covariate.
  # model_data contains an element "A" which is a list; names(A) is the vector of names of
  # covariates to be modelled smoothly. A$covariate is itself a list containing Amat, the 
  # random effects design matrix, and u, the vector containing the UNIQUE values of the covariate. The columns 
  # of A$covariate$Amat should already be ordered according to order(unique(u)); u can thus be either
  # sorted or not. It will be sorted inside this function
  # u should ALREADY be unique and sorted. This just makes this explicit in the code.
  theta1 <- theta[1]
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
  
  exp(theta1) * forceSymmetric(crossprod(H,crossprod(AA,H)))
}


Q_matrix_rw2 <- function(theta,model_data,tau = exp(7), buffer = 1/exp(7)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
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
  
  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(10)))
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }
  
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Ad),
    cbind(-tau*t(Ad)%*% model_data$lambdainv,Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad)))
  )
}

# Q matrix for model containing both linear and rw2 terms
Q_matrix_both_rw2_linear <- function(theta,model_data,tau = exp(7) , buffer = 1/exp(7)) {

  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
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
  
  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(12)))
  
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
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

Q_matrix_both_rw2_iid <- function(theta,model_data,tau = exp(7),buffer = 1/exp(7)) {


  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
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
  
  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(12)))
  
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }
  
  # IID
  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)
  
  howmanyiid <- length(whichiid)
  # The thetas are in order.
  
  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)
  
  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% Bd
  Q22 <- Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * crossprod(Ad,crossprod(model_data$lambdainv,Bd))
  Q33 <- Suinv2 + tau*crossprod(crossprod(model_data$lambdainv,Bd),Bd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}


Q_matrix_both_iid_linear <- function(theta,model_data,tau = exp(7)) {

  
  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)
  
  howmanyiid <- length(whichiid)
  # The thetas are in order.
  
  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)
  
  
  
  # linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))
  
  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Bd
  Q13 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv2 + tau * crossprod(Bd,crossprod(model_data$lambdainv,Bd))
  Q23 <- tau * crossprod(Bd,crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}


# Q matrix for model containing all of linear, rw2 and iid terms
Q_matrix_all <- function(theta,model_data,tau = exp(7)) {
  
  if (is.null(model_data$A)) stop("no rw2 components in model")
  
  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)
  
  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.
  
  
  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
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
  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(12)))
  
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }
  
  # Linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))
  
  #IID random effect
  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)
  
  howmanyiid <- length(whichiid)
  # The thetas are in order.
  
  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()
  
  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)
  
  
  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% Bd
  Q14 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv + tau * crossprod(Ad,crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * crossprod(Ad,crossprod(model_data$lambdainv,Bd))
  Q24 <- tau * crossprod(Ad,crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Suinv2 + tau * crossprod(Bd,crossprod(model_data$lambdainv,Bd))
  Q34 <- tau * crossprod(Bd,crossprod(model_data$lambdainv,model_data$Xd))
  Q44 <- Sbinv + tau*crossprod(crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13,Q14),
    cbind(t(Q12),Q22,Q23,Q24),
    cbind(t(Q13),t(Q23),Q33,Q34),
    cbind(t(Q14),t(Q24),t(Q34),Q44)
  )
}



# Make a general function to compute the Q matrix for a model
# This can be written to just call one of the above... somehow. Maybe as an option inside model_data
Q_matrix <- function(theta,model_data,tau = exp(7),forcesymm = TRUE) {
  # theta: vector of hyperparameters. The structure of this will depend on the model,
  # as specified by model_data
  if (is.null(model_data$A)) {
    if (model_data$p == 0) {stop("both X and A are null...")}
    else if(is.null(model_data$B)){
      mat <- Q_matrix_linear(theta,model_data,tau)
    }
    else{
      mat <- Q_matrix_both_iid_linear(theta,model_data,tau)
    }
  }
  else {
    if (model_data$p > 0) {
      if(is.null(model_data$B)){
        mat <- Q_matrix_both_rw2_linear(theta,model_data,tau)
      }
      else{
        mat <- Q_matrix_all(theta,model_data,tau)
      }
    }
    else {
      if(is.null(model_data$B)){
        mat <- Q_matrix_rw2(theta,model_data,tau)
      }
      # warning("random walk-only models are rank deficient for case crossover. adding fudge factor. consider adding a linear term.")
      else{
        mat <- Q_matrix_both_rw2_iid(theta,model_data,tau)
      }
      # mat <- mat + (1/tau) * Diagonal(dim(mat)[1])    
    }
  }
  if (forcesymm) return(forceSymmetric(mat))
  mat
}

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

# Hessian of log posterior.
# Option to pass in Q as well
hessian_log_posterior_W <- function(W,theta = NULL,Q = NULL,model_data) {
  if (is.null(theta) & is.null(Q)) stop("One of Q or theta must be provided")
  if (is.null(Q)) Q <- Q_matrix(theta,model_data)
  A <- as.matrix(-(Q + hessian_log_likelihood(W,model_data)))
  A <- as(A,'dgCMatrix')
  A
}






### HYPERPARAMETERS ###

# Function to compute the log posterior approximation for theta
# model_data will have a function theta_logprior() which takes vector theta
# and evaluates the log of the joint prior
log_posterior_theta_withbuffer <- function(theta,W,model_data,Q = NULL, buffer = 1/exp(10)) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  Qe <- eigen(Q,symmetric = TRUE,only.values = TRUE)$values
  Qe[length(Qe)] <- Qe[length(Qe)] + buffer
  term2_det <- (1/2) * sum(log(Qe))
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  term1 <- log_likelihood(W,model_data)
  # dt <- determinant(Q,logarithm = TRUE) # For this, we DO need the determinant (original)
  # term2_det <- (1/2) * as.numeric(dt$modulus * dt$sign)
  # term2_det <- (1/2) * as.numeric(dt$modulus) (original)
  term2 <- logprior_W(W,theta,model_data) # Doesn't contain the determinant
  term3 <- model_data$theta_logprior(theta)
  qcdet <- determinant(Q_p_C,logarithm = TRUE)
  # term4 <- -(1/2)*as.numeric(qcdet$modulus * qcdet$sign) # The gaussian approx evaluated at conditional mode
  term4 <- -(1/2) * as.numeric(qcdet$modulus) # (original)
  as.numeric(term1 + term2_det + term2 + term3 + term4)
}

log_posterior_theta <- function(theta,W,model_data,Q = NULL) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  term1 <- log_likelihood(W,model_data)
  dt <- determinant(Q,logarithm = TRUE) # For this, we DO need the determinant (original)
  # term2_det <- (1/2) * as.numeric(dt$modulus * dt$sign)
  term2_det <- (1/2) * as.numeric(dt$modulus) #(original)
  term2 <- logprior_W(W,theta,model_data) # Doesn't contain the determinant
  term3 <- model_data$theta_logprior(theta)
  qcdet <- determinant(Q_p_C,logarithm = TRUE)
  # term4 <- -(1/2)*as.numeric(qcdet$modulus * qcdet$sign) # The gaussian approx evaluated at conditional mode
  term4 <- -(1/2) * as.numeric(qcdet$modulus) # (original)
  as.numeric(term1 + term2_det + term2 + term3 + term4)
}


# Function to compute the log posterior for sigma, the standard deviation
# sigma = exp(-.5 * theta)
# jacobian is 2/sigma
log_posterior_sigma <- function(sigma,W,model_data,Q = NULL) {
  length(sigma)* log(2) - sum(log(sigma)) +
    log_posterior_theta(theta = -2 * log(sigma),
                        W = W,
                        model_data = model_data,
                        Q = Q)
}













### OPTIMIZATION ###
optlist_to_tibble <- function(lst) {
  # If lst is a single element of opt, wrap it in a list so the below works
  if ("optimizer" %in% names(lst)) lst <- list(lst)
  dplyr::tibble(
    theta = purrr::map(lst,"theta"),
    optimizer = purrr::map(lst,"optimizer"),
    starting = purrr::map(lst,"starting"),
    solution = purrr::map(lst,"solution"),
    function_value = purrr::map_dbl(lst,"function_value"),
    iterations = purrr::map_dbl(lst,"iterations")
  )
}

#' Optimize the conditional posterior of the latent Gaussian variables for one hyperparameter configuration
#'
#' @description Compute the conditional mode W-hat for one value of the hyperparameters theta. Uses
#' trust region methods implemented in trustOptim::trust.optim. Makes use of the sparsity of the Hessian.
#'
#' @param theta Vector, representing one configuration of hyperparameters (optimize_latentfield_trustoptim). List of vectors, each
#' representing one configuration of hyperparameters (optimize_all_thetas_parallel).
#' @param model_data ccmodeldata object as output by model_setup().
#' @param Q Optional, provide pre-computed Q matrix for this theta. Will be computed if not provided.
#' @param optcontrol Optional. Override the casecrossover default control parameters for trust.optim(). See cc_control()$opt_control.
#'
#' @export
#'
optimize_latentfield_trustoptim <- function(theta,model_data,Q = NULL,optcontrol = NULL) {
  
  # Zero is the prior mean. It's a reasonable starting value.
  # My experience has been this particular optimization problem is pretty robust to this.
  startingvals <- rep(0,model_data$Wd)
  
  # Set up the functions. Negated, functions of W only with other arguments fixed.
  optfunc <- function(W) -1 * log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) -1 * grad_log_posterior_W(W,theta,model_data,Q)
  
  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  # Get the hessian
  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  
  
  
  # Set default control arguments if not provided.
  # Note: this is really for custom development use only. In the main casecrossover() function, these are set using
  # the cc_control() function.
  if (is.null(optcontrol)) {
    optcontrol <- cc_control()$opt_control
  }
  
  # Perform the optimization
  
  opt <- trustOptim::trust.optim(
    x = startingvals,
    fn = optfunc,
    gr = optfuncgrad,
    method = "SR1",
    control = optcontrol
  )
  
  
  # Return a custom-formatted list with optimization results.
  # This is done to enable easy stacking of one of these per theta value in a dataframe.
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

#' Find the conditional model for a user-defined list of theta values
#'
#' @description This is a wrapper for optimize_latentfield_trustoptim() that computes the
#' conditional mode in parallel over a user-specified grid of theta values.
#'
#' Parallelization only works on linux and mac; the "parallel" package needs to work. If you're on windows, you'll just
#' have to wait a little longer (long enough to go to the store and buy a serious computer maybe?).
#'
#' @rdname optimize_latentfield_trustoptim
#' @param thetagrid A NIGrid object which contains the grid of theta values to be optimized for.
#' @param doparallel Logical. Use parallel::mclapply (TRUE) or lapply (FALSE) for execution? Former executes in parallel, latter in series. Default TRUE.
#'
#' @export
#'
optimize_all_thetas_parallel <- function(thetagrid,model_data,optcontrol = NULL,doparallel = TRUE) {
  
  # Check thetagrid is formatted correctly
  if (!inherits(thetagrid,"NIGrid")) stop("theta should be a NIGrid object returned by mvQuad::createNIgrid()")
  # Create the theta list
  theta <- split(mvQuad::getNodes(thetagrid),rep(1:nrow(mvQuad::getNodes(thetagrid)),ncol(mvQuad::getNodes(thetagrid))))
  
  if (!all(purrr::map_lgl(theta,is.numeric))) stop("theta should be a list of numeric vectors")
  thetalengths <- purrr::map_dbl(theta,length)
  if (length(unique(thetalengths)) != 1) stop("Make sure all thetas are the same length.")
  
  if (is.null(optcontrol)) {
    optcontrol <- model_data$control$opt_control
  }
  
  do_opt <- function(theta) {
    cat("Optimizing with theta =",theta,"\n")
    optimize_latentfield_trustoptim(theta = theta,
                                    model_data = model_data,
                                    optcontrol = optcontrol)
  }
  
  cat("Performing optimization, start time: ", format(Sys.time(),"%H:%M:%S"),"\n")
  
  tm <- proc.time()
  if (length(model_data$modelspec$model) > 0 & doparallel) {
    opt <- parallel::mclapply(theta,do_opt)
  } else {
    opt <- lapply(theta,do_opt)
  }
  opt_time <- proc.time() - tm
  opt_time <- unname(opt_time["elapsed"])
  cat("Time taken for optimization with ",length(theta)," values of theta:",opt_time,"seconds.\n")
  
  out <- optlist_to_tibble(opt)
  attr(out,"thetagrid") <- thetagrid
  out
}








### COMPUTATION OF FINAL RESULT AND SUMMARIES ###
# Take in output of optimize_all_thetas(). Add the value of the log_posterior_theta().
# This will be used for posterior summary and density computation; but we also need to plot it directly.
add_log_posterior_values <- function(optresults,model_data) {
  optresults <- dplyr::ungroup(optresults)
  # Log posterior for theta
  logposttheta <- optresults %>%
    purrr::pmap(~log_posterior_theta(unlist(..1),unlist(..4),model_data)) %>%
    as.numeric()
  optresults$theta_logposterior <- logposttheta
  
  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
                  sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])
  
  attr(out,"thetagrid") <- attributes(optresults)$thetagrid
  out
}

normalize_optresults_logpost <- function(optresults) {
  # Get the nodes and weights corresponding to thetas in optresults
  nodesandweights <- cbind(mvQuad::getNodes(model_data$thetagrid),mvQuad::getWeights(model_data$thetagrid))
  K <- ncol(nodesandweights) - 1
  colnames(nodesandweights) <- c(stringr::str_c("theta",1:K),"weights")
  nodesandweights <- as.data.frame(nodesandweights)
  
  thetaopt <- cbind(purrr::reduce(optresults$theta,rbind),optresults$theta_logposterior)
  colnames(thetaopt) <- c(stringr::str_c("theta",1:K),"theta_logposterior")
  thetaopt <- as.data.frame(thetaopt)
  
  suppressMessages(thetaoptmerged <- dplyr::left_join(thetaopt,nodesandweights))
  
  ww <- thetaoptmerged$weights
  pp <- thetaoptmerged$theta_logposterior
  
  # thetanormconst <- matrixStats::logSumExp(pp + log(ww))
  sp <- matrixStats::logSumExp(pp[ww > 0] + log(ww[ww > 0]))
  sn <- matrixStats::logSumExp(pp[ww < 0] + log(abs(ww)[ww < 0]))
  thetanormconst <- sp - log(1 + exp(sn-sp))
  optresults$theta_logposterior <- optresults$theta_logposterior - thetanormconst
  
  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
                  sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])
  
  out$weights <- ww
  
  out
}

normalize_log_posterior <- function(pp,tt) {
  # tt: grid returned by mvQuad::createNIGrid
  # pp: log posterior evaluated at these points
  ww <- mvQuad::getWeights(tt)
  matrixStats::logSumExp(log(ww) + pp)
}


# Get the marginal distributions of theta and sigma
marginal_hyperparameter_posterior <- function(j,optresults,quantiles = c(2.5,97.5)/100) {
  thetagridfull <- model_data$thetagrid
  S <- thetagridfull$dim
  # If it's already one-dimensional, don't need to do anything new, but do compute quantiles
  if (S == 1) {
    outmat <- dplyr::tibble(
      theta = purrr::reduce(optresults$theta,rbind),
      thetalogmargpost = optresults$theta_logposterior,
      sigma = purrr::reduce(optresults$sigma,rbind),
      sigmalogmargpost = optresults$sigma_logposterior
    )
    
    thetacumsum <- cumsum(mvQuad::getWeights(thetagridfull) * exp(outmat$thetalogmargpost))
    thetaquantiles <- purrr::map(quantiles,~outmat$theta[which(thetacumsum == min(thetacumsum[thetacumsum > .x]))]) %>% as.numeric()
    sigmaquantiles <- rev(exp(-.5 * thetaquantiles))
    
    return(
      list(
        margpost = outmat,
        quantiles = dplyr::tibble(whichmarginal = rep(1,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
      )
    )
  }
  # Get the reduced grid
  thetagridreduced <- mvQuad::createNIGrid(
    dim = thetagridfull$dim - 1,
    type = thetagridfull$type[-j],
    level = as.numeric(thetagridfull$level[ ,-j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagridreduced,domain = thetagridfull$features$domain[-j, ])
  
  # Get a 1-d grid, for computing quantiles at the end.
  thetagrid1d <- mvQuad::createNIGrid(
    dim = 1,
    type = thetagridfull$type[j],
    level = as.numeric(thetagridfull$level[ ,j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagrid1d,domain = thetagridfull$features$domain[j, ])
  
  # Return the marginal posterior and its evaluation points
  # In the optimization results, we have a matrix of theta values which matches the full grid,
  # and the log posterior evaluated at these values.
  nodesfull <- purrr::reduce(optresults$theta,rbind)
  nodesfull <- cbind(nodesfull,optresults$theta_logposterior) # Theta logposterior is now the last column
  nodesfull <- nodesfull[order(nodesfull[ ,j]), ]
  colnames(nodesfull) <- c(paste0("theta",1:S),"thetalogpost")
  nodesfull <- as.data.frame(nodesfull)
  
  # Now we have a matrix of thetavalues, nodes, and function values-- add on the weights
  nodesmulti <- mvQuad::getNodes(thetagridreduced)
  nodesmulti <- cbind(nodesmulti,mvQuad::getWeights(thetagridreduced))
  colnames(nodesmulti) <- c(paste0("theta",(1:S)[-j]),"weights")
  nodesmulti <- as.data.frame(nodesmulti)
  
  suppressMessages({# It prints what it's joining by, which is all columns, and I don't want to see this printed
    thetamargposts <- dplyr::left_join(nodesfull,nodesmulti) %>%
      dplyr::group_by(.data[[stringr::str_c("theta",j)]]) %>%
      dplyr::summarize(thetalogmargpost = matrixStats::logSumExp(.data[["thetalogpost"]] + log(.data[["weights"]])))
  })
  
  thetamargposts$whichmarginal <- rep(j,nrow(thetamargposts))
  # Now add on the sigmas
  outmat <- thetamargposts %>%
    dplyr::mutate(sigma = exp(-.5 * .data[[paste0("theta",j)]]),
                  sigmalogmargpost = log(2/.data[["sigma"]]) + .data[["thetalogmargpost"]]
    ) %>%
    dplyr::rename(theta = .data[[paste0("theta",j)]])
  
  # Mean and sd
  ww <- mvQuad::getWeights(thetagridfull)[ ,1]
  thetamean <- apply(ww * purrr::reduce(optresults$theta,rbind) * exp(optresults$theta_logposterior),2,sum)
  sigmamean <- apply(ww * purrr::reduce(optresults$sigma,rbind) * exp(optresults$theta_logposterior),2,sum)
  thetasd <- sqrt(apply(ww * (purrr::reduce(optresults$theta,rbind) - thetamean)^2 * exp(optresults$theta_logposterior),2,sum))
  sigmasd <- sqrt(apply(ww * (purrr::reduce(optresults$sigma,rbind) - sigmamean)^2 * exp(optresults$theta_logposterior),2,sum))
  
  
  # Quantiles
  thetacumsum <- cumsum(mvQuad::getWeights(thetagrid1d) * exp(outmat$thetalogmargpost))
  thetaquantiles <- purrr::map(quantiles,~outmat$theta[min(which(thetacumsum == min(thetacumsum[thetacumsum > .x])))]) %>% purrr::reduce(c)
  sigmaquantiles <- rev(exp(-.5 * thetaquantiles))
  
  list(
    margpost = outmat,
    margmoments = dplyr::tibble(moment = c("mean","sd"),theta = c(thetamean[j],thetasd[j]),sigma = c(sigmamean[j],sigmasd[j])),
    quantiles = dplyr::tibble(whichmarginal = rep(j,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
  )
}




compute_marginal_means_and_variances <- function(i=NULL,model_results,model_data,constrA = NULL,lincomb = NULL) {
  if (nrow(model_results) > 1) {
    # Add log posterior values for theta if not present
    if (!("theta_logposterior" %in% names(model_results))) {
      model_results <- add_log_posterior_values(model_results,model_data)
    }
    # Normalize
    thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,attributes(model_results)$thetagrid)
    model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
    # Get the integration weights
    intweights <- mvQuad::getWeights(attributes(model_results)$thetagrid)[ ,1]
  }
  # Compute the precision matrices for each theta
  precision_matrices <- model_results %>% 
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  # Compute the hessians for each theta
  hessians <- list()
  myhes <- mclapply(model_results$solution, hessian_log_likelihood,model_data = model_data,mc.cores = cores)
  for (j in 1:length(model_results$theta)) {
    hessians[[j]] <- list(C=myhes[[j]],theta=model_results$theta[j])
  }
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
      result <- c()
      for (j in 1:ncol(lincomb)){
        AA[[j]] <- as(lincomb[,j],"sparseVector")
        result[j] <- compute_var_one_lincomb(AA[[j]],Q)
      }
      result
    }
    compute_one_lincomb_correction <- function(a,WW,VV) {
      as.numeric(crossprod(crossprod(VV,a),crossprod(WW,a)))
    }
    compute_all_lincomb_correction <- function(A,WW,VV) {
      AA <- list()
      result <- c()
      for (j in 1:ncol(lincomb)){
        AA[[j]] <- as(lincomb[,j],"sparseVector")
        result[j] <- compute_one_lincomb_correction(AA[[j]],WW,VV)
      }
      result
    }
  }
  # If no linear constraints, compute the marginal means and variances as normal
  if (is.null(constrA)) {
    margmeans <- model_results %>%
      purrr::pmap(~..4) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)
    # Marginal variances: add the precision and the hessian and get diagOfInv
    if (is.null(i)) i <- 1:model_data$Wd
    margvars <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
      purrr::map(~diag(solve(.x))[i]) %>%
      purrr::reduce(rbind)
    # If there are linear combinations, compute their variances separately from diagOfInv
    if (!is.null(lincomb)) {
      # lincomb is a column matrix. Change to list and map over the columns
      lincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
  } 
  else {
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
    corrected_hessians <- list()
    myhes2 <- mclapply(model_results$corrected_mean, hessian_log_likelihood,model_data = model_data,mc.cores = cores)
    for (j in 1:length(model_results$theta)) {
      corrected_hessians[[j]] <- list(C=myhes2[[j]],theta=model_results$theta[j])
    }
    # Get the corrected precision matrix of the GMRF- Q + C_correct
    QpC_corrected <- purrr::map2(precision_matrices,corrected_hessians,~list(QpC = as((as.matrix(.x[["Q"]] + .y[["C"]])),"sparseMatrix"),theta = .x[["theta"]]))
    # uncorrectedvariances <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = NULL,i = i))
    margvars <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = constrA,i = i)) %>%
      purrr::reduce(rbind)
    
    
    if (!is.matrix(margvars)) margvars <- matrix(margvars,nrow = 1)
    # If we require marginal variances for linear combinations, compute them separately
    if (!is.null(lincomb)) {
      uncorrectedlincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]]))
      
      # Compute the corrections
      WW <- purrr::map(QpC_corrected,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
      VV <- purrr::map(WW,~list(VV = solve(t(.x[["WW"]]) %*% constrA,t(.x[["WW"]])),theta = .x[["theta"]])) %>% 
        purrr::map(~list(VV = t(.x[["VV"]]),theta = .x[["theta"]]))
      
      lincombvarcorrections <- list()
      for (jj in 1:length(WW)) {
        lincombvarcorrections[[jj]] <- list(lincombvarcorrection = compute_all_lincomb_correction(lincomb,WW[[jj]]$WW,VV[[jj]]$VV),theta = WW[[jj]]$theta)
      }
      
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
    if (!is.null(lincomb)) 
      finallincombvars <- as.numeric(lincombvars)
  } 
  else {
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
