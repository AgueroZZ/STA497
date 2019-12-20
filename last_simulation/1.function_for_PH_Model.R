# Get the diagonal of Q^-1 via inverting sparse cholesky decomp
diagOfInv = function(x, verbose=FALSE,constrA = NULL,i = NULL) {
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





### MODEL SETUP ###
# Function to create the difference matrix D
create_diff_matrix <- function(n) {
  # n is the total # sample size.
  cbind(Matrix(1,n-1,1),Diagonal(n-1,-1))
}
# Function to create the crossproduct's inverse,(DD^T)^(-1)
create_full_dtcp_matrix <- function(n) {
  Diagonal(n-1,1) - Matrix(1/(n+1),n-1,n-1)
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

# Compute the log-likelihood for data y and parameters W. It is fast but needs to store a lot huge Risk sets.(It is obsolete)
old_log_likelihood <- function(W,model_data) {
  # Arguments:
  #   W: the latent field, a vector of dimension Wd. The first Nd elements are the
  #      "deltas", and are what's actually used to compute the likelihood.
  #   model_data: a list containing the output of model_setup; A and X don't need to be
  #               included (for efficiency). The dimensions of everything are what's needed.
  delta <- prep_data_for_log_lik(W,model_data)
  delta <- c(0,delta)
  ob <- Get_Observed_index(model_data)
  Risk <- Get_all_Risk_Set(model_data)
  compute_loglik_for_ith_Risk <- function(i) {
    # deltavec is a vector of deltas for each Risk Set.
    deltavec <- delta[Risk[[i]]]
    # Note we are adding both minus signs in here.
    -log(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec)))
  }
  # Now apply that function to the whole list and sum the result to get the answer
  lapply(1:length(ob),compute_loglik_for_ith_Risk) %>% reduce(sum)
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

# (obsolete) Function that compute the value of the gradient of log-likelihood at a specific value of W. Fast but requires risk sets.(This is obsolete)
old_grad_log_likelihood <- function(W,model_data) {
  # Arguments: see log_likelihood
  # Returns: numeric vector of the same length as W representing the gradient.
  #          It will have zeroes on the end of it (see appendix of paper).
  
  # Prepare the deltas that we will be using for computation of gradient.
  delta <- prep_data_for_log_lik(W,model_data)
  delta <- c(0,delta)
  # Get the index of observed times and corresponding Risk sets.
  ob <- Get_Observed_index(model_data)
  Risk <- Get_all_Risk_Set(model_data)
  # Helper to compute the ith gradient term.
  compute_i_gradient_term <- function(i) {
    grad <- rep(0,model_data$n)
    deltavec <- delta[Risk[[i]]]
    grad[Risk[[i]]] <- exp(delta[ob[i]]-deltavec)/(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec)))
    grad[ob[i]] <- -(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec))-1)/(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec)))
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

# Fastest function that computes the gradient, but does not correct for ties.
grad_log_likelihood_without_ties <- function(W,model_data) {
  # Arguments: see log_likelihood
  # Returns: numeric vector of the same length as W representing the gradient.
  #          It will have zeroes on the end of it (see appendix of paper).
  # Prepare the deltas that we will be using for computation of gradient.
  delta <- prep_data_for_log_lik(W,model_data)
  delta <- c(0,delta)
  # Get the index of observed times and corresponding Risk sets.
  ob <- Get_Observed_index(model_data)
  # Helper to compute the ith gradient term.
  compute_i_gradient_term <- function(i) {
    grad <- rep(0,model_data$n)
    deltavec <- delta[ob[i]:length(delta)]
    grad[ob[i]:length(delta)] <- exp(delta[ob[i]]-deltavec)/(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec)))
    grad[ob[i]] <- -(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec))-1)/(exp(matrixStats::logSumExp(delta[ob[i]]-deltavec)))
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

#Function that computes the hessian matrix of log-likelihood function with respect to latent field at a specific value.
#W is the values of latent field that we evaluate. (Please use this one)
hessian_log_likelihood_correct <- function(W,model_data) {
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
      Risk <- Get_all_Risk_Set(model_data)
      ob <- Get_Observed_index(model_data)
      #Function 1: to be used
      ComputeC_ij <- function(M,i,r,Risk,delta,model_data){
        deltavec <- delta[Risk[[r]]]
        vec <- M[i,]
        j <- as.numeric(vec[1])
        k <- as.numeric(vec[2])
        C <- 0
        case <- 0
        if(j==k & k==ob[r]){
          case <- 1
          C <- -(exp(matrixStats::logSumExp(delta[ob[r]]-deltavec))-1)/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
        }
        else if(j==k & k!=ob[r]){
          case <- 2
          deltavec1 <- deltavec[Risk[[r]]!=j & Risk[[r]]!=ob[r]]
          if(length(deltavec1)==0){
            C <- -(exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
          }
          if(length(deltavec1)!=0){
            C <- -((exp(matrixStats::logSumExp(2*delta[ob[r]]-delta[j]-deltavec1)))+exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
          }
        }
        else if(j==ob[r] & k!=ob[r]){
          case <- 3
          C <- (exp(delta[ob[r]]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
        }
        else if(j!=ob[r] & k==ob[r]){
          case <- 4
          C <- (exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
        }
        else{
          case <- 5
          C <- (exp(2*delta[ob[r]]-delta[j]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
        }
        C
      }
      # Function 2:to be used
      Compute_full_C_i <- function(M,r,model_data,Risk,delta){
        C_i <- matrix(data = 0,nrow = model_data$n, ncol = model_data$n)
        if(nrow(M)==0){
          C_i
        }
        else{
          values <- lapply(1:nrow(M),ComputeC_ij, M=M,r=r,Risk=Risk,model_data=model_data, delta=delta)
          for (i in 1:nrow(M)) {
            #Updating entries of C_i
            v <- M[i,]
            C_i[v[1],v[2]] <- values[[i]]
          }
          C_i
        }
      }
      #Used the two function above to get C_i
      M <- as.matrix(expand.grid(Risk[[which(ob %in% i)]],Risk[[which(ob %in% i)]]))
      M <- M[M[,1]<=M[,2],,drop = FALSE]
      C_i <- Compute_full_C_i(M,r=which(ob %in% i),model_data = model_data,delta = delta,Risk = Risk)
      C_i
    }
  }
  R <- lapply(ob, GetCi, model_data=model_data,delta=delta)
  C_final <- Reduce('+',R)
  C_final <- bdiag(C_final[2:nrow(C_final),2:ncol(C_final)],diag(rep(0,model_data$Wd-model_data$Nd),nrow =model_data$Wd-model_data$Nd))
  C_final <- as(C_final,"dgCMatrix")
  return(-forceSymmetric(C_final))
}

#Note that it uses a very big storage for the n of n*n matrices.(This is obsolete)
old_hessian_log_likelihood <- function(W,model_data) {
  Risk <- Get_all_Risk_Set(model_data)
  ob <- Get_Observed_index(model_data)
  delta <- prep_data_for_log_lik(W,model_data)
  #Temporally fill in Delta_11 into delta, to make the future computation easier.
  delta <- c(0,delta)
  #Compute the i-th observed time's Hessian matrix, C_i:
  #First Step: Define a function to compute a single entry of C_i.(r here is the index of Risk set that we are working on)
  #So we want to make sure that r-th risk set is matched with index of its corresponding observation, which is ob[r].
  ComputeC_ij <- function(M,i,r,Risk,delta,model_data){
    deltavec <- delta[Risk[[r]]]
    vec <- M[i,]
    j <- as.numeric(vec[1])
    k <- as.numeric(vec[2])
    C <- 0
    case <- 0
    if(j==k & k==ob[r]){
      case <- 1
      C <- -(exp(matrixStats::logSumExp(delta[ob[r]]-deltavec))-1)/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else if(j==k & k!=ob[r]){
      case <- 2
      deltavec1 <- deltavec[Risk[[r]]!=j & Risk[[r]]!=ob[r]]
      if(length(deltavec1)==0){
        C <- -(exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
      }
      if(length(deltavec1)!=0){
        C <- -((exp(matrixStats::logSumExp(2*delta[ob[r]]-delta[j]-deltavec1)))+exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
      }
    }
    else if(j==ob[r] & k!=ob[r]){
      case <- 3
      C <- (exp(delta[ob[r]]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else if(j!=ob[r] & k==ob[r]){
      case <- 4
      C <- (exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else{
      case <- 5
      C <- (exp(2*delta[ob[r]]-delta[j]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    C
  }
  #Second Step: Use ComputeC_ij and M matrix to define function that can get the full C_i matrix:
  #Here only these entries (j,k) of C_i where both j and k are in R_i are updated, because otherwise it will just be zeros.
  Compute_full_C_i <- function(M,r,model_data,Risk,delta){
    C_i <- matrix(data = 0,nrow = model_data$n, ncol = model_data$n)
    if(nrow(M)==0){
      C_i
    }
    else{
      values <- lapply(1:nrow(M),ComputeC_ij, M=M,r=r,Risk=Risk,model_data=model_data, delta=delta)
      for (i in 1:nrow(M)) {
        #Updating entries of C_i
        v <- M[i,]
        C_i[v[1],v[2]] <- values[[i]]
      }
      C_i
    }
  }
  #With these two functions, let's define M and use them to compute the full C matrix:
  listofC <- list(NULL)
  for (i in 1:length(ob)) {
    M <- as.matrix(expand.grid(Risk[[i]],Risk[[i]]))
    M <- M[M[,1]<=M[,2],,drop = FALSE]
    listofC[[i]] <- Compute_full_C_i(M,r=i,model_data = model_data,delta = delta,Risk = Risk)
  }
  #Sum all C_i's for observed times, and get the final C matrix.
  C_final <- Reduce('+',listofC)
  #Filling out the missing values of C using its symmetry.
  C_final <- forceSymmetric(C_final)
  #Removing all the entries that are differentiated with respect to Delta_11, as Delta_11 should be 
  # a constant zero. 
  MyC <- bdiag(C_final[2:nrow(C_final),2:ncol(C_final)],diag(rep(0,model_data$Wd-model_data$Nd),nrow =model_data$Wd-model_data$Nd))
  MyC <- as(MyC,"dgCMatrix")
  -MyC
}

#hessian_log_likelihood_alternative is an alternative to hessian_log_likelihood that improves the
#storage problem, as it only store one n*n C matrix each time instead of n of n*n C matrix. But has a 
#longer computation time. For a 100 data with approximately 40% censoring: (This is obsolete as well)

hessian_log_likelihood_alternative_old <- function(W,model_data) {
  Risk <- Get_all_Risk_Set(model_data)
  ob <- Get_Observed_index(model_data)
  delta <- prep_data_for_log_lik(W,model_data)
  #Temporally fill in Delta_11 into delta, to make the future computation easier.
  delta <- c(0,delta)
  #Compute the i-th observed time's Hessian matrix, C_i:
  #First Step: Define a function to compute a single entry of C_i.(r here is the index of Risk set that we are working on)
  #So we want to make sure that r-th risk set is matched with index of its corresponding observation, which is ob[r].
  ComputeC_ij <- function(M,i,r,Risk,delta,model_data){
    deltavec <- delta[Risk[[r]]]
    vec <- M[i,]
    j <- as.numeric(vec[1])
    k <- as.numeric(vec[2])
    C <- 0
    case <- 0
    if(j==k & k==ob[r]){
      case <- 1
      C <- -(exp(matrixStats::logSumExp(delta[ob[r]]-deltavec))-1)/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else if(j==k & k!=ob[r]){
      case <- 2
      deltavec1 <- deltavec[Risk[[r]]!=j & Risk[[r]]!=ob[r]]
      if(length(deltavec1)==0){
        C <- -(exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
      }
      if(length(deltavec1)!=0){
        C <- -((exp(matrixStats::logSumExp(2*delta[ob[r]]-delta[j]-deltavec1)))+exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
      }
    }
    else if(j==ob[r] & k!=ob[r]){
      case <- 3
      C <- (exp(delta[ob[r]]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else if(j!=ob[r] & k==ob[r]){
      case <- 4
      C <- (exp(delta[ob[r]]-delta[j]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    else{
      case <- 5
      C <- (exp(2*delta[ob[r]]-delta[j]-delta[k]))/((exp(matrixStats::logSumExp(delta[ob[r]]-deltavec)))^2)
    }
    C
  }
  #Second Step: Use ComputeC_ij and M matrix to define function that can get the full C_i matrix:
  #Here only these entries (j,k) of C_i where both j and k are in R_i are updated, because otherwise it will just be zeros.
  Compute_full_C_i <- function(M,r,model_data,Risk,delta){
    C_i <- matrix(data = 0,nrow = model_data$n, ncol = model_data$n)
    if(nrow(M)==0){
      C_i
    }
    else{
      values <- lapply(1:nrow(M),ComputeC_ij, M=M,r=r,Risk=Risk,model_data=model_data, delta=delta)
      for (i in 1:nrow(M)) {
        #Updating entries of C_i
        v <- M[i,]
        C_i[v[1],v[2]] <- values[[i]]
      }
      C_i
    }
  }
  #With these two functions, let's define M and use them to compute the full C matrix:
  #Sum all C_i's for observed times, and get the final C matrix:
  C <- 0
  for (i in 1:length(ob)) {
    M <- as.matrix(expand.grid(Risk[[i]],Risk[[i]]))
    M <- M[M[,1]<=M[,2],,drop = FALSE]
    C <- Reduce('+',list(C,Compute_full_C_i(M,r=i,model_data = model_data,delta = delta,Risk = Risk)))
  }
  #Filling out the missing values of C using its symmetry.
  C_final <- forceSymmetric(C)
  #Removing all the entries that are differentiated with respect to Delta_11, as Delta_11 should be 
  # a constant zero. 
  MyC <- bdiag(C_final[2:nrow(C_final),2:ncol(C_final)],diag(rep(0,model_data$Wd-model_data$Nd),nrow =model_data$Wd-model_data$Nd))
  MyC <- as(MyC,"dgCMatrix")
  MyC
}

#hessian_log_likelihood if no ties occured in your dataset.Every fast.(Use this one if you are sure that there are no ties)
hessian_log_likelihood <- function(W,model_data) {
  ob <- Get_Observed_index(model_data)
  delta <- prep_data_for_log_lik(W,model_data)
  #Temporally fill in Delta_11 into delta, to make the future computation easier.
  delta <- c(0,delta)
  #Define a function to compute one C_i matrix:
  GetCi <- function(i,model_data=model_data,delta=delta){
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
  R <- lapply(ob, GetCi, model_data=model_data,delta=delta)
  C_final <- Reduce('+',R)
  C_final <- bdiag(C_final[2:nrow(C_final),2:ncol(C_final)],diag(rep(0,model_data$Wd-model_data$Nd),nrow =model_data$Wd-model_data$Nd))
  C_final <- as(C_final,"dgCMatrix")
  return(-forceSymmetric(C_final))
}





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
    } 
    else {
      # warning("random walk-only models are rank deficient for case crossover. adding fudge factor. consider adding a linear term.")
      mat <- Q_matrix_rw2(theta,model_data,tau)
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
log_posterior_theta <- function(theta,W,model_data,Q = NULL) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
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

# Function to compute the log posterior for sigma, the standard deviation
# sigma = exp(-.5 * theta)
# jacobian is 2/sigma
log_posterior_sigma <- function(sigma,W,model_data,Q = NULL) {
  log(2/sigma) + log_posterior_theta(theta = -2 * log(sigma),
                                     W = W,
                                     model_data = model_data,
                                     Q = Q)
}







### OPTIMIZATION ###
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

#Using ipopt:
optimize_latentfield_ipopt <- function(theta,model_data,startingvals=NULL,Q = NULL,ipcontrol = NULL,boxconstr = c(-Inf,Inf)) {
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
  optfunchess_full <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
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

#Using Newton-Raphson:
optimize_latentfield_newton <- function(theta,model_data,startingvals=NULL,random_start_sd = .1,Q = NULL,maxiter = 100,eps = 1e-04,debug = FALSE) {
  # Get the Q matrix if not provided
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  if (length(startingvals) != model_data$Wd) stop(stringr::str_c("Length of starting values: ",
                                                                 length(startingvals),
                                                                 " but length of latent field: ",
                                                                 model_data$Wd))
  
  optfunc <- function(W) log_posterior_W(W,theta,model_data,Q)
  optfuncgrad <- function(W) grad_log_posterior_W(W,theta,model_data,Q)
  
  optfunchess <- function(W) -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  
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

optimize_all_thetas <- function(theta,model_data,num_starts,startingvals=NULL,random_start_sd = .1,starting_opt = "trust.optim",inner_opt = "trust.optim",debug = FALSE,doitinparallel = TRUE,...) {
  # theta is a LIST of theta vectors.
  # num_starts is the number of times to run the initial optimization with random starting values
  # Then, the best value from all these runs is taken and used as the starting value for the next value
  # of theta, which is only run once. Then this solution is used for the NEXT theta, and the optimization
  # proceeds in this way in series.
  
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  # Do the first optimization
  if (starting_opt == "ipoptr") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 5
    do_initial_opt <- function(x) optimize_latentfield_ipopt(theta[[1]],
                                                             model_data,
                                                             startingvals,
                                                             print_level = howmuchtoprint)
  }
  else if (starting_opt == "trust.optim") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 1
    do_initial_opt <- function(x) optimize_latentfield_trustoptim(theta[[1]],
                                                                  model_data,
                                                                  startingvals,
                                                                  report_freq = howmuchtoprint)
  }
  else if (starting_opt == "newton") {
    warning("I wouldn't recommend using newton as the starting optimizer. Try trust.optim instead")
    do_initial_opt <- function(x) optimize_latentfield_newton(theta[[1]],
                                                              model_data,
                                                              startingvals,
                                                              debug = debug)
  }
  else {
    stop(stringr::str_c("Unknown initial optimizer: ",optimizer))
  }
  cat("Performing initial optimization...\n")
  
  tm <- proc.time()
  if (num_starts > 1 & doitinparallel) {
    initial_opt <- mclapply(1:num_starts,do_initial_opt, mc.cores = cores) %>% optlist_to_tibble()
  }
  else {
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
                                                                            print_level = howmuchtoprint
    )
  }
  else if (inner_opt == "trust.optim") {
    do_inner_opt <- function(theta,startingvals) optimize_latentfield_trustoptim(theta,
                                                                                 model_data,
                                                                                 startingvals,
                                                                                 report_freq = howmuchtoprint)
  }
  else if (inner_opt == "newton") {
    do_inner_opt <- function(theta,startingvals) optimize_latentfield_newton(theta,
                                                                             model_data,
                                                                             startingvals,
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

# Now do a version that just does them all in parallel
optimize_all_thetas_parallel <- function(theta,model_data,startingvals = NULL,random_start_sd = 0,optimizer = "trust.optim",optcontrol = NULL,debug = FALSE,doparallel = TRUE) {
  # Trick: to use 0 as a starting value, set random_start_sd = 0. rnorm(1,0,sd = 0) returns 0 wp1
  # Startingvals
  if (is.null(startingvals)) startingvals <- rnorm(model_data$Wd,sd = random_start_sd)
  if (optimizer == "ipopt") {
    howmuchtoprint <- 0
    if (debug) howmuchtoprint <- 5
    do_opt <- function(theta) {
      optimize_latentfield_ipopt(theta,
                                 model_data,
                                 startingvals,
                                 rnorm(model_data$Wd,sd = random_start_sd),
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
                                      trcontrol = optcontrol)
    }
  }
  cat("Performing optimization, start time: ", format(Sys.time(),"%H:%M:%S"),"\n")
  tm <- proc.time()
  if (doparallel) {
    opt <- mclapply(theta,do_opt, mc.cores = cores)
  } 
  else {
    opt <- lapply(theta,do_opt)
  }
  opt_time <- proc.time() - tm
  opt_time <- unname(opt_time["elapsed"])
  cat("Time taken for optimization with ",length(theta)," values of theta:",opt_time,"seconds.\n")
  # Take the BEST mode for each theta
  opt %>% optlist_to_tibble() %>% group_by(theta) %>% arrange(function_value) %>% slice(1)
}








### COMPUTATION OF FINAL RESULT AND SUMMARIES ###
# Take in output of optimize_all_thetas(). Add the value of the log_posterior_theta().
# This will be used for posterior summary and density computation; but we also need to plot it directly.
add_log_posterior_values <- function(optresults,model_data) {
  optresults <- ungroup(optresults)
  optresults$sigma <- exp(-.5 * optresults$theta)
  # Log posterior for theta
  logposttheta <- c()
  for (i in 1:length(optresults$theta)) {
     logposttheta[i] <- log_posterior_theta(optresults$theta[i],optresults$solution[[i]],model_data)
     }
  optresults$theta_logposterior <- logposttheta
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
  # Compute values of the Gaussian approximation for each theta, for each value
  # of W supplied. 
  # Marginal variances
  precision_matrices <- model_results %>% 
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data),
      theta = ..1)
    )
  QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]] + exp(correction) * Diagonal(dim(.x[["Q"]])[1]),theta = .x[["theta"]]))
  margvar <- QpC %>%
    purrr::map(~solve(.x[["QpC"]])) %>%
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

# Compute marginal means and variances only
compute_marginal_means <- function(i,model_results,model_data,constrA = NULL) {
  # i: vector of indices. which elements of the latent field do you want the marginal for?
  # model_results: tibble containing output of optimization
  # Have to deal with the fact that the means can be negative. This function basically does
  # fancy averaging.
  # constrA is a matrix of linear constraints, AW = 0 where W is the latent field.
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
  } 
  else {
    # Correct for linear constraints
    # List of uncorrected means, one for each theta
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))
    # Compute the Q^-1A^T term. Q here is really Q + C of course...
    precision_matrices <- model_results %>% 
      purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                        theta = ..1)
      )
    hessians <- model_results %>% 
      purrr::pmap(~list(
        C = hessian_log_likelihood(W = ..4,model_data = model_data),
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
  postvals <- exp(model_results$theta_logposterior + log(intweights))
  sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
}


#Currently, this is not efficient enough, not compatible with super large dataset as it is 
#computing the inverse of a large dense matrix by brutal force.
compute_marginal_variances_without_RandomEffect <- function(i,model_results,model_data,constrA = NULL,correction = -Inf) {
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
    purrr::pmap(~list(Q = as(Q_matrix(theta = ..1,model_data = model_data),Class = "dsCMatrix"),
                      theta = ..1)
    )
  Q <- precision_matrices[[1]]$Q
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data),
      theta = ..1)
    )
  C <- hessians[[1]]$C
  # Marginal variances: add the precision and the hessian and get diagOfInv
  margvar <- diag(diag(solve(C+Q))[i])
  if (nrow(model_results) == 1) return(margvar)
  postvals <- exp(model_results$theta_logposterior)
  sweep(margvar,1,postvals,"*") %>% apply(2,sum)
}



compute_marginal_variances <- function(i,model_results,model_data,constrA = NULL,correction = -Inf) {
  # W: vector of values of W to compute log pi(W_i|y) at
  # i: index. which element of the latent field do you want the marginal for?
  # model_results: tibble containing output of optimization
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
  hessians <- model_results %>% 
    purrr::pmap(~list(
      C = hessian_log_likelihood(W = ..4,model_data = model_data),
      theta = ..1)
    )
  # Marginal variances: add the precision and the hessian and get diagOfInv
  margvar <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
    purrr::map(~diag(solve(.x))[i]) %>%
    purrr::reduce(rbind)
  if (nrow(model_results) == 1) return(margvar)
  postvals <- exp(model_results$theta_logposterior + log(intweights))
  sweep(margvar,1,postvals,"*") %>% apply(2,sum)
}



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
  hessians <- list()
  myhes <- mclapply(model_results$solution, hessian_log_likelihood,model_data = model_data,mc.cores = detectCores())
  for (i in 1:length(length(model_results$theta))) {
    hessians[[i]] <- list(C=myhes[[i]],theta=model_results$theta[i])
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
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[,j],"sparseVector")
      AA %>% purrr::map(~compute_var_one_lincomb(.x,Q)) %>% purrr::reduce(c)
    }
    compute_one_lincomb_correction <- function(a,WW,VV) {
      as.numeric(crossprod(crossprod(VV,a),crossprod(WW,a)))
    }
    compute_all_lincomb_correction <- function(A,WW,VV) {
      AA <- list()
      for (j in 1:ncol(lincomb)) AA[[j]] <- as(lincomb[,j],"sparseVector")
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
    corrected_hessians <- model_results %>%
      purrr::pmap(~list(
        C = hessian_log_likelihood(W = ..12,model_data = model_data),
        theta = ..1)
      )
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


