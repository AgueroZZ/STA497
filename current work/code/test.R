library(parallel)

lis

#$n
#[1] 4
#$W
#[1] 2 3 4 1 1 1 1
#$ID
#[1] 1 2 3 4
#$times
#[1] 3.0 3.2 1.9 2.0
#$entry
#[1] 0 0 0 0
#$censoring
#[1] 1 1 1 0

grad_log_likelihood(c(lis$W,1,3,2,5,1),lis)

#0.23166014 -0.95862930  0.01521943  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000

hessian_log_likelihood(c(lis$W,1,3,2,5,1),lis)

#3 x 3 Matrix of class "dsyMatrix"




