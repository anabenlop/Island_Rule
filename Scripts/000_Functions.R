# load libraries
library(corpcor)
library(Matrix)

# this function is to make sure that a is matrix positive-definitive
# it uses the function nearPD in the Matrix package to compute the
# nearest positive definite matrix to an approximate one
# Function reused from: https://osf.io/ayxrt/
PDfunc <- function(matrix){
  require(Matrix)
  if(corpcor::is.positive.definite(matrix) == T){
    x <- corDiag
  }else{
    x <- Matrix::nearPD(matrix)
  }
  mat <- x$mat
  return(mat)
}

# Variance-covariance matrix for studies with multiple treatments-common control
# for response ratios
calc.v <- function(x) {
  v <- matrix(x$sd_m[1]^2 / (x$N_m[1]*x$Mean_m[1]^2), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$var
  v
}

# R2 - Nakagawa & Schielzeth 2013

mR2.func <- function(model){
  fix<-var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  R2m<-fix/(fix+sum(model$sigma2))
  x <- round(100*R2m, 2)
  return(x)
}


cR2.func <- function(model){
  fix<-var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  R2c<-(sum(model$sigma2) - model$sigma2[2])/(fix+sum(model$sigma2))
  x <- round(100*R2c, 2)
  return(x)
}

