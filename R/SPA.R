#' @name SPA
#' @title Successive projection algorithm for separable NMF
#'
#' @description At each step of the algorithm, the column of R maximizing ||.||_2 is
#' extracted, and R is updated by projecting its columns onto the orthogonal
#' complement of the extracted column. The residual R is initializd with X.
#'
#' See N. Gillis and S.A. Vavasis, Fast and Robust Recursive Algorithms
#' for Separable Nonnegative Matrix Factorization,  IEEE Trans. on Pattern
#' Analysis and Machine Intelligence 36 (4), pp. 698-714, 2014.
#'
#' @param X an m-by-n matrix.
#' Ideally admitting a near-separable factorization, that is, X = WH + N, where
#'
#' conv(\[W, 0]) has r vertices,
#'
#' H = \[I,H']P where I is the identity matrix, H'>= 0 and its
#' columns sum to at most one, P is a permutation matrix, and
#'
#' N is sufficiently small.
#'
#' @param r number of columns to be extracted.
#' @param options list of string
#' if options contain
#'
#' "normalize" := 1, will scale the columns of X so that they sum to one,
#'              hence matrix H will satisfy the assumption above for any
#'              nonnegative separable matrix X.
#'
#' "normalize" := 0, is the default value for which no scaling is
#'              performed. For example, in hyperspectral imaging, this
#'              assumption is already satisfied and normalization is not
#'              necessary.
#'
#' "precision" : stops the algorithm when
#'               max_k ||R(:,k)||_2 <= relerror * max_k ||X(:,k)||_2
#'               where R is the residual, that is, R = X-X(:,K)H.
#'               default: 1e-6
#'
#' "display"   : = 1, displays the iteration count (default)
#'
#' @return `K` index set of the extracted columns.
#'
#' @export
#'
#' @examples
#' SPA(matrix(1:16,4),2,c('normalize'))
#'
SPA <- function(X,r,options=list()){

  m <- dim(X)[1]
  n <- dim(X)[2]

  if ('normalize' %in% options){
    options.normalize <- 1
  }
  if (!'display' %in% options){
    options.display <- 1
  }
  if (!'precision' %in% options){
    options.precision <- 1e-6
  }
  if (options.normalize == 1){
    # Normalization of the columns of M so that they sum to one
    D <- diag(c(colSums(X)^(-1)))
    X <- X %*% D
  }

  normX0 <- colSums(X**2)
  nXmax <- max(normX0)
  normR <- normX0

  i <- 1
  # Perform r recursion steps (unless the relative approximation error is
  # smaller than 10^-9)
  if (options.display == 1){
    cat('Extraction of the indices by SPA: \n')
  }
  K <- list()
  while (i <= r && sqrt(max(normR)/nXmax) > options.precision){
    # Select the column of M with largest l2-norm
    a <- max(normR)
    # Norm of the columns of the input matrix M
    # Check ties up to 1e-6 precision

    ckeck_precision <- ((a - normR) / a <= 1e-6)

    # In case of a tie, select column with largest norm of the input matrix M
    if (sum(ckeck_precision) > 1){
      b <- which.max(normX0[which(ckeck_precision)])
    }else{
      b <- which(ckeck_precision)[1]
    }
    # Update the index set, and extracted column
    K[i] <- b
    if (i == 1){
      U <- matrix(X[,b])
    }else{
      U <- cbind(U, matrix(X[,b]))
    }
    # Compute (I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T) U(:,i), that is,
    # R^(i)(:,J(i)), where R^(i) is the ith residual (with R^(1) = M).
    j <- 1
    while(j <= i-1){
      U[,i] <- U[,i] - U[,j] %*% (t(U[,j]) %*% U[,i])
      j <- j + 1
    }
    # Normalize U(:,i)
    U[,i] <- U[,i] / norm(U[,i], type='2')
    # Update the norm of the columns of M after orhogonal projection using
    # the formula ||r^(i)_k||^2 = ||r^(i-1)_k||^2 - ( U(:,i)^T m_k )^2 for all k.
    normR <- normR - (t(U[,i]) %*% X)^2
    if (options.display == 1){
      sprintf('%2.0f...', i)
      if (i %% 10 == 0){
        cat('\n')
      }
    }
    i = i + 1
  }
  if (options.display == 1){
    cat("\n")
  }
  return(K)
}
