#' Perform the cyclic shift procedure on the columns of a matrix
#'
#'
#' @param X a matrix or a data frame of copy number data. The rows and columns
#'   of X correspond to genes and subjects, respectively.
#'
#' @param randomSeed a random seed.  Default = NULL.
#'
#' @return a matrix Z whose dimensions are the same as X.  Each column of Z
#'   is obtained by perform a cyclic shift of the corresponding column of X.
#'
#' @examples
#' test = matrix(c(1:50), 10, 5)
#' cyclicShiftColR(test, randomSeed = NULL)
#'
#' @export

cyclicShiftColR = function(X, randomSeed = NULL)
	{
	n = nrow(X)
	m = ncol(X)
	if (!is.null(randomSeed))
		{
		set.seed(randomSeed)
		}
	Z = X
	shiftVec = sample(1:n, m, replace = TRUE)
	for (i in 1:m)
		{
		k = shiftVec[i]
		Z[,i] = ((k == 1) * Z[,i]) + ((k != 1) * Z[c(c(k:n), c(1:(k - 1))), i][c(1:n)])
		}
	return(Z)	
	}

