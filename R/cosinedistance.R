#' Pairwise cosine similarity
#'
#' Compare each pair of rows
#'
#' @param x Matrix
#' @return A symmetric matrix of cosine similarity values.
cossimX <- function(x) {
  y <- x / sqrt(rowSums(x * x))
  return(y %*% t(y))
}

#' Cosine similarity between a pair of matrices
#'
#' Matrices \code{x} and \code{y} must be of equal dimension. The i-th row in
#' \code{x} is compared to the i-th row in \code{y}.
#'
#' @param x Matrix
#' @param y Matrix
#' @return A similarity value for each row.
cossimXYMatMat <- function(x, y) {
  assertthat::are_equal(dim(x), dim(y))
  a <- rowSums(x^2)
  b <- rowSums(y^2)
  e <- rowSums(x * y)
  return(e / sqrt(a * b))
}

#' Cosine similarity between a pair of vectors
#'
#' Vectors \code{x} and \code{y} must be equal length
#'
#' @param x Vector
#' @param y Vector
#' @return A single similarity value
cossimXYVecVec <- function(x, y) {
  assertthat::are_equal(length(x), length(y))
  a <- sum(x^2)
  b <- sum(y^2)
  e <- sum(x * y)
  return(e / sqrt(a * b))
}

#' Cosine similarity between multiple rows and single vector
#'
#' @param x Matrix
#' @param y Vector
#' @return A similarity value for each row in \code{x}
cossimXYMatVec <- function(x, y) {
  assertthat::are_equal(ncol(x), length(y))
  assertthat::are_equal(nrow(x) %% length(y), 0)
  a <- rowSums(x^2)
  b <- sum(y^2)
  e <- colSums(t(x) * y)
  return(e / sqrt(a * b))
}

#' Compute cosine similarity
#'
#' Cosine similarity is the cosine of the angle between the two rays that each
#' pass through the origin and another point. Cosine is computed as the inner
#' product of the two vectors divided by the product of their norms.
#'
#' @param x A vector or matrix.
#' @param y An optional second vector or matrix
#' @return The cosine similarity
#'
#' @export
cosine_similarity <- function(x, y = NULL) {
  if (is.null(y))
    return(cossimX(x))

  if (is.matrix(x) && is.matrix(y))
    return(cossimXYMatMat(x, y))

  if (is.matrix(x) && is.vector(y))
    return(cossimXYMatVec(x, y))

  if (is.vector(x) && is.matrix(y))
    return(cossimXYMatVec(y, x))

  if (is.vector(x) && is.vector(y))
    return(cossimXYVecVec(x, y))
}

#' Compute cosine distances
#'
#' Cosine distance is simply 1 - cosine_similarity
#'
#' @param x A vector or matrix.
#' @param y An optional second vector or matrix
#' @return The cosine distance
#'
#' @export
cosine_distance <- function(x, y = NULL) {
  return(1 - cosine_similarity(x, y))
}
