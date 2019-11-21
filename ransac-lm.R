select_inliers <- function(min_points, num_rows) {
  # draw a vector of observations and return it as a logical index vector
  sample <- seq_len(num_rows) %in% sample(num_rows, min_points)
  sample
}

find_model <- function(formula, data) {
  model <- lm(formula, data)
  model
}

find_inliers <- function(model, formula, data, error_threshold) {
  # find elements with small absolute error of prediction
  response <- all.vars(as.formula(formula))[[1]]
  predictions <- predict(model, data[, -which(names(data) == response)])
  inliers <- as.logical(abs(data[, response] - predictions) <= error_threshold)
  inliers
}

extract_error <- function(model) {
  # mean square error
  mean(model[["residuals"]]^2)
}

check_inputs <- function(formula, data, error_threshold, inlier_threshold,
                         iterations, seed) {
  # lm can use formulas as well as objects coercible to a formula
  tryCatch(as.formula(formula),
    error = function(e) stop("Couldn't coerce given 
                                    formula to formula object")
  )
  checkmate::assert_data_frame(
    data,
    types = c(
      "logical", "integer", "integerish", "double", "numeric",
      "complex", "character", "factor"
    ),
    any.missing = FALSE, null.ok = FALSE
  )
  checkmate::assert_number(error_threshold, lower = 0)
  checkmate::assert_count(inlier_threshold)
  checkmate::assert_count(iterations)
  checkmate::assert_count(seed)
  response <- all.vars(formula)[[1]]
  covariates <- attr(terms.formula(formula, data = data), "term.labels")
  # make sure all the terms in the formula are columns in the dataframe
  checkmate::assert_true(all(c(covariates, response) %in% colnames(data)))
  # make sure the linear model isn't underspecified
  checkmate::assert_true(length(covariates) < nrow(data))
}

#' Fit a linear model with RANSAC, see
#' https://en.wikipedia.org/wiki/Random_sample_consensus
#'
#' @param formula - formula object or string coercible to formula
#' @param data - dataframe containing the variables from the formula
#' @param error_threshold - maximum absolute error for finding inliers
#' @param inlier_threshold - minimal number of inliers required for a model
#' @param iterations - number of iterations
#' @param seed - seed
#'
#' @return List containing the best linear model and data with consensus set
#' @export
ransaclm <- function(formula, data, error_threshold,
                     inlier_threshold, iterations = 1000, seed) {
  check_inputs(formula, data, error_threshold, 
               inlier_threshold, iterations, seed)
  set.seed(seed)
  best_fit <- NULL
  best_error <- Inf
  consensus_set <- NULL
  # extract the number of covariates (can deal with ~ . type formulas)
  min_points <- length(attr(
    terms.formula(formula, data = data), "term.labels")) + 1

  for (i in seq_len(iterations)) {

    # warnings can occur if the data is rank deficient - ignore those samples
    tryCatch({
        current_sample <- select_inliers(min_points, nrow(data))
        current_model <- find_model(formula, data[current_sample, ])
        inliers <- find_inliers(current_model, formula, data, error_threshold)
      },
      warning = function(w) {
        inliers <- 0
      }
    )

    if (sum(inliers) >= inlier_threshold) {
      better_model <- find_model(formula, data[(current_sample | inliers), ])
      current_error <- extract_error(better_model)

      # save best model and its parameters
      if (current_error < best_error) {
        best_fit <- better_model
        best_error <- current_error
        consensus_set <- inliers | current_sample
      }
    }
  }
  data[[".consensus_set"]] <- consensus_set
  list(model = best_fit, data = data)
}


#' Parallelized version of ransaclm
#'
#' @param formula - formula object or string coercible to formula
#' @param data - dataframe containing the variables from the formula
#' @param error_threshold - maximum absolute error for finding inliers
#' @param inlier_threshold - minimal number of inliers required for a model
#' @param iterations - number of iterations
#' @param seed - seed
#'
#' @return List containing the best linear model and data with consensus set
#' @export
ransaclm_par <- function(formula, data, error_threshold, 
                         inlier_threshold, iterations = 1000, seed) {
  check_inputs(formula, data, error_threshold, 
               inlier_threshold, iterations, seed)
  set.seed(seed)
  # extract the number of covariates (can deal with ~ . type formulas)
  min_points <- length(attr(
    terms.formula(formula, data = data), "term.labels")) + 1

  single_ransac_loop <- function() {

    # warnings can occur if the data is rank deficient - ignore those samples
    tryCatch({
        current_sample <- select_inliers(min_points, nrow(data))
        current_model <- find_model(formula, data[current_sample, ])
        inliers <- find_inliers(current_model, formula, data, error_threshold)
      },
      warning = function(w) {
        inliers <- 0
      }
    )

    # early stopping if model doesn't produce enough inliers
    if (sum(inliers) < inlier_threshold) {
      return(list(model = NULL, error = Inf, consensus = NULL))
    }

    # save model, error and consensus set for good models
    current_model <- find_model(formula, data[(current_sample | inliers), ])
    current_error <- extract_error(current_model)
    consensus_set <- inliers | current_sample
    return(list(model = current_model, 
                error = current_error, 
                consensus = consensus_set))
  }

  # setup the cluster
  local_cluster <- parallel::makeCluster(3)
  needed_vars <- c(
    "select_inliers", "min_points", "data", "find_model",
    "formula", "error_threshold", "inlier_threshold",
    "find_inliers", "extract_error"
  )
  parallel::clusterExport(
    cl = local_cluster, varlist = needed_vars,
    envir = environment()
  )

  # obtain a list of results in parallelized fashion
  result <- parallel::parSapply(
    cl = local_cluster, seq_len(iterations),
    function(i) {
      single_ransac_loop()
    }
  )
  parallel::stopCluster(local_cluster)

  # we can't compare the models while doing parallel computations - 
  # evaluate results and find best model afterwards
  best_model_idx <- which(unlist(result[2, ] == min(unlist(result[2, ]))))[[1]]
  data[[".consensus_set"]] <- result[, best_model_idx][["consensus"]]
  list(model = result[, best_model_idx][["model"]], data = data)
}
