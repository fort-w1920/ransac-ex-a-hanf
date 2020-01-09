#' Fit a random sample consensus (RANSAC) linear model for outlier-contaminated data
#     Details: https://en.wikipedia.org/wiki/Random_sample_consensus
#' input: formula: valid lm formula (univariate response only, no offsets)
#'        data: data.frame containing the variables in formula
#'        error_threshold: observations with absolute prediction error below
#           this are considered for the consensus set
#'        inlier_threshold: minimum size of consensus set
#'        iterations: how many random draws of candidate sets to try (
#'          defaults to half the number of observations in data)
#'        seed: RNG seed
#' output: list with model (best fit lm-object on consensus set),
#'   data (modified input data (complete cases), with additional boolean variable
#'   ".consensus_set" indicating the best set that was found)
#' imports checkmate, future.lapply. use future::plan("multiprocess") for parallelizing.
ransaclm <- function(formula, data, error_threshold, inlier_threshold,
                     iterations = nrow(data) / 2, seed = NULL) {
  checkmate::assert_class(formula, "formula")
  checkmate::assert_class(data, "data.frame")
  checkmate::assert_number(error_threshold, lower = 0)
  checkmate::assert_integerish(inlier_threshold, lower = 1, upper = nrow(data) - 1)
  checkmate::assert_count(iterations, positive = TRUE)
  checkmate::assert_int(seed, null.ok = TRUE, lower = 1)

  if (!is.null(seed)) set.seed(as.integer(seed))

  model_all <- try(lm(formula, data), silent = TRUE)
  if (inherits(model_all, "try-error")) {
    stop("Could not fit model on the supplied data, `lm` exited with error:\n",
         model_all)
  }
  # multivariate response not allowed, easiest way to check is to disallow
  #   matrix-coefficients
  if (is.matrix(coefficients(model_all))) {
    stop("Multivariate response models are not implemented.")
  }
  if (!is.null(model_all[["offset"]])) {
    stop("Models with offset are not implemented.")
  }

  # model.frame gets rid of incomplete observations in <data>
  data_all <- model.frame(model_all)
  design <- model.matrix(model_all)
  response <- model.response(data_all)

  # do the work:
  candidates <- future.apply::future_replicate(
    n = iterations,
    expr = ransac_once(design, response, error_threshold, inlier_threshold),
    simplify = FALSE
  )
  best_error <- which.min(vapply(candidates, `[[`, "error",
                                 FUN.VALUE = numeric(1)))
  if (!length(best_error)) {
    warning("No RANSAC model satisfied criteria, try lower inlier or error tresholds.\n")
    return(list(model = NULL, data = data_all))
  }
  best_set <- candidates[[best_error]][["set"]]

  # use update.lm so return object is a real lm-object (not a list from .lm.fit)
  best_model <- update(model_all, data = data_all[best_set, ])
  data_all$.consensus_set <- seq_len(nrow(data_all)) %in% best_set
  list(model = best_model, data = data_all)
}

ransac_once <- function(design, response, error_threshold, inlier_threshold) {
  # 1. generate candidate sets
  candidate_set <- get_candidate_set(design, response, error_threshold)
  # 2. filter out too small candidate sets
  keep_set <- length(candidate_set) > inlier_threshold
  if (!keep_set) {
    return(list(error = NA_real_, set = numeric(0)))
  }
  # 3. refit & return error on candidate set
  error <- refit_candidate_set(candidate_set,
                               design = design, response = response)
  list(error = error, set = candidate_set)
}

# return indices of all observations which <design>%*%<coefficients> predicts well
get_inliers <- function(coefficients, design, response, error_threshold) {
  predictions <- design %*% coefficients
  errors <- response - predictions
  which(abs(errors) < error_threshold)
}

# candidate for the consensus set based on a random draw of a minimal dataset:
get_candidate_set <- function(design, response, error_threshold) {
  # draw row indices of candidate observations for minimal dataset
  use <- sample(x = seq_len(nrow(design)), size = ncol(design), replace = FALSE)
  # use .lm.fit instead of lm for better performance
  candidate_model <- try(.lm.fit(x = design[use, ], y = response[use]))
  if (inherits(candidate_model, "try-error") ||
      any(is.na(candidate_model[["coefficients"]]))) {
    warning("Model fit failed on subsample.")
    return(numeric(0))
  }
  # get set of observations which this_model predicts well
  get_inliers(candidate_model[["coefficients"]],
              design = design, response = response,
              error_threshold = error_threshold)
}

# get error for candidate consensus set
refit_candidate_set <- function(candidate_set, design, response) {
  # use .lm.fit instead of lm for better performance
  candidate_set_model <- try(.lm.fit(
    x = design[candidate_set, ],
    y = response[candidate_set]
  ))
  if (inherits(candidate_set_model, "try-error")) {
    warning("Model fit failed on candidate set.")
    return(NA_real_)
  }
  mean(residuals(candidate_set_model)^2)
}
