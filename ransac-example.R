# immer set.seed() um Ergebnisse reproduzierbar zu machen...
set.seed(1874374111)
source("ransac-utils.R")
source("ransac-lm.R")

data_simple <- make_ransac_data(n_obs = 100, n_coef = 1, inlier_fraction = 0.7)

# univariate example:
ransac_simple <- ransaclm(y ~ . - inlier,
                          data = data_simple, error_threshold = 2,
                          inlier_threshold = 50, iterations = 100, seed = 20171111
)

validate_ransac(ransac_simple)


# parallel version:
ransac_simple <- ransaclm_par(y ~ . - inlier,
                           data = data_simple, error_threshold = 2,
                           inlier_threshold = 50, iterations = 100, seed = 20171111
)

validate_ransac(ransac_simple)