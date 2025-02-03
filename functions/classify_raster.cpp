// Generated with ChatGPT. Modified as necessary.
#include <Rcpp.h>
#include <cmath> // For std::isnan

// [[Rcpp::export]]
Rcpp::NumericMatrix applyThreshold(Rcpp::NumericMatrix matrix, Rcpp::NumericVector thresholds) {
    int nrows = matrix.nrow();
    int ncols = matrix.ncol();

    // Check if the number of thresholds matches the number of columns
    if (thresholds.size() != ncols) {
        Rcpp::stop("Length of thresholds vector must match the number of columns in the matrix.");
    }

    // Iterate over each element of the matrix
    for (int col = 0; col < ncols; ++col) {
        for (int row = 0; row < nrows; ++row) {
            double& value = matrix(row, col);
            if (std::isnan(value)) {
                // If the value is NaN (NA in R), leave it unchanged
                continue;
            } else if (value < thresholds[col]) {
                // If the value is less than the threshold, set it to 0
                value = 0.0;
            }
            // Otherwise, retain the original value
        }
    }

    return matrix;
}
