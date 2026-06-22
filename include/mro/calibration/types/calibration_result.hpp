#pragma once

/**
 * @file calibration_result.hpp
 * @brief Output types shared by every damping-constant calibration method.
 *
 * Every calibration method (MASE-, MAD- and MSE-minimizing) returns one
 * of these two structs: SingleParamCalibResult for forecasting methods
 * driven by a single damping constant (Croston, SBA, SBJ — alpha only),
 * and DualParamCalibResult for methods driven by two (TSB — alpha and
 * beta). The forecasting methods themselves are unchanged: they still
 * take alpha/beta as plain double parameters. Calibration is just a
 * search procedure that finds good values for them; the caller plugs
 * the resulting struct's alpha (and beta, if applicable) back into the
 * same forecast_by_*() call they would have used anyway.
 */

namespace mro::calibration::types {

	/**
	 * @brief Calibrated damping constant for single-parameter forecasting
	 *        methods (Croston, SBA, SBJ).
	 */
	struct SingleParamCalibResult {
		double alpha = 0.0;          ///< Calibrated damping constant.
		double minimum_error = 0.0;  ///< Value of the statistical metric (MASE, MAD or MSE) attained at alpha.
	};

	/**
	 * @brief Calibrated damping constants for dual-parameter forecasting
	 *        methods (TSB).
	 */
	struct DualParamCalibResult {
		double alpha = 0.0;          ///< Calibrated occurrence/size-smoothing constant.
		double beta = 0.0;           ///< Calibrated probability-smoothing constant.
		double minimum_error = 0.0;  ///< Value of the statistical metric (MASE, MAD or MSE) attained at (alpha, beta).
	};

} // namespace mro::calibration::types
