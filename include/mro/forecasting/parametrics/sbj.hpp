#pragma once

/**
 * @file sbj.hpp
 * @brief Shale-Boylan-Johnston (SBJ) intermittent demand forecasting method.
 *
 * Reference:
 *   Shale, E.A., Boylan, J.E. & Johnston, F.R. (2006). Forecasting for
 *   intermittent demand: the estimation of an unbiased average. Journal
 *   of the Operational Research Society, 57(5), 588-592.
 */

#include <vector>

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/types/initialization_results.hpp"

namespace mro::forecasting::parametrics {

	/**
	 * @brief Output of the SBJ forecasting method.
	 *
	 * SBJ is, like SBA, a bias-corrected variant of Croston's method. It
	 * uses the same underlying z_hat, p_hat and sigma2_hat from Croston's
	 * smoothing, but applies the 'exact' deflating factor
	 * (1 - alpha/(2-alpha)) derived by Shale, Boylan & Johnston (2006),
	 * instead of the approximate factor used by SBA.
	 */
	struct SbjResult : public mro::forecasting::types::ForecastResult {
		double z_hat;       ///< Smoothed estimate of demand size (from Croston).
		double p_hat;       ///< Smoothed estimate of demand interval (from Croston).
		double forecast;    ///< Bias-corrected forecast = (1 - alpha/(2-alpha)) * Croston forecast.
		double sigma2_hat;  ///< Sample variance of the non-zero demand sizes observed.
		double alpha;       ///< Smoothing constant used to obtain this result.

		/// @copydoc ForecastResult::point_forecast
		double point_forecast() const override { return forecast; }

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * The bias-correction factor (1 - alpha/(2-alpha)) is squared,
		 * since Var(c*X) = c^2 * Var(X), and applied to Croston's
		 * underlying demand variance (see CrostonResult::demand_variance()).
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Forecasts intermittent demand using the Shale-Boylan-Johnston
	 *        (SBJ) 'exact' bias correction.
	 *
	 * Runs Croston's method internally and deflates both the forecast and
	 * its variance by (1 - alpha/(2-alpha)).
	 *
	 * @param time_series  Chronological sequence of period demands (zeros allowed).
	 * @param init         Initial values for z_hat and p_hat, from any
	 *                     initialization method (e.g. calculate_by_moments()).
	 * @param alpha        Smoothing constant, 0 < alpha <= 1.
	 * @param skip_periods Number of periods at the beginning of the series to skip
	 *                     (typically the periods used during initialization).
	 * @return SbjResult with the bias-corrected forecast and variance.
	 * @throws std::invalid_argument if @p alpha is not in (0, 1], @p init.p1 <= 0, or
	 *         @p skip_periods is greater than or equal to @p time_series.size().
	 */
	SbjResult forecast_by_sbj(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			size_t skip_periods = 0);

} // namespace mro::forecasting::parametrics
