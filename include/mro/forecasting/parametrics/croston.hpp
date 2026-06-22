#pragma once

/**
 * @file croston.hpp
 * @brief Croston (1972) intermittent demand forecasting method.
 *
 * Reference:
 *   Croston, J.D. (1972). Forecasting and stock control for intermittent
 *   demands. Operational Research Quarterly, 23(3), 289-303.
 */

#include <vector>

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/types/initialization_results.hpp"

namespace mro::forecasting::parametrics {

	/**
	 * @brief Output of the Croston forecasting method.
	 *
	 * Implements mro::forecasting::ForecastResult so it can be consumed
	 * generically by mro::metrics, alongside any other forecasting
	 * method's result (SBA, TSB, SBJ, ...).
	 */
	struct CrostonResult : public mro::forecasting::types::ForecastResult {
		double z_hat;       ///< Smoothed estimate of demand size.
		double p_hat;       ///< Smoothed estimate of demand interval.
		double forecast;    ///< Forecast per period = z_hat / p_hat.
		double sigma2_hat;  ///< Sample variance of the non-zero demand sizes observed.
		double alpha;       ///< Smoothing constant used to obtain this result.

		/// @copydoc ForecastResult::point_forecast
		double point_forecast() const override { return forecast; }

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Based on the variance of the underlying compound demand process
		 * (Bernoulli arrival with probability 1/p, Normal(mu, sigma^2)
		 * demand sizes), using z_hat as the estimate of mu, sigma2_hat as
		 * the estimate of sigma^2, and p_hat as the estimate of p:
		 *
		 *   Var = alpha/(2-alpha) * [ (p-1)/p^2 * mu^2 + sigma^2/p ]
		 *
		 * Note this is the variance of the per-period demand-generating
		 * process implied by the model, not the sample variance of the
		 * raw historical series (which conflates zero and non-zero
		 * periods and does not reflect the model's current state).
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Forecasts intermittent demand using Croston's (1972) method.
	 *
	 * Smoothing updates occur only on periods with non-zero demand:
	 *   z_hat is updated toward the new demand size,
	 *   p_hat is updated toward the interval since the last non-zero demand.
	 * On zero-demand periods, both estimates remain unchanged.
	 *
	 * Accepts any InitializationParameters-derived struct, so the
	 * initialization method (moments, simple average, etc.) is decoupled
	 * from the forecasting method.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param init        Initial values for z_hat and p_hat, from any
	 *                     initialization method (e.g. calculate_by_moments()).
	 * @param alpha       Smoothing constant, 0 < alpha <= 1.
	 * @return CrostonResult with the final z_hat, p_hat, and the resulting forecast.
	 * @param skip_periods Number of periods at the beginning of the series to skip
	 * (typically the periods used during initialization to prevent overfitting).
	 * @throws std::invalid_argument if @p alpha is not in (0, 1], @p init.p1 <= 0, or 
	 * @p skip_periods is greater than or equal to @p time_series.size().
	 */
	CrostonResult forecast_by_croston(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			size_t skip_periods = 0);

} // namespace mro::forecasting::parametrics
