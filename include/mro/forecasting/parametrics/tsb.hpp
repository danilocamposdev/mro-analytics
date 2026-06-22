#pragma once

/**
 * @file tsb.hpp
 * @brief Teunter-Syntetos-Babai (TSB) intermittent demand forecasting method.
 *
 * Reference:
 *   Teunter, R.H., Syntetos, A.A. & Babai, M.Z. (2011). Intermittent demand:
 *   Linking forecasting to inventory obsolescence. European Journal of
 *   Operational Research, 214(3), 606-615.
 */

#include <vector>

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/types/initialization_results.hpp"

namespace mro::forecasting::parametrics {

	/**
	 * @brief Output of the TSB forecasting method.
	 *
	 * Unlike Croston/SBA/SBJ, TSB does not implement a Croston wrapper:
	 * it smooths a per-period demand *probability* (updated every period,
	 * including zero-demand periods) instead of a demand *interval*
	 * (updated only on demand-occurring periods). This makes TSB
	 * responsive to obsolescence: probability decays toward zero during
	 * extended runs of zero demand, which a demand-interval-based method
	 * cannot represent.
	 */
	struct TsbResult : public mro::forecasting::types::ForecastResult {
		double z_hat;       ///< Smoothed estimate of demand size.
		double p_hat;       ///< Smoothed estimate of demand probability per period (0, 1].
		double forecast;    ///< Forecast per period = z_hat * p_hat.
		double sigma2_hat;  ///< Sample variance of the non-zero demand sizes observed.
		double alpha;       ///< Smoothing constant for demand size.
		double beta;        ///< Smoothing constant for demand probability.

		/// @copydoc ForecastResult::point_forecast
		double point_forecast() const override { return forecast; }

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Based on the variance of a Bernoulli(p_hat) occurrence combined
		 * with Normal(z_hat, sigma2_hat) demand sizes when demand occurs:
		 *
		 *   Var = p_hat * sigma2_hat + p_hat * (1 - p_hat) * z_hat^2
		 *
		 * This is the variance of the per-period compound process
		 * implied by the model's current state, analogous to
		 * CrostonResult::demand_variance() but built from a per-period
		 * Bernoulli probability instead of a geometric inter-arrival time.
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Forecasts intermittent demand using the Teunter-Syntetos-Babai
	 *        (TSB) method.
	 *
	 * Demand size (z_hat) is updated only on periods with non-zero demand,
	 * as in Croston. Demand probability (p_hat) is updated every period,
	 * including zero-demand periods:
	 *
	 *   z_hat  = alpha * y_t + (1 - alpha) * z_hat            (only if y_t > 0)
	 *   p_hat  = beta * 1{y_t > 0} + (1 - beta) * p_hat       (every period)
	 *   forecast = z_hat * p_hat
	 *
	 * @param time_series  Chronological sequence of period demands (zeros allowed).
	 * @param init         Initial values for z_hat and the demand
	 *                     probability (init.prob1), from any initialization
	 *                     method (e.g. calculate_by_moments()).
	 * @param alpha        Smoothing constant for demand size, 0 < alpha <= 1.
	 * @param beta         Smoothing constant for demand probability, 0 < beta <= 1.
	 * @param skip_periods Number of periods at the beginning of the series to skip
	 *                     (typically the periods used during initialization).
	 * @return TsbResult with the final z_hat, p_hat, and the resulting forecast.
	 * @throws std::invalid_argument if @p alpha is not in (0, 1], @p beta is not
	 *         in (0, 1], @p init.prob1 is not in (0, 1], or @p skip_periods is
	 *         greater than or equal to @p time_series.size().
	 */
	TsbResult forecast_by_tsb(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			double beta,
			size_t skip_periods = 0);

} // namespace mro::forecasting::parametrics
