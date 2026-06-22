#pragma once

/**
 * @file gamma_demand_size.hpp
 * @brief Gamma-distributed demand size forecasting method (per period).
 *
 * Models the size of an individual non-zero demand as Gamma(shape, scale),
 * combined with the per-period probability of occurrence (analogous to
 * Croston's z_hat * (1/p_hat) decomposition, but with z_hat replaced by
 * the mean of a fitted Gamma distribution instead of a smoothed average).
 *
 * For modelling lead-time-aggregated demand directly with a Gamma
 * distribution instead (e.g. fitting Gamma to historical per-lead-time
 * totals), see gamma_lead_time.hpp.
 */

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/parameter_estimators/mle_gamma.hpp"

namespace mro::forecasting::distributions {

	/**
	 * @brief Output of the Gamma demand-size forecasting method.
	 */
	struct GammaDemandResult : public mro::forecasting::types::ForecastResult {
		double shape;  ///< Gamma shape parameter, k.
		double scale;  ///< Gamma scale parameter, theta.
		double prob1;  ///< Per-period probability of demand occurrence.

		/**
		 * @brief Point forecast for a single period ahead.
		 *
		 * E[X] = prob1 * (shape * scale)
		 */
		double point_forecast() const override;

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Variance of a Bernoulli(prob1)-gated Gamma(shape, scale):
		 *   Var = prob1 * shape * scale^2
		 *       + prob1 * (1 - prob1) * (shape * scale)^2
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Builds a GammaDemandResult from a fitted Gamma model and a
	 *        demand occurrence probability.
	 *
	 * This is a thin wrapper: shape/scale estimation is performed by
	 * mro::initialization::parameter_estimators::calculate_by_mle_gamma(),
	 * and prob1 typically comes from the same initialization window
	 * (e.g. mro::initialization::types::SeedGeneratorsResults::prob1).
	 *
	 * @param params Fitted (shape, scale), e.g. from calculate_by_mle_gamma().
	 * @param prob1  Per-period probability of demand occurrence, in (0, 1].
	 * @return GammaDemandResult with point_forecast() and demand_variance()
	 *         ready to use, including via mro::metrics::evaluate_policy().
	 * @throws std::invalid_argument if @p prob1 is not in (0, 1].
	 */
	GammaDemandResult forecast_by_gamma_demand_size(
			const mro::initialization::parameter_estimators::GammaParameters& params,
			double prob1);

} // namespace mro::forecasting::distributions
