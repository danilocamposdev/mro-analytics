#pragma once

/**
 * @file negative_binomial.hpp
 * @brief Negative Binomial intermittent demand forecasting method (per period).
 *
 * The Negative Binomial NB(r, mu) models each period's demand as an
 * overdispersed count: same mean mu as a Poisson(mu), but with extra
 * variance controlled by the dispersion parameter r (Var = mu + mu^2/r,
 * which collapses to the Poisson Var = mu as r -> infinity). This is a
 * common alternative to ZIP for MRO demand that shows overdispersion
 * without necessarily having an excess of structural zeros.
 *
 * Reference:
 *   Hilbe, J. M. (2011). Negative Binomial Regression (2nd ed.).
 *   Cambridge University Press.
 */

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/parameter_estimators/mle_negative_binomial.hpp"

namespace mro::forecasting::distributions {

	/**
	 * @brief Output of the Negative Binomial forecasting method.
	 *
	 * Wraps a fitted NB(r, mu) model as a per-period demand process, so
	 * it can be consumed by mro::metrics like any other ForecastResult.
	 */
	struct NegativeBinomialResult : public mro::forecasting::types::ForecastResult {
		double r;   ///< Dispersion parameter (r -> infinity recovers Poisson).
		double mu;  ///< Mean demand per period.

		/**
		 * @brief Point forecast for a single period ahead.
		 *
		 * E[X] = mu
		 */
		double point_forecast() const override;

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Var[X] = mu + mu^2 / r
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Builds a NegativeBinomialResult from a fitted NB model.
	 *
	 * This is a thin wrapper: parameter estimation itself is performed by
	 * mro::initialization::parameter_estimators::calculate_by_mle_negative_binomial(),
	 * which is passed in already fitted.
	 *
	 * @param params Fitted (r, mu), e.g. from calculate_by_mle_negative_binomial().
	 * @return NegativeBinomialResult with point_forecast() and demand_variance()
	 *         ready to use, including via mro::metrics::evaluate_policy().
	 */
	NegativeBinomialResult forecast_by_negative_binomial(
			const mro::initialization::parameter_estimators::NegativeBinomialParameters& params);

} // namespace mro::forecasting::distributions
