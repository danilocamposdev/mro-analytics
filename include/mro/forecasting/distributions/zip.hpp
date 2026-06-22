#pragma once

/**
 * @file zip.hpp
 * @brief Zero-Inflated Poisson (ZIP) intermittent demand forecasting method.
 *
 * The ZIP model treats each period's demand as a mixture: with
 * probability pi the period is a structural zero, and with probability
 * (1-pi) the demand is drawn from a Poisson(lambda). This captures
 * "more zeros than a plain Poisson would predict", which is typical of
 * intermittent MRO demand.
 *
 * Reference:
 *   Lambert, D. (1992). Zero-inflated Poisson regression, with an
 *   application to defects in manufacturing. Technometrics, 34(1), 1-14.
 */

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/parameter_estimators/mle_zip.hpp"

namespace mro::forecasting::distributions {

	/**
	 * @brief Output of the ZIP forecasting method.
	 *
	 * Wraps a fitted ZIP(pi, lambda) model as a per-period demand
	 * process, so it can be consumed by mro::metrics like any other
	 * ForecastResult.
	 */
	struct ZipResult : public mro::forecasting::types::ForecastResult {
		double pi;      ///< Probability of a structural zero.
		double lambda;  ///< Poisson rate for non-structural periods.

		/**
		 * @brief Point forecast for a single period ahead.
		 *
		 * E[X] = (1 - pi) * lambda
		 */
		double point_forecast() const override;

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Var[X] = (1 - pi) * lambda * (1 + pi * lambda)
		 */
		double demand_variance() const override;
	};

	/**
	 * @brief Builds a ZipResult from a fitted ZIP model.
	 *
	 * This is a thin wrapper: parameter estimation itself is performed by
	 * mro::initialization::parameter_estimators::calculate_by_mle_zip(),
	 * which is passed in already fitted.
	 *
	 * @param params Fitted (pi, lambda), e.g. from calculate_by_mle_zip().
	 * @return ZipResult with point_forecast() and demand_variance() ready
	 *         to use, including via mro::metrics::evaluate_policy().
	 */
	ZipResult forecast_by_zip(
			const mro::initialization::parameter_estimators::ZipParameters& params);

} // namespace mro::forecasting::distributions
