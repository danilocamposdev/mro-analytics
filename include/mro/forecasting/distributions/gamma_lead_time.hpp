#pragma once

/**
 * @file gamma_lead_time.hpp
 * @brief Gamma-distributed lead-time demand (LTD) forecasting method.
 *
 * Models the total demand aggregated over a lead-time window directly
 * as Gamma(shape, scale), rather than building it up from a per-period
 * process. This is appropriate when the Gamma was fit to historical
 * observations already aggregated over windows of the relevant lead
 * time length (e.g. summing demand over rolling L-period windows
 * before fitting), as is common practice for lead-time demand modelling
 * in MRO inventory control (see e.g. Silver, Pyke & Thomas, 2017, on
 * fitting the LTD distribution directly).
 *
 * Because the model already describes the lead-time-aggregated demand,
 * it implements LeadTimeDistributionResult rather than ForecastResult:
 * there is no meaningful single-period decomposition to scale.
 */

#include "mro/forecasting/types/forecast_result.hpp"
#include "mro/initialization/parameter_estimators/mle_gamma.hpp"

namespace mro::forecasting::distributions {

	/**
	 * @brief Output of the Gamma lead-time-demand forecasting method.
	 */
	struct GammaLeadTimeResult : public mro::forecasting::types::LeadTimeDistributionResult {
		double shape;  ///< Gamma shape parameter, k, fit to lead-time-aggregated demand.
		double scale;  ///< Gamma scale parameter, theta.

		/**
		 * @brief Mean and standard deviation of the lead-time demand.
		 *
		 * mu_L    = shape * scale
		 * sigma_L = sqrt(shape) * scale
		 */
		mro::forecasting::types::LeadTimeDemand lead_time_demand() const override;
	};

	/**
	 * @brief Builds a GammaLeadTimeResult from a Gamma model fit to
	 *        lead-time-aggregated demand.
	 *
	 * This is a thin wrapper: estimation is performed by
	 * mro::initialization::parameter_estimators::calculate_by_mle_gamma(),
	 * called on a series of historical lead-time-window totals rather
	 * than the raw per-period series.
	 *
	 * @param params Fitted (shape, scale), e.g. from calculate_by_mle_gamma()
	 *               applied to lead-time-aggregated demand observations.
	 * @return GammaLeadTimeResult with lead_time_demand() ready to use,
	 *         including via mro::metrics::evaluate_policy().
	 */
	GammaLeadTimeResult forecast_by_gamma_lead_time(
			const mro::initialization::parameter_estimators::GammaParameters& params);

} // namespace mro::forecasting::distributions
