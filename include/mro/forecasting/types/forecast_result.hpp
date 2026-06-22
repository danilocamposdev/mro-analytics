#pragma once

/**
 * @file forecast_result.hpp
 * @brief Common interfaces implemented by every forecasting method's result.
 *
 * Decouples mro::metrics from any specific forecasting method: a
 * performance metric only needs a point forecast and a per-period
 * demand variance estimate (or, for methods where the i.i.d.-periods
 * assumption doesn't hold, a custom lead_time_demand() override) — it
 * doesn't care whether those came from Croston, SBA, TSB, SBJ, a
 * statistical-distribution model, or a bootstrap simulation.
 *
 * Some models (e.g. a Gamma fit to lead-time demand directly) describe
 * the lead-time-aggregated demand itself rather than a per-period
 * process, so a per-period point_forecast()/demand_variance() is not
 * meaningful for them. LeadTimeDistributionResult covers that case.
 */

#include <cmath>

namespace mro::forecasting::types {

	/**
	 * @brief Mean and standard deviation of demand aggregated over a
	 *        lead time window.
	 */
	struct LeadTimeDemand {
		double mean;     ///< mu_L: expected demand during the lead time.
		double std_dev;  ///< sigma_L: standard deviation of demand during the lead time.
	};

	/**
	 * @brief Polymorphic base for the result of any per-period
	 *        demand-forecasting method.
	 */
	struct ForecastResult {
		virtual ~ForecastResult() = default;

		/// Point forecast for a single period ahead.
		virtual double point_forecast() const = 0;

		/// Estimated variance of demand for a single period ahead.
		virtual double demand_variance() const = 0;

		/**
		 * @brief Mean and standard deviation of demand over a lead time
		 *        spanning @p lead_time_periods periods.
		 *
		 * Default implementation assumes i.i.d. periods, so both mean
		 * and variance scale linearly with the number of periods:
		 *   mu_L    = point_forecast() * L
		 *   sigma_L = sqrt(demand_variance() * L)
		 *
		 * Override this for methods where that assumption doesn't hold
		 * (e.g. simulation/bootstrap-based forecasts, where the
		 * lead-time distribution is built empirically instead).
		 *
		 * @param lead_time_periods Lead time, in the same period unit
		 * as the demand series (e.g. quarters), > 0.
		 */
		virtual LeadTimeDemand lead_time_demand(double lead_time_periods) const {
			return LeadTimeDemand{
				point_forecast() * lead_time_periods,
					std::sqrt(demand_variance() * lead_time_periods)};
		}
	};

	/**
	 * @brief Polymorphic base for models that describe the lead-time-
	 *        aggregated demand directly, rather than a per-period process.
	 *
	 * Some statistical models (e.g. a Gamma distribution fit to historical
	 * demand observed over lead-time-length windows) produce mu_L and
	 * sigma_L directly, with no meaningful single-period decomposition.
	 * For these, lead_time_demand() does not scale a per-period estimate;
	 * it is the estimate.
	 *
	 * mro::metrics::evaluate_policy() accepts either this or
	 * ForecastResult, so such models still plug into the same pipeline.
	 */
	struct LeadTimeDistributionResult {
		virtual ~LeadTimeDistributionResult() = default;

		/**
		 * @brief Mean and standard deviation of demand over the lead time
		 *        the model was fit to.
		 *
		 * Unlike ForecastResult::lead_time_demand(), this does not take a
		 * lead_time_periods argument: the lead time is intrinsic to how
		 * the underlying distribution was fit (e.g. a Gamma fit to
		 * historical per-lead-time-window aggregated demand), not a
		 * multiplier applied afterward.
		 */
		virtual LeadTimeDemand lead_time_demand() const = 0;
	};

} // namespace mro::forecasting::types
