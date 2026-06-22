#include "mro/forecasting/parametrics/sba.hpp"
#include "mro/forecasting/parametrics/croston.hpp"

namespace mro::forecasting::parametrics {

	namespace {
		constexpr double bias_factor(double alpha) {
			return 1.0 - alpha / 2.0;
		}
	} // namespace

	double SbaResult::demand_variance() const {
		// Var(c * X) = c^2 * Var(X), where c = (1 - alpha/2) and Var(X) is
		// Croston's underlying demand variance (see CrostonResult::demand_variance()).
		const double c = bias_factor(alpha);
		const double croston_variance =
			(alpha / (2.0 - alpha)) *
			(((p_hat - 1.0) / (p_hat * p_hat)) * z_hat * z_hat + sigma2_hat / p_hat);
		return c * c * croston_variance;
	}

	SbaResult forecast_by_sba(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			size_t skip_periods)
	{
		const CrostonResult croston =
			forecast_by_croston(time_series, init, alpha, skip_periods);

		SbaResult result;
		result.z_hat = croston.z_hat;
		result.p_hat = croston.p_hat;
		result.sigma2_hat = croston.sigma2_hat;
		result.alpha = alpha;
		result.forecast = bias_factor(alpha) * croston.forecast;

		return result;
	}

} // namespace mro::forecasting::parametrics
