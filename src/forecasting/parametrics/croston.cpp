#include "mro/forecasting/parametrics/croston.hpp"

#include <stdexcept>

namespace mro::forecasting::parametrics {

	double CrostonResult::demand_variance() const {
		// Var = alpha/(2-alpha) * [ (p-1)/p^2 * mu^2 + sigma^2/p ]
		// (variance of the compound Bernoulli-arrival / Normal-size demand
		// process implied by the model; see Croston (1972), eq. for Var(y)).
		return (alpha / (2.0 - alpha)) *
			(((p_hat - 1.0) / (p_hat * p_hat)) * z_hat * z_hat + sigma2_hat / p_hat);
	}

	CrostonResult forecast_by_croston(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			size_t skip_periods)
	{
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		if (skip_periods >= time_series.size()) {
			throw std::invalid_argument("skip_periods cannot be greater than or equal to time_series size.");
		}
		if (alpha <= 0.0 || alpha > 1.0) {
			throw std::invalid_argument("alpha must be in (0, 1].");
		}
		if (init.p1 <= 0.0) {
			throw std::invalid_argument("init.p1 must be positive.");
		}

		double z_hat = init.z1;
		double p_hat = init.p1;

		size_t periods_since_demand = 0;

		// Running stats over the non-zero demands actually processed
		// (post skip_periods), used to estimate sigma^2 for demand_variance().
		size_t n_positive = 0;
		double sum_positive = 0.0;
		double sum_sq_positive = 0.0;

		for (size_t i = skip_periods; i < time_series.size(); ++i) {
			double val = time_series[i];
			++periods_since_demand;

			if (val > 0.0) {
				z_hat = alpha * val + (1.0 - alpha) * z_hat;
				p_hat = alpha * static_cast<double>(periods_since_demand) +
					(1.0 - alpha) * p_hat;

				periods_since_demand = 0;

				++n_positive;
				sum_positive += val;
				sum_sq_positive += val * val;
			}
		}
		double forecast = z_hat / p_hat;

		CrostonResult result;
		result.z_hat = z_hat;
		result.p_hat = p_hat;
		result.forecast = forecast;
		result.alpha = alpha;

		if (n_positive >= 2) {
			const double mean_positive = sum_positive / static_cast<double>(n_positive);
			result.sigma2_hat =
				(sum_sq_positive - static_cast<double>(n_positive) * mean_positive * mean_positive) /
				static_cast<double>(n_positive - 1);
		} else {
			result.sigma2_hat = 0.0;
		}

		return result;
	}

} // namespace mro::forecasting::parametrics
