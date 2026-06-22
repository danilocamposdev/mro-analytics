#include "mro/forecasting/parametrics/tsb.hpp"

#include <stdexcept>

namespace mro::forecasting::parametrics {

	double TsbResult::demand_variance() const {
		// Var = p*sigma^2 + p*(1-p)*mu^2
		// (variance of a per-period Bernoulli(p) occurrence combined with
		// Normal(mu, sigma^2) demand sizes when demand occurs).
		return p_hat * sigma2_hat + p_hat * (1.0 - p_hat) * z_hat * z_hat;
	}

	TsbResult forecast_by_tsb(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			double beta,
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
		if (beta <= 0.0 || beta > 1.0) {
			throw std::invalid_argument("beta must be in (0, 1].");
		}
		if (init.prob1 <= 0.0 || init.prob1 > 1.0) {
			throw std::invalid_argument("init.prob1 must be in (0, 1].");
		}

		double z_hat = init.z1;
		double p_hat = init.prob1;

		// Running stats over the non-zero demands actually processed
		// (post skip_periods), used to estimate sigma^2 for demand_variance().
		size_t n_positive = 0;
		double sum_positive = 0.0;
		double sum_sq_positive = 0.0;

		for (size_t i = skip_periods; i < time_series.size(); ++i) {
			const double val = time_series[i];
			const bool occurred = val > 0.0;

			if (occurred) {
				z_hat = alpha * val + (1.0 - alpha) * z_hat;

				++n_positive;
				sum_positive += val;
				sum_sq_positive += val * val;
			}

			// p_hat is updated every period, including zero-demand ones.
			p_hat = beta * (occurred ? 1.0 : 0.0) + (1.0 - beta) * p_hat;
		}

		TsbResult result;
		result.z_hat = z_hat;
		result.p_hat = p_hat;
		result.forecast = z_hat * p_hat;
		result.alpha = alpha;
		result.beta = beta;

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
