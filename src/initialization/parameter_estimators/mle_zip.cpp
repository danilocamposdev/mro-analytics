#include "mro/initialization/parameter_estimators/mle_zip.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mro::initialization::parameter_estimators {

	ZipParameters calculate_by_mle_zip(
			const std::vector<double>& time_series,
			std::size_t initial_periods,
			std::size_t max_iterations,
			double tolerance)
	{
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		if (initial_periods <= 1) {
			throw std::invalid_argument("initial_periods must be greater than 1.");
		}
		if (initial_periods > time_series.size()) {
			throw std::invalid_argument("initial_periods exceeds time series size.");
		}

		const std::vector<double> window(
				time_series.begin(), time_series.begin() + static_cast<long>(initial_periods));

		const double n = static_cast<double>(window.size());
		double sum_x = 0.0;
		std::size_t n_zero = 0;
		double sum_nonzero = 0.0;
		std::size_t n_nonzero = 0;

		for (double val : window) {
			sum_x += val;
			if (val == 0.0) {
				++n_zero;
			} else {
				sum_nonzero += val;
				++n_nonzero;
			}
		}

		if (sum_x <= 0.0) {
			throw std::invalid_argument(
					"calculate_by_mle_zip: every demand in the initialization window is "
					"zero; lambda is undefined.");
		}

		// Initial guess: lambda from the non-zero observations, pi from the
		// raw proportion of zeros (EM refines both from here).
		double lambda = sum_nonzero / static_cast<double>(n_nonzero);
		double pi = std::clamp(
				static_cast<double>(n_zero) / n, 1e-6, 1.0 - 1e-6);

		for (std::size_t iter = 0; iter < max_iterations; ++iter) {
			// E-step: every zero observation shares the same responsibility
			// under this i.i.d. model, so it is computed once.
			const double exp_neg_lambda = std::exp(-lambda);
			const double resp_zero = (n_zero > 0)
				? pi / (pi + (1.0 - pi) * exp_neg_lambda)
				: 0.0;
			const double sum_resp = resp_zero * static_cast<double>(n_zero);

			// M-step.
			const double pi_new = sum_resp / n;
			const double denom = n - sum_resp;
			const double lambda_new = (denom > 1e-12) ? (sum_x / denom) : lambda;

			const bool converged =
				std::abs(pi_new - pi) < tolerance && std::abs(lambda_new - lambda) < tolerance;

			pi = pi_new;
			lambda = lambda_new;

			if (converged) {
				break;
			}
		}

		double log_likelihood = 0.0;
		const double exp_neg_lambda = std::exp(-lambda);
		for (double val : window) {
			if (val == 0.0) {
				log_likelihood += std::log(pi + (1.0 - pi) * exp_neg_lambda);
			} else {
				log_likelihood += std::log(1.0 - pi) - lambda + val * std::log(lambda) -
					std::lgamma(val + 1.0);
			}
		}

		ZipParameters result;
		result.pi = pi;
		result.lambda = lambda;
		result.log_likelihood = log_likelihood;

		return result;
	}

} // namespace mro::initialization::parameter_estimators
