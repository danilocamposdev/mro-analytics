#include "mro/initialization/parameter_estimators/mle_gamma.hpp"

#include <cmath>
#include <stdexcept>

namespace mro::initialization::parameter_estimators {

	namespace {

		/// Digamma function psi(x) = d/dx ln(Gamma(x)), via recurrence to
		/// push x >= 6 followed by the standard asymptotic series.
		double digamma(double x) {
			double result = 0.0;
			while (x < 6.0) {
				result -= 1.0 / x;
				x += 1.0;
			}
			const double inv_x2 = 1.0 / (x * x);
			result += std::log(x) - 0.5 / x -
					inv_x2 * (1.0 / 12.0 -
							inv_x2 * (1.0 / 120.0 -
									inv_x2 * (1.0 / 252.0 - inv_x2 * (1.0 / 240.0))));
			return result;
		}

		/// Trigamma function psi'(x) = d/dx psi(x), via the analogous
		/// recurrence/asymptotic-series approach.
		double trigamma(double x) {
			double result = 0.0;
			while (x < 6.0) {
				result += 1.0 / (x * x);
				x += 1.0;
			}
			const double inv_x = 1.0 / x;
			const double inv_x2 = inv_x * inv_x;
			result += inv_x + 0.5 * inv_x2 +
					inv_x2 * inv_x *
							(1.0 / 6.0 -
									inv_x2 * (1.0 / 30.0 -
											inv_x2 * (1.0 / 42.0 - inv_x2 * (1.0 / 30.0))));
			return result;
		}

	} // namespace

	GammaParameters calculate_by_mle_gamma(
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

		std::vector<double> positive_values;
		for (std::size_t i = 0; i < initial_periods; ++i) {
			if (time_series[i] > 0.0) {
				positive_values.push_back(time_series[i]);
			}
		}

		if (positive_values.size() < 2) {
			throw std::invalid_argument(
					"calculate_by_mle_gamma: at least 2 non-zero demands are required "
					"in the initialization window.");
		}

		const double n = static_cast<double>(positive_values.size());
		double sum_x = 0.0;
		double sum_log_x = 0.0;
		for (double val : positive_values) {
			sum_x += val;
			sum_log_x += std::log(val);
		}

		const double mean_x = sum_x / n;
		const double mean_log_x = sum_log_x / n;
		const double s = std::log(mean_x) - mean_log_x;

		if (s <= 1e-9) {
			throw std::invalid_argument(
					"calculate_by_mle_gamma: non-zero demand sizes have (near) zero "
					"variance; the gamma shape MLE is undefined (diverges to infinity).");
		}

		// Minka's closed-form initial approximation for the shape parameter.
		double shape = (3.0 - s + std::sqrt((s - 3.0) * (s - 3.0) + 24.0 * s)) / (12.0 * s);

		// Newton-Raphson refinement of log(k) - psi(k) = s.
		for (std::size_t iter = 0; iter < max_iterations; ++iter) {
			const double g = std::log(shape) - digamma(shape) - s;
			const double g_prime = 1.0 / shape - trigamma(shape);
			double next_shape = shape - g / g_prime;

			if (next_shape <= 0.0) {
				next_shape = shape / 2.0; // safeguard against Newton overshoot
			}

			const bool converged = std::abs(next_shape - shape) < tolerance;
			shape = next_shape;

			if (converged) {
				break;
			}
		}

		const double scale = mean_x / shape;
		const double log_likelihood = (shape - 1.0) * sum_log_x - sum_x / scale -
				n * shape * std::log(scale) - n * std::lgamma(shape);

		GammaParameters result;
		result.shape = shape;
		result.scale = scale;
		result.log_likelihood = log_likelihood;

		return result;
	}

} // namespace mro::initialization::parameter_estimators
