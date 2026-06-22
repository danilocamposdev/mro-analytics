#include "mro/initialization/parameter_estimators/mle_negative_binomial.hpp"

#include <cmath>
#include <stdexcept>

namespace mro::initialization::parameter_estimators {

	namespace {

		/// Digamma function psi(x) = d/dx ln(Gamma(x)), via recurrence to
		/// push x >= 6 followed by the standard asymptotic series.
		/// (Mirrors the implementation in mle_gamma.cpp; kept local here
		/// since there is currently no shared math-utilities header.)
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

	NegativeBinomialParameters calculate_by_mle_negative_binomial(
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
		for (double val : window) {
			sum_x += val;
		}

		if (sum_x <= 0.0) {
			throw std::invalid_argument(
					"calculate_by_mle_negative_binomial: every demand in the "
					"initialization window is zero; mu is undefined.");
		}

		// mu is fixed at its profile-MLE value (the sample mean) for every
		// candidate r, exactly as in the Poisson case.
		const double mu = sum_x / n;

		double sum_sq_dev = 0.0;
		for (double val : window) {
			sum_sq_dev += (val - mu) * (val - mu);
		}
		const double sample_variance = sum_sq_dev / n;

		if (sample_variance <= mu + 1e-9) {
			throw std::invalid_argument(
					"calculate_by_mle_negative_binomial: the initialization window "
					"is not overdispersed relative to the Poisson (sample variance "
					"<= sample mean); the negative binomial shape MLE for r diverges "
					"to infinity. Consider a Poisson model instead.");
		}

		// Method-of-moments initial guess from Var = mu + mu^2 / r,
		// i.e. r = mu^2 / (Var - mu).
		double r = (mu * mu) / (sample_variance - mu);

		// Newton-Raphson refinement of the profile log-likelihood's
		// first-order condition in r:
		//   l'(r) = sum_i [ psi(x_i + r) - psi(r) ] + n * log(r / (r + mu))
		//         + n * mu / (r + mu)
		// Differentiating again:
		//   l''(r) = sum_i [ psi'(x_i + r) - psi'(r) ] + n * mu / (r * (r + mu))
		//          - n * mu / (r + mu)^2
		for (std::size_t iter = 0; iter < max_iterations; ++iter) {
			double sum_digamma_diff = 0.0;
			double sum_trigamma_diff = 0.0;
			for (double val : window) {
				sum_digamma_diff += digamma(val + r) - digamma(r);
				sum_trigamma_diff += trigamma(val + r) - trigamma(r);
			}

			const double g = sum_digamma_diff + n * std::log(r / (r + mu));
			const double g_prime = sum_trigamma_diff + n * mu / (r * (r + mu) * (r + mu));

			double next_r = r - g / g_prime;

			if (next_r <= 0.0) {
				next_r = r / 2.0; // safeguard against Newton overshoot
			}

			const bool converged = std::abs(next_r - r) < tolerance;
			r = next_r;

			if (converged) {
				break;
			}
		}

		double log_likelihood = 0.0;
		for (double val : window) {
			log_likelihood += std::lgamma(val + r) - std::lgamma(r) - std::lgamma(val + 1.0) +
					r * std::log(r / (r + mu)) + val * std::log(mu / (r + mu));
		}

		NegativeBinomialParameters result;
		result.r = r;
		result.mu = mu;
		result.log_likelihood = log_likelihood;

		return result;
	}

} // namespace mro::initialization::parameter_estimators
