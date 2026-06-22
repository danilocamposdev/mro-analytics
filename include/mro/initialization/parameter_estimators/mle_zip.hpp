#pragma once

/**
 * @file zip.hpp
 * @brief Maximum-likelihood estimation of a Zero-Inflated Poisson (ZIP)
 *        model for intermittent demand.
 *
 * The ZIP model assumes each period's demand X is generated as:
 *   - with probability pi:       X = 0 (a "structural" zero)
 *   - with probability (1-pi):   X ~ Poisson(lambda)
 *
 * so that P(X = 0) = pi + (1-pi) e^{-lambda} and, for k >= 1,
 * P(X = k) = (1-pi) e^{-lambda} lambda^k / k!.
 *
 * There is no closed-form MLE for (pi, lambda) because the zero count
 * is a mixture of structural and Poisson-generated zeros, so the
 * estimate is obtained via Expectation-Maximization (EM):
 *
 *   E-step: for each zero observation, the probability that it is a
 *           structural zero (as opposed to a Poisson-generated zero):
 *             resp = pi / (pi + (1-pi) e^{-lambda})
 *   M-step: pi     <- mean(resp) over all periods (resp = 0 for non-zero periods)
 *           lambda <- sum(x) / (n - sum(resp))
 *
 * iterated to convergence.
 */

#include <cstddef>
#include <vector>

namespace mro::initialization::parameter_estimators {

	/**
	 * @brief Output of the ZIP maximum-likelihood estimation.
	 */
	struct ZipParameters {
		double pi;              ///< Estimated probability of a structural zero.
		double lambda;          ///< Estimated Poisson rate for non-structural periods.
		double log_likelihood;  ///< Log-likelihood attained at (pi, lambda).
	};

	/**
	 * @brief Estimates (pi, lambda) of a Zero-Inflated Poisson model via
	 *        Expectation-Maximization on the first @p initial_periods
	 *        observations.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of leading periods used for the estimate.
	 * @param max_iterations Maximum number of EM iterations.
	 * @param tolerance Convergence tolerance on both pi and lambda between iterations.
	 * @return ZipParameters with pi, lambda and the attained log-likelihood.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, @p initial_periods > time_series.size(),
	 *         or every demand in the initialization window is zero
	 *         (lambda would be undefined).
	 */
	ZipParameters calculate_by_mle_zip(
			const std::vector<double>& time_series,
			std::size_t initial_periods,
			std::size_t max_iterations = 200,
			double tolerance = 1e-9);

} // namespace mro::initialization::parameter_estimators
