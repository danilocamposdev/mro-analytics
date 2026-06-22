#pragma once

/**
 * @file gamma.hpp
 * @brief Maximum-likelihood estimation of a Gamma(shape, scale) model
 *        for the size of non-zero demands.
 *
 * The Gamma distribution has support x > 0, so it is fit to the
 * non-zero demands found in the initialization window (analogous to
 * how calculate_by_classic_ad_hoc() computes z1 from the non-zero
 * demands only).
 *
 * There is no closed-form MLE for the shape parameter k. Writing
 * mean(x) and mean(log(x)) for the sample moments, the MLE satisfies:
 *
 *   log(k) - psi(k) = log(mean(x)) - mean(log(x))   =: s
 *
 * where psi is the digamma function. This is solved for k via
 * Newton-Raphson, starting from Minka's closed-form approximation,
 * then the scale follows from theta = mean(x) / k.
 *
 * Reference:
 *   Minka, T. (2002). Estimating a Gamma distribution. Microsoft Research.
 */

#include <cstddef>
#include <vector>

namespace mro::initialization::parameter_estimators {

	/**
	 * @brief Output of the Gamma maximum-likelihood estimation.
	 */
	struct GammaParameters {
		double shape;            ///< Estimated shape parameter, k.
		double scale;            ///< Estimated scale parameter, theta (mean = k * theta).
		double log_likelihood;   ///< Log-likelihood attained at (shape, scale).
	};

	/**
	 * @brief Estimates (shape, scale) of a Gamma model from the non-zero
	 *        demands in the first @p initial_periods observations.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of leading periods used for the estimate.
	 * @param max_iterations Maximum number of Newton-Raphson iterations.
	 * @param tolerance Convergence tolerance on the shape parameter between iterations.
	 * @return GammaParameters with shape, scale and the attained log-likelihood.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, @p initial_periods > time_series.size(),
	 *         fewer than 2 non-zero demands are found in the window, or
	 *         the non-zero demands have (near) zero variance (the shape
	 *         MLE diverges to infinity in that case).
	 */
	GammaParameters calculate_by_mle_gamma(
			const std::vector<double>& time_series,
			std::size_t initial_periods,
			std::size_t max_iterations = 100,
			double tolerance = 1e-10);

} // namespace mro::initialization::parameter_estimators
