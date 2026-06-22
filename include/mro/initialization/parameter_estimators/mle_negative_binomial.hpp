#pragma once

/**
 * @file mle_negative_binomial.hpp
 * @brief Maximum-likelihood estimation of a Negative Binomial NB(r, mu)
 *        model for intermittent/overdispersed demand.
 *
 * The Negative Binomial is parameterized here by the mean mu and the
 * dispersion parameter r (sometimes called the "size"), so that for a
 * count X ~ NB(r, mu):
 *
 *   P(X = k) = Gamma(k + r) / (Gamma(r) * k!) * (r / (r + mu))^r * (mu / (r + mu))^k
 *
 *   E[X]   = mu
 *   Var[X] = mu + mu^2 / r
 *
 * As r -> infinity the model approaches Poisson(mu); finite r captures
 * the overdispersion (Var > mean) typical of MRO demand, including the
 * extra zeros that a plain Poisson cannot reproduce.
 *
 * There is no closed-form MLE for (r, mu) jointly. However, for any
 * fixed r the mu that maximizes the log-likelihood is exactly the
 * sample mean (the same identity that makes mu the MLE in a Poisson
 * model). Substituting mu = mean(x) into the log-likelihood yields a
 * profile log-likelihood in r alone:
 *
 *   l(r) = sum_i [ lgamma(x_i + r) - lgamma(r) ]
 *          - n * (r * log(r + mu) - r * log(r))   [grouped over i]
 *          + sum_i [ x_i * log(mu) - x_i * log(r + mu) ] - sum_i lgamma(x_i + 1)
 *
 * which is solved for r via Newton-Raphson on l'(r) = 0, using the
 * digamma function (the same dependency already used by the Gamma
 * estimator), starting from a method-of-moments initial guess.
 *
 * Reference:
 *   Hilbe, J. M. (2011). Negative Binomial Regression (2nd ed.).
 *   Cambridge University Press.
 */

#include <cstddef>
#include <vector>

namespace mro::initialization::parameter_estimators {

	/**
	 * @brief Output of the Negative Binomial maximum-likelihood estimation.
	 */
	struct NegativeBinomialParameters {
		double r;                ///< Estimated dispersion parameter, r (r -> infinity recovers Poisson).
		double mu;                ///< Estimated mean, mu.
		double log_likelihood;    ///< Log-likelihood attained at (r, mu).
	};

	/**
	 * @brief Estimates (r, mu) of a Negative Binomial model from the first
	 *        @p initial_periods observations (zeros included, analogous
	 *        to calculate_by_mle_zip()).
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of leading periods used for the estimate.
	 * @param max_iterations Maximum number of Newton-Raphson iterations on r.
	 * @param tolerance Convergence tolerance on r between iterations.
	 * @return NegativeBinomialParameters with r, mu and the attained log-likelihood.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, @p initial_periods > time_series.size(),
	 *         every demand in the initialization window is zero (mu would
	 *         be zero and the model degenerates), or the sample is
	 *         under-dispersed relative to the Poisson (sample variance <=
	 *         sample mean), in which case the Negative Binomial MLE for r
	 *         diverges to infinity and a Poisson model should be used instead.
	 */
	NegativeBinomialParameters calculate_by_mle_negative_binomial(
			const std::vector<double>& time_series,
			std::size_t initial_periods,
			std::size_t max_iterations = 200,
			double tolerance = 1e-9);

} // namespace mro::initialization::parameter_estimators
