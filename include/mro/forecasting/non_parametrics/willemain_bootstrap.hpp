#pragma once

/**
 * @file willemain_bootstrap.hpp
 * @brief Willemain-Smart-Schwarz (WSS) bootstrap method for intermittent
 *        demand, with jittering.
 *
 * Unlike the parametric methods (Croston/SBA/TSB/SBJ, ZIP, Gamma, ...),
 * this method makes no distributional assumption about demand. Instead,
 * it simulates the lead-time demand (LTD) directly, many times over, by
 * combining:
 *
 *   1. A two-state (zero / non-zero) Markov chain fit to the historical
 *      occurrence indicator, to generate an occurrence sequence I*(1..L)
 *      for the L periods of the lead-time window. This captures
 *      autocorrelation between consecutive occurrences that a method
 *      resampling periods independently (plain bootstrap) would miss.
 *   2. Bootstrap resampling (with replacement) of historical non-zero
 *      demand sizes for every period where I*(t) = 1.
 *   3. Jittering: replacing each resampled size Z_k by Z_k + eps, with
 *      eps ~ U(-h, h), to turn the discrete empirical distribution of
 *      sizes into a continuous one (avoids ties/steps in the simulated
 *      LTD distribution, which matters for percentiles/safety stock).
 *   4. Summing X*(t) = I*(t) * (Z_k + eps) over the L periods to obtain
 *      one simulated lead-time-demand draw, D_LT.
 *
 * Steps 1-4 are repeated n_simulations times (Willemain et al. suggest
 * on the order of 10,000) to build an empirical distribution of D_LT,
 * from which the mean, standard deviation, and any percentile of
 * interest can be read off directly — no Normality assumption is made
 * anywhere in the pipeline.
 *
 * Because the LTD is simulated directly over the L-period horizon
 * (rather than derived by scaling a per-period i.i.d. forecast),
 * WillemainBootstrapResult overrides ForecastResult::lead_time_demand()
 * to return the empirical mean/std-dev of the simulated D_LT sample
 * instead of the default sqrt(L)-scaling behavior.
 *
 * Three rules are offered for the jitter half-width h, all from the
 * accompanying course material:
 *   - kAlphaSigma (default, "original" jitter): h = alpha * sigma_Z,
 *     alpha in [0.1, 0.3] (Willemain's rule of thumb).
 *   - kIqr: h = c * IQR(Z), more robust to extreme outliers/heavy tails.
 *   - kSilverman: h = 0.9 * min(sigma_Z, IQR(Z) / 1.34) * n^(-1/5),
 *     Silverman's robust bandwidth rule.
 *
 * Reference:
 *   Willemain, T.R., Smart, C.N., Schwarz, H.F. (2004). A new approach
 *   to forecasting intermittent demand for service parts inventories.
 *   International Journal of Forecasting, 20(3), 375-387.
 */

#include <cstddef>
#include <vector>

#include "mro/forecasting/types/forecast_result.hpp"

namespace mro::forecasting::non_parametrics {

	/**
	 * @brief Rule used to set the jitter half-width h, where the
	 *        simulated size is Z_k + eps, eps ~ U(-h, h).
	 */
	enum class JitterStrategy {
		kAlphaSigma,  ///< h = jitter_alpha * sigma_Z (original WSS rule; default).
		kIqr,         ///< h = jitter_alpha * IQR(Z) (robust to outliers/heavy tails).
		kSilverman    ///< h = 0.9 * min(sigma_Z, IQR(Z) / 1.34) * n^(-1/5).
	};

	/**
	 * @brief Output of the Willemain bootstrap forecasting method.
	 *
	 * Wraps the empirical distribution of simulated lead-time demand
	 * (D_LT) so it can be consumed by mro::metrics like any other
	 * ForecastResult, while still exposing the full simulated sample
	 * for percentile-based uses (e.g. empirical P95, empirical fill
	 * rate) that the Normal-based metrics in mro::metrics don't cover.
	 */
	struct WillemainBootstrapResult : public mro::forecasting::types::ForecastResult {
		std::vector<double> simulated_lead_time_demand;  ///< The n_simulations draws of D_LT.
		double lead_time_periods;                        ///< L used to generate the simulation.

		/**
		 * @brief Point forecast for a single period ahead.
		 *
		 * Recovered from the simulated lead-time demand as
		 * mean(simulated_lead_time_demand) / lead_time_periods, i.e. the
		 * average per-period rate implied by the simulated LTD sample
		 * (consistent with mro::forecasting::types::ForecastResult's
		 * per-period contract; lead_time_demand() below is what should
		 * be used directly for LTD-based metrics).
		 */
		double point_forecast() const override;

		/**
		 * @brief Estimated variance of demand for a single period ahead.
		 *
		 * Recovered from the simulated lead-time demand as
		 * var(simulated_lead_time_demand) / lead_time_periods, assuming
		 * (only for this per-period decomposition) that lead-time periods
		 * contribute independently and identically to the total variance.
		 * As with point_forecast(), prefer lead_time_demand() directly
		 * for LTD-based metrics, since it uses the simulated sample
		 * as-is with no such assumption.
		 */
		double demand_variance() const override;

		/**
		 * @brief Mean and standard deviation of the simulated lead-time
		 *        demand sample.
		 *
		 * Overrides the i.i.d.-scaling default in ForecastResult: the
		 * lead-time demand here was simulated directly over the
		 * @p lead_time_periods horizon used to build this result (via
		 * the Markov-chain-driven occurrence sequence, so consecutive
		 * periods are NOT treated as independent), so it is returned
		 * as-is rather than recomputed by scaling point_forecast()/
		 * demand_variance() by the periods argument.
		 *
		 * @param lead_time_periods Ignored if it matches the
		 * lead_time_periods this result was built with; if it differs,
		 * the result was simulated for a different horizon and the
		 * sample does not apply, so the same construction must be
		 * re-run for the new horizon (the i.i.d. scaling fallback would
		 * silently misrepresent the Markov-driven autocorrelation
		 * captured in the simulation).
		 * @throws std::invalid_argument if @p lead_time_periods does not
		 * match (within floating-point tolerance) the lead_time_periods
		 * this result was built with.
		 */
		mro::forecasting::types::LeadTimeDemand lead_time_demand(
				double lead_time_periods) const override;

		/**
		 * @brief Empirical quantile of the simulated lead-time demand.
		 *
		 * Provided because, per the underlying method, percentiles, ROP,
		 * CSL and fill rate are meant to be read directly off the
		 * empirical distribution of simulated D_LT rather than assumed
		 * Normal (unlike mro::metrics, which assumes Normal lead-time
		 * demand throughout). Sorts a copy of
		 * simulated_lead_time_demand internally.
		 *
		 * @param probability Quantile level, in [0, 1] (e.g. 0.95 for P95).
		 * @throws std::invalid_argument if @p probability is not in [0, 1],
		 * or if simulated_lead_time_demand is empty.
		 */
		double empirical_quantile(double probability) const;
	};

	/**
	 * @brief Runs the WSS bootstrap-with-jitter simulation and builds a
	 *        WillemainBootstrapResult from it.
	 *
	 * Self-contained: the two-state Markov transition matrix (P00, P01,
	 * P10, P11) is estimated internally from the occurrence indicator of
	 * @p time_series, and the bootstrap pool of demand sizes is the set
	 * of non-zero values in @p time_series. No separate initialization
	 * step is required (unlike the Gamma/ZIP/Negative-Binomial methods),
	 * since resampling and jittering are themselves the estimation.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param lead_time_periods Lead time, in the same period unit as the
	 * demand series (e.g. quarters), > 0. Need not be an integer: the
	 * occurrence simulation runs for floor(lead_time_periods) full
	 * periods plus, if the fractional remainder is non-zero, one final
	 * partial period whose simulated size (if it occurs) is scaled by
	 * that remainder.
	 * @param n_simulations Number of Monte Carlo repetitions used to
	 * build the empirical D_LT distribution (Willemain et al. suggest
	 * on the order of 10,000).
	 * @param jitter_strategy Rule used to compute the jitter half-width h
	 * (see JitterStrategy). Defaults to kAlphaSigma, the original WSS rule.
	 * @param jitter_alpha Multiplier used by kAlphaSigma (h = jitter_alpha
	 * * sigma_Z) and kIqr (h = jitter_alpha * IQR(Z)). Ignored by
	 * kSilverman. Willemain's rule of thumb is jitter_alpha in [0.1, 0.3]
	 * (defaults to 0.2); the course material additionally suggests
	 * narrowing to [0.05, 0.10] for small/discrete-item samples and
	 * [0.15, 0.30] for continuous items when using kAlphaSigma.
	 * @param random_seed Seed for the random number generator, for
	 * reproducibility. If 0, a non-deterministic seed is used.
	 * @return WillemainBootstrapResult with point_forecast(),
	 * demand_variance(), lead_time_demand() and empirical_quantile()
	 * ready to use, including via mro::metrics::evaluate_policy()
	 * (though, per the note on lead_time_demand(), prefer reading
	 * percentiles directly from empirical_quantile() over the
	 * Normal-based mro::metrics pipeline when possible).
	 * @throws std::invalid_argument if @p time_series is empty,
	 * @p time_series has fewer than 2 observations (a transition cannot
	 * be estimated from a single observation), @p lead_time_periods <= 0,
	 * @p n_simulations == 0, @p jitter_alpha < 0, or fewer than 2 non-zero
	 * demands are found in @p time_series (the bootstrap pool of sizes
	 * would be degenerate).
	 */
	WillemainBootstrapResult forecast_by_willemain_bootstrap(
			const std::vector<double>& time_series,
			double lead_time_periods,
			std::size_t n_simulations = 10000,
			JitterStrategy jitter_strategy = JitterStrategy::kAlphaSigma,
			double jitter_alpha = 0.2,
			unsigned int random_seed = 0);

} // namespace mro::forecasting::non_parametrics
