#pragma once

/**
 * @file negative_binomial_overdispersion.hpp
 * @brief Negative-Binomial-based overdispersion check and outlier
 *        treatment for MRO demand series.
 *
 * Companion to adjusted_boxplot.hpp in mro::preprocessing, following
 * the same shape: a single function takes a raw time series and
 * returns both the diagnostic statistics and an already-treated
 * series, ready to feed into the forecasting methods.
 *
 * Unlike the adjusted boxplot (which only looks at the shape of the
 * empirical distribution via the medcouple), this method first checks
 * *whether* an item's demand is overdispersed relative to a Poisson
 * process (Var > mean) by fitting a Negative Binomial NB(r, mu) via
 * maximum likelihood -- reusing
 * mro::initialization::parameter_estimators::calculate_by_mle_negative_binomial()
 * directly, so the estimation logic lives in exactly one place in the
 * library. The fitted dispersion is then used to choose how
 * aggressively to Winsorize the series: items that are MORE
 * overdispersed than a caller-supplied threshold get a tighter
 * percentile band (more capping), everyone else gets a more
 * conservative band.
 *
 * The dispersion threshold is intentionally a parameter rather than
 * computed internally: deciding what counts as "more overdispersed
 * than usual" requires comparing across many items (e.g. the median
 * dispersion ratio among all overdispersed items in a portfolio),
 * which only the caller -- iterating over the whole dataset -- is in a
 * position to compute. This function treats one item at a time, by
 * design, exactly like treat_outliers_adjusted_boxplot().
 */

#include <utility>
#include <vector>

namespace mro::preprocessing {

	/**
	 * @brief Output of the Negative-Binomial-based overdispersion check
	 *        and outlier treatment.
	 */
	struct NegativeBinomialOverdispersionResult {
		double r;                  ///< Fitted NB dispersion parameter; +infinity if not overdispersed (Var <= mu), NaN only if mu <= 0 (see is_overdispersed/dispersion_ratio).
		double mu;                 ///< Fitted NB mean (the sample mean of the series).
		double dispersion_ratio;   ///< Var/mu (the Poisson boundary is 1.0); NaN if mu <= 0.
		bool is_overdispersed;     ///< True if Var > mu (a finite, meaningful r exists).
		bool used_aggressive_band; ///< True if the tighter percentile band was applied (see aggressive_dispersion_threshold).

		std::vector<double> treated_series;  ///< Winsorized series (same length/order as the input).
		std::vector<bool> is_outlier;        ///< Per-period flag: true where treated_series differs from the input.
	};

	/**
	 * @brief Fits a Negative Binomial to @p time_series, then Winsorizes
	 *        it using a percentile band chosen by how overdispersed the
	 *        fit turned out to be.
	 *
	 * Mirrors treat_outliers_adjusted_boxplot(): by default, zeros are
	 * excluded both from the NB fit's sample mean/variance and from the
	 * Winsorization percentiles (an intermittent-demand zero reflects
	 * absence of demand, not a low observation of a continuous
	 * variable), and are passed through to treated_series unchanged.
	 *
	 * Band selection: if the item is overdispersed (Var > mu) AND its
	 * dispersion_ratio is >= @p aggressive_dispersion_threshold, the
	 * tighter @p aggressive_band is used; otherwise @p normal_band is
	 * used (this includes every non-overdispersed item). Pass
	 * @p aggressive_dispersion_threshold = +infinity to always use
	 * @p normal_band (i.e. disable the aggressive band entirely).
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param aggressive_dispersion_threshold Dispersion-ratio (Var/mu)
	 * cutoff at or above which @p aggressive_band is used instead of
	 * @p normal_band. Typically the median dispersion_ratio among the
	 * overdispersed items in the caller's full dataset (see file-level
	 * note above), but any positive threshold is accepted.
	 * @param normal_band (low, high) percentile pair in [0, 100] used
	 * for items below the threshold, low < high. Defaults to (1.0, 99.0).
	 * @param aggressive_band (low, high) percentile pair in [0, 100]
	 * used for items at/above the threshold, low < high. Defaults to
	 * (5.0, 95.0).
	 * @param exclude_zeros If true (default), zeros are excluded from
	 * both the NB fit and the Winsorization percentiles, and are never
	 * altered in treated_series.
	 * @return NegativeBinomialOverdispersionResult with the NB fit
	 * (r = +infinity when the sample is not overdispersed -- this is a
	 * normal, expected outcome, not an error: most demand portfolios
	 * have some Poisson-like items, and they still get treated using
	 * @p normal_band), the dispersion ratio, which band was used, and
	 * the treated series.
	 * @throws std::invalid_argument if @p time_series is empty; if every
	 * value is zero and @p exclude_zeros leaves an empty sample to fit/
	 * treat; or if either percentile band is invalid (low/high not in
	 * [0, 100], or low >= high). Unlike calling
	 * mro::initialization::parameter_estimators::calculate_by_mle_negative_binomial()
	 * directly, a non-overdispersed sample does NOT throw here -- see
	 * the r = +infinity note above.
	 */
	NegativeBinomialOverdispersionResult treat_outliers_negative_binomial(
			const std::vector<double>& time_series,
			double aggressive_dispersion_threshold,
			std::pair<double, double> normal_band = {1.0, 99.0},
			std::pair<double, double> aggressive_band = {5.0, 95.0},
			bool exclude_zeros = true);

} // namespace mro::preprocessing
