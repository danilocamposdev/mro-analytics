#pragma once

/**
 * @file adjusted_boxplot.hpp
 * @brief Adjusted boxplot outlier treatment for skewed distributions
 *        (Hubert & Vandervieren, 2008), applied to MRO demand series.
 *
 * Tukey's classic boxplot flags x as an outlier when:
 *   x < Q1 - 1.5*IQR   or   x > Q3 + 1.5*IQR
 *
 * This works well for roughly symmetric data, but intermittent/MRO
 * demand sizes are typically right-skewed, which makes the classic
 * rule flag too many points on the long (upper) tail as outliers.
 * Hubert & Vandervieren (2008) adjust the fences using the medcouple
 * (MC), a robust measure of skewness in [-1, 1]:
 *
 *   if MC >= 0: [Q1 - 1.5*e^(-4*MC)*IQR,  Q3 + 1.5*e^(3*MC)*IQR]
 *   if MC <  0: [Q1 - 1.5*e^(-3*MC)*IQR,  Q3 + 1.5*e^(4*MC)*IQR]
 *
 * which widens the fence on the long tail and narrows it on the short
 * one, recovering Tukey's original rule when MC = 0 (symmetric data).
 *
 * References:
 *   Hubert, M., Vandervieren, E. (2008). An adjusted boxplot for
 *   skewed distributions. Computational Statistics & Data Analysis,
 *   52(12), 5186-5201.
 *   Brys, G., Hubert, M., Struyf, A. (2004). A robust measure of
 *   skewness. Journal of Computational and Graphical Statistics, 13(4).
 *
 * Note: this method is purely empirical (order statistics + medcouple,
 * no fitted parametric distribution), so it is implemented with the
 * standard library only — it doesn't draw on StatsLib. StatsLib would
 * be the natural fit for a *parametric* outlier rule instead (e.g.
 * flagging based on quantiles of a fitted Gamma/ZIP model), which is
 * a different, complementary approach to this one.
 */

#include <vector>

namespace mro::preprocessing {

	/**
	 * @brief Output of the adjusted-boxplot outlier treatment.
	 */
	struct AdjustedBoxplotResult {
		double q1;               ///< First quartile of the sample used for the fences.
		double q3;                ///< Third quartile of the sample used for the fences.
		double iqr;               ///< Interquartile range, q3 - q1.
		double medcouple_value;   ///< Medcouple (robust skewness), in [-1, 1].
		double lower_fence;       ///< Skew-adjusted lower fence.
		double upper_fence;       ///< Skew-adjusted upper fence.

		/// Same length and period order as the input series. Values
		/// outside the fences are winsorized (clipped to the nearest
		/// fence); everything else (including any periods excluded
		/// from the fence calculation, e.g. zeros) is left unchanged.
		/// Ready to be fed directly into a forecasting method.
		std::vector<double> treated_series;

		/// Per-period flag, true where the original value was outside
		/// the fences (and was therefore winsorized in treated_series).
		std::vector<bool> is_outlier;
	};

	/**
	 * @brief Robust measure of skewness (medcouple) of a sample.
	 *
	 * MC = median over all pairs (xi, xj), xi <= median <= xj, xi != xj,
	 * of h(xi, xj) = ((xj - median) - (median - xi)) / (xj - xi), with
	 * the standard tie-breaking rule (Brys, Hubert & Struyf, 2004) for
	 * pairs where xi = xj = median.
	 *
	 * MC = 0 for symmetric data, MC > 0 for right-skewed data, MC < 0
	 * for left-skewed data, MC in [-1, 1].
	 *
	 * @param data Sample (need not be sorted).
	 * @return The medcouple; 0.0 if data.size() < 3 (skewness undefined).
	 */
	double medcouple(const std::vector<double>& data);

	/**
	 * @brief Detects and winsorizes outliers in a demand series using
	 *        the adjusted (skew-corrected) boxplot rule.
	 *
	 * @param time_series Chronological sequence of period demands.
	 * @param exclude_zeros If true (default), periods with zero demand
	 * are excluded from the fence computation and left untouched in
	 * the output — for intermittent MRO demand, zeros represent the
	 * absence of a transaction, not an extreme observation. If false,
	 * the full series (including zeros) is used.
	 * @param whisker_factor Whisker multiplier (Tukey's classic value
	 * is 1.5); larger values flag fewer points as outliers.
	 * @return AdjustedBoxplotResult with the fences and the winsorized series.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p whisker_factor <= 0, or fewer than 3 values remain
	 *         after excluding zeros (if requested) to compute the fences.
	 */
	AdjustedBoxplotResult treat_outliers_adjusted_boxplot(
			const std::vector<double>& time_series,
			bool exclude_zeros = true,
			double whisker_factor = 1.5);

} // namespace mro::preprocessing
