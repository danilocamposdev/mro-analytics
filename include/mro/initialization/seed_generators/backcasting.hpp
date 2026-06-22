#pragma once

/**
 * @file backcasting.hpp
 * @brief Backcasting initialization method for intermittent-demand
 *        forecasting methods (Croston, SBA, TSB, SBJ, ...).
 *
 * Backcasting "warms up" the smoothed state before forward forecasting
 * begins, instead of seeding z1/p1 from a single, unweighted average
 * over the initialization window (as calculate_by_classic_ad_hoc()
 * does):
 *
 *   1. A naive seed is computed from the initialization window
 *      (mean of the non-zero demands for z, and window_length / count
 *      of non-zero demands for p — the same idea used by
 *      calculate_by_moments() for p1).
 *   2. The window is reversed (most recent period first) and the
 *      standard Croston smoothing recursion is run once over it,
 *      using this naive seed as the starting point.
 *   3. The resulting (z, p) after that backward pass become z1/p1.
 *
 * This lets the initial estimate reflect the dynamics of the whole
 * initialization window, with the earliest period in the window
 * (the one adjacent to where forward forecasting will start) carrying
 * the most weight in the final smoothed seed.
 */

#include <cstddef>
#include <vector>

#include "mro/initialization/types/initialization_results.hpp"

namespace mro::initialization::seed_generators {

	/**
	 * @brief Output of the backcasting initialization method.
	 *
	 * No additional fields beyond SeedGeneratorsResults are needed for
	 * this method; the struct exists to keep the same pattern used by
	 * the other initialization methods (ClassicAdHocParameters,
	 * MomentsParameters, ...).
	 */
	struct BackcastingParameters : mro::initialization::types::SeedGeneratorsResults {};

	/**
	 * @brief Computes z1 and p1 using the backcasting initialization method.
	 *
	 * @param time_series     Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of periods to use for initialization.
	 * @param alpha           Smoothing constant used for the backward warm-up pass,
	 *                        0 < alpha <= 1 (typically the same alpha used afterward
	 *                        for forward forecasting, so the seed is consistent with it).
	 * @return BackcastingParameters with z1 and p1.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, @p initial_periods > time_series.size(),
	 *         @p alpha is not in (0, 1], or there are no non-zero demands in the
	 *         initialization window (z1 undefined).
	 */
	BackcastingParameters calculate_by_backcasting(
			const std::vector<double>& time_series,
			size_t initial_periods,
			double alpha);

} // namespace mro::initialization::seed_generators
