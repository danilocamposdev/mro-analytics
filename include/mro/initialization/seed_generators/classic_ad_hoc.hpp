#pragma once

/**
 * @file classic_ad_hoc.hpp
 * @brief Classic ad-hoc initialization method for Croston's (1972) method.
 *
 * This is the initialization scheme originally described by Croston (1972)
 * and commonly used in textbook/teaching examples: the initial demand size
 * (z1) is the arithmetic mean of the non-zero demands observed in the
 * initialization window, and the initial inter-demand interval (p1) is the
 * arithmetic mean of the gaps (in periods) between successive non-zero
 * demands within that same window.
 *
 * Note this differs from calculate_by_moments() (see moments.hpp), which
 * estimates p1 indirectly from the proportion of non-zero periods
 * (p1 = 1 / (nonzero / initial_periods)). Both are valid initialization
 * strategies found in the literature; they generally produce different
 * values of p1 for the same window.
 */

#include <vector>

#include "mro/initialization/types/initialization_results.hpp"

namespace mro::initialization::seed_generators {

	/**
	 * @brief Output of the classic ad-hoc initialization method.
	 *
	 * No additional fields beyond InitializationParameters are needed for
	 * this method; the struct exists to keep the same pattern used by other
	 * initialization methods (e.g. MomentsParameters) and to allow future
	 * extension without breaking callers.
	 */
	struct ClassicAdHocParameters : mro::initialization::types::SeedGeneratorsResults {};

	/**
	 * @brief Computes z1 and p1 using the classic ad-hoc initialization
	 *        method described by Croston (1972).
	 *
	 * z1 is the mean of the non-zero demands within the first
	 * @p initial_periods periods.
	 *
	 * p1 is the mean of the gaps (in periods) between successive non-zero
	 * demands within the same window. At least two non-zero demands are
	 * required within the window to compute at least one gap.
	 *
	 * @param time_series     Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of periods to use for initialization.
	 * @return ClassicAdHocParameters with z1 and p1.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, @p initial_periods > time_series.size(),
	 *         there are no non-zero demands in the window (z1 undefined), or
	 *         there is only one non-zero demand in the window (no gap to
	 *         compute p1 from).
	 */
	ClassicAdHocParameters calculate_by_classic_ad_hoc(const std::vector<double>& time_series,
			size_t initial_periods);

} // namespace mro::initialization::seed_generators
