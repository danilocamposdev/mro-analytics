#pragma once

/**
 * @file moments.hpp
 * @brief Moments method for forecast initialization.
 */

#include <vector>

#include "mro/initialization/types/initialization_results.hpp"

namespace mro::initialization::seed_generators {

	/**
	 * @brief Output of the moments initialization method.
	 */
	struct MomentsParameters : mro::initialization::types::SeedGeneratorsResults {};

	/**
	 * @brief Computes z1 and p1 using the moments initialization method.
	 *
	 * @param time_series     Chronological sequence of period demands (zeros allowed).
	 * @param initial_periods Number of periods to use for initialization.
	 * @return MomentsParameters with z1 and p1.
	 * @throws std::invalid_argument if @p time_series is empty,
	 *         @p initial_periods <= 1, or @p initial_periods > time_series.size().
	 */
	MomentsParameters calculate_by_moments(const std::vector<double>& time_series,
			size_t initial_periods);

} // namespace mro::initialization::seed_generators
