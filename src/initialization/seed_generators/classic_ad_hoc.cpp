#include "mro/initialization/seed_generators/classic_ad_hoc.hpp"

#include <stdexcept>

namespace mro::initialization::seed_generators {

	ClassicAdHocParameters calculate_by_classic_ad_hoc(const std::vector<double>& time_series,
			size_t initial_periods)
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

		double sum_demand   = 0.0;
		size_t nonzero      = 0;
		double sum_gaps     = 0.0;
		size_t num_gaps     = 0;
		bool   has_previous = false;
		size_t previous_idx = 0;

		for (size_t i = 0; i < initial_periods; ++i) {
			if (time_series[i] > 0.0) {
				sum_demand += time_series[i];
				++nonzero;

				if (has_previous) {
					sum_gaps += static_cast<double>(i - previous_idx);
					++num_gaps;
				}
				previous_idx = i;
				has_previous = true;
			}
		}

		if (nonzero == 0) {
			throw std::invalid_argument(
					"No non-zero demand found in the initialization window; cannot compute z1.");
		}
		if (num_gaps == 0) {
			throw std::invalid_argument(
					"At least two non-zero demands are required in the initialization "
					"window to compute p1.");
		}

		ClassicAdHocParameters result;
		result.z1 = sum_demand / static_cast<double>(nonzero);
		result.p1 = sum_gaps / static_cast<double>(num_gaps);
		result.prob1 = 1.0 / result.p1;

		return result;
	}

} // namespace mro::initialization::seed_generators
