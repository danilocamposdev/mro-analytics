#include "mro/initialization/seed_generators/moments.hpp"

#include <stdexcept>

namespace mro::initialization::seed_generators {

	MomentsParameters calculate_by_moments(const std::vector<double>& time_series,
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

		double sum     = 0.0;
		size_t nonzero = 0;

		for (size_t i = 0; i < initial_periods; ++i) {
			if (time_series[i] > 0.0) {
				sum += time_series[i];
				++nonzero;
			}
		}

		MomentsParameters result;
		result.z1 = (nonzero > 0) ? sum / nonzero : 0.0;
		result.prob1 = static_cast<double>(nonzero) / initial_periods;
		result.p1 = 1.0 / result.prob1;

		return result;
	}

} // namespace mro::initialization::seed_generators
