#include "mro/initialization/seed_generators/backcasting.hpp"

#include <algorithm>
#include <stdexcept>

namespace mro::initialization::seed_generators {

	BackcastingParameters calculate_by_backcasting(
			const std::vector<double>& time_series,
			size_t initial_periods,
			double alpha)
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
		if (alpha <= 0.0 || alpha > 1.0) {
			throw std::invalid_argument("alpha must be in (0, 1].");
		}

		// Initialization window, reversed so the warm-up pass runs from
		// the most recent period in the window back to the earliest one.
		std::vector<double> window(
				time_series.begin(), time_series.begin() + static_cast<long>(initial_periods));
		std::reverse(window.begin(), window.end());

		// Naive seed for the warm-up pass itself: mean of the non-zero
		// demands for z, and window_length / nonzero_count for p (same
		// idea as calculate_by_moments()'s p1 = 1 / prob1).
		double sum_demand = 0.0;
		size_t nonzero    = 0;
		for (double val : window) {
			if (val > 0.0) {
				sum_demand += val;
				++nonzero;
			}
		}

		if (nonzero == 0) {
			throw std::invalid_argument(
					"No non-zero demand found in the initialization window; cannot backcast.");
		}

		double z = sum_demand / static_cast<double>(nonzero);
		double p = static_cast<double>(initial_periods) / static_cast<double>(nonzero);

		// Single backward warm-up pass: the standard Croston smoothing
		// recursion, applied to the reversed window.
		size_t periods_since_demand = 0;
		for (double val : window) {
			++periods_since_demand;

			if (val > 0.0) {
				z = alpha * val + (1.0 - alpha) * z;
				p = alpha * static_cast<double>(periods_since_demand) + (1.0 - alpha) * p;

				periods_since_demand = 0;
			}
		}

		BackcastingParameters result;
		result.z1 = z;
		result.p1 = p;
		result.prob1 = 1.0 / p;

		return result;
	}

} // namespace mro::initialization::seed_generators
