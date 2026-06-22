#pragma once

/**
 * @file demand_stats.hpp
 * @brief Internal demand statistics shared across classification methods.
 * @note This header is not part of the public API.
 */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mro::classification::detail {

	struct DemandStats {
		size_t total_periods   = 0;
		size_t nonzero_periods = 0;
		double mean_nonzero    = 0.0;
		double var_nonzero     = 0.0;
		double p               = 0.0;
		double adi             = 0.0;
		double cv2             = 0.0;
	};

	inline DemandStats compute_stats(const std::vector<double>& series) {
		if (series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}

		DemandStats s;
		s.total_periods = series.size();

		// Pass 1: count non-zero periods and accumulate their sum
		double sum_nonzero = 0.0;
		for (double val : series) {
			if (val > 0.0) {
				++s.nonzero_periods;
				sum_nonzero += val;
			}
		}

		if (s.nonzero_periods == 0) {
			s.p   = 0.0;
			s.adi = std::numeric_limits<double>::infinity();
			return s;
		}

		s.mean_nonzero = sum_nonzero / static_cast<double>(s.nonzero_periods);

		// Pass 2: population variance over non-zero periods only
		double var_sum = 0.0;
		for (double val : series) {
			if (val > 0.0) {
				double diff = val - s.mean_nonzero;
				var_sum += diff * diff;
			}
		}
		s.var_nonzero = var_sum / static_cast<double>(s.nonzero_periods);

		double cv = (s.mean_nonzero > 0.0)
			? (std::sqrt(s.var_nonzero) / s.mean_nonzero)
			: 0.0;
		s.cv2 = cv * cv;

		s.p   = static_cast<double>(s.nonzero_periods) /
			static_cast<double>(s.total_periods);
		s.adi = 1.0 / s.p;

		return s;
	}

} // namespace mro::classification::detail
