#include "mro/preprocessing/adjusted_boxplot.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace mro::preprocessing {

	namespace {

		/// Type-7 (linear-interpolation) quantile, matching R's and
		/// NumPy's default. Assumes @p sorted_data is already sorted.
		double quantile_of_sorted(const std::vector<double>& sorted_data, double p) {
			const std::size_t n = sorted_data.size();
			if (n == 1) {
				return sorted_data[0];
			}

			const double pos = p * static_cast<double>(n - 1);
			const std::size_t lower = static_cast<std::size_t>(std::floor(pos));
			const std::size_t upper = static_cast<std::size_t>(std::ceil(pos));

			if (lower == upper) {
				return sorted_data[lower];
			}

			const double frac = pos - static_cast<double>(lower);
			return sorted_data[lower] + frac * (sorted_data[upper] - sorted_data[lower]);
		}

	} // namespace

	double medcouple(const std::vector<double>& data) {
		const std::size_t n = data.size();
		if (n < 3) {
			return 0.0;
		}

		std::vector<double> sorted_data = data;
		std::sort(sorted_data.begin(), sorted_data.end());

		const double med = quantile_of_sorted(sorted_data, 0.5);

		std::vector<double> z_plus;   // values >= median, ascending
		std::vector<double> z_minus;  // values <= median, ascending
		for (double v : sorted_data) {
			if (v >= med) z_plus.push_back(v);
			if (v <= med) z_minus.push_back(v);
		}

		const std::size_t n_plus = z_plus.size();
		const std::size_t n_minus = z_minus.size();

		// Size of the tied group exactly at the median (0 if the median
		// falls strictly between two distinct values).
		std::size_t k = 0;
		for (double v : sorted_data) {
			if (v == med) ++k;
		}

		std::vector<double> h_values;
		h_values.reserve(n_plus * n_minus);

		for (std::size_t i = 0; i < n_plus; ++i) {
			for (std::size_t j = 0; j < n_minus; ++j) {
				const double xi = z_plus[i];
				const double xj = z_minus[j];

				if (xi != xj) {
					h_values.push_back(((xi - med) - (med - xj)) / (xi - xj));
				} else if (k > 0 && i < k && j >= n_minus - k) {
					// Both values equal the median: standard tie-breaking
					// rule (Brys, Hubert & Struyf, 2004). i (0-indexed,
					// 0..k-1) is this value's rank within the tied group
					// as seen from z_plus; j_rank is its rank as seen
					// from z_minus (counting from the median outward).
					const double i_rank = static_cast<double>(i) + 1.0;          // 1..k
					const double j_rank = static_cast<double>(n_minus - j);      // 1..k
					const double rank_sum = i_rank + j_rank - 1.0;

					if (rank_sum < static_cast<double>(k)) {
						h_values.push_back(-1.0);
					} else if (rank_sum > static_cast<double>(k)) {
						h_values.push_back(1.0);
					} else {
						h_values.push_back(0.0);
					}
				}
			}
		}

		std::sort(h_values.begin(), h_values.end());
		return quantile_of_sorted(h_values, 0.5);
	}

	AdjustedBoxplotResult treat_outliers_adjusted_boxplot(
			const std::vector<double>& time_series,
			bool exclude_zeros,
			double whisker_factor)
	{
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		if (whisker_factor <= 0.0) {
			throw std::invalid_argument("whisker_factor must be > 0.");
		}

		// Sample used to estimate the fences: by default, only the
		// non-zero demands (zeros reflect absence of demand, not an
		// extreme observation, for intermittent MRO series).
		std::vector<double> sample;
		for (double v : time_series) {
			if (!exclude_zeros || v != 0.0) {
				sample.push_back(v);
			}
		}

		if (sample.size() < 3) {
			throw std::invalid_argument(
					"treat_outliers_adjusted_boxplot: at least 3 values (after "
					"excluding zeros, if requested) are required to compute the "
					"adjusted boxplot fences.");
		}

		std::sort(sample.begin(), sample.end());

		const double q1 = quantile_of_sorted(sample, 0.25);
		const double q3 = quantile_of_sorted(sample, 0.75);
		const double iqr = q3 - q1;
		const double mc = medcouple(sample);

		double lower_fence;
		double upper_fence;
		if (mc >= 0.0) {
			lower_fence = q1 - whisker_factor * std::exp(-4.0 * mc) * iqr;
			upper_fence = q3 + whisker_factor * std::exp(3.0 * mc) * iqr;
		} else {
			lower_fence = q1 - whisker_factor * std::exp(-3.0 * mc) * iqr;
			upper_fence = q3 + whisker_factor * std::exp(4.0 * mc) * iqr;
		}

		AdjustedBoxplotResult result;
		result.q1 = q1;
		result.q3 = q3;
		result.iqr = iqr;
		result.medcouple_value = mc;
		result.lower_fence = lower_fence;
		result.upper_fence = upper_fence;
		result.treated_series.reserve(time_series.size());
		result.is_outlier.reserve(time_series.size());

		for (double v : time_series) {
			if (exclude_zeros && v == 0.0) {
				result.treated_series.push_back(v);
				result.is_outlier.push_back(false);
				continue;
			}

			if (v < lower_fence) {
				result.treated_series.push_back(lower_fence);
				result.is_outlier.push_back(true);
			} else if (v > upper_fence) {
				result.treated_series.push_back(upper_fence);
				result.is_outlier.push_back(true);
			} else {
				result.treated_series.push_back(v);
				result.is_outlier.push_back(false);
			}
		}

		return result;
	}

} // namespace mro::preprocessing
