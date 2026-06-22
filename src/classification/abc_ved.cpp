#include "mro/classification/abc_ved.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace mro::classification {

	namespace {

		char abc_letter(AbcCategory c) {
			switch (c) {
				case AbcCategory::A: return 'A';
				case AbcCategory::B: return 'B';
				case AbcCategory::C: return 'C';
			}
			return '?';
		}

		char ved_letter(VedCategory v) {
			switch (v) {
				case VedCategory::Vital:     return 'V';
				case VedCategory::Essential: return 'E';
				case VedCategory::Desirable: return 'D';
			}
			return '?';
		}

	} // namespace

	std::vector<AbcResult> classify_by_abc(
			const std::vector<double>& annual_values,
			const AbcThresholds& thresholds)
	{
		if (annual_values.empty()) {
			throw std::invalid_argument("annual_values cannot be empty.");
		}
		for (double v : annual_values) {
			if (v < 0.0) {
				throw std::invalid_argument("annual_values must be non-negative.");
			}
		}
		if (thresholds.a_cutoff <= 0.0 || thresholds.a_cutoff > 1.0 ||
				thresholds.b_cutoff <= 0.0 || thresholds.b_cutoff > 1.0) {
			throw std::invalid_argument("a_cutoff and b_cutoff must be in (0, 1].");
		}
		if (thresholds.a_cutoff >= thresholds.b_cutoff) {
			throw std::invalid_argument("a_cutoff must be less than b_cutoff.");
		}

		const size_t n = annual_values.size();
		const double total = std::accumulate(annual_values.begin(), annual_values.end(), 0.0);

		// Rank items by descending value; stable_sort keeps deterministic
		// tie-breaking (by original input order) when values are equal.
		std::vector<size_t> order(n);
		for (size_t i = 0; i < n; ++i) {
			order[i] = i;
		}
		std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
				return annual_values[a] > annual_values[b];
				});

		std::vector<AbcResult> result(n);
		double cumulative = 0.0;

		for (size_t rank = 0; rank < n; ++rank) {
			size_t idx = order[rank];
			double value = annual_values[idx];
			double share = (total > 0.0) ? (value / total) : 0.0;
			cumulative += share;

			AbcCategory cat;
			if (total <= 0.0) {
				// All values are zero: there is no meaningful Pareto curve;
				// nothing differentiates the items, so default to C.
				cat = AbcCategory::C;
			} else if (cumulative <= thresholds.a_cutoff) {
				cat = AbcCategory::A;
			} else if (cumulative <= thresholds.b_cutoff) {
				cat = AbcCategory::B;
			} else {
				cat = AbcCategory::C;
			}

			result[idx] = AbcResult{ cat, value, share, cumulative };
		}

		return result;
	}

	std::string combine_abc_ved(AbcCategory abc, VedCategory ved) {
		return std::string{ abc_letter(abc), ved_letter(ved) };
	}

} // namespace mro::classification
