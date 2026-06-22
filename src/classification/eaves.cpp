#include "mro/classification/eaves.hpp"
#include "detail/demand_stats.hpp"

namespace mro::classification {

	namespace {

		const char* label(EavesCategory cat) {
			switch (cat) {
				case EavesCategory::Regular:    return "Regular (Smooth)";
				case EavesCategory::Erratic:    return "Erratic";
				case EavesCategory::SlowMovers: return "Slow Movers";
				case EavesCategory::Lumpy:      return "Lumpy";
				default:                        return "Unknown";
			}
		}

		EavesClassificationResult make_result(EavesCategory cat,
				double adi,
				double cv2,
				double p,
				double mean_nonzero) {
			return EavesClassificationResult{ cat, label(cat), adi, cv2, p, mean_nonzero };
		}

		// Simple 2x2 matrix: Eaves & Kingsman do not split Lumpy further
		// into Lumpy/HighlyLumpy the way Williams does.
		EavesCategory evaluate_bounds(double adi, double cv2) {
			if (adi <= eaves_thresholds::ADI_INTERMITTENCY) {
				return (cv2 <= eaves_thresholds::CV2_LUMPINESS)
					? EavesCategory::Regular
					: EavesCategory::Erratic;
			}
			return (cv2 <= eaves_thresholds::CV2_LUMPINESS)
				? EavesCategory::SlowMovers
				: EavesCategory::Lumpy;
		}

	} // namespace

	EavesClassificationResult classify_by_eaves(const std::vector<double>& time_series) {
		detail::DemandStats s = detail::compute_stats(time_series);

		if (s.nonzero_periods == 0) {
			return make_result(
					EavesCategory::SlowMovers, s.adi, s.cv2, s.p, s.mean_nonzero);
		}

		EavesCategory cat = evaluate_bounds(s.adi, s.cv2);
		return make_result(cat, s.adi, s.cv2, s.p, s.mean_nonzero);
	}

} // namespace mro::classification
