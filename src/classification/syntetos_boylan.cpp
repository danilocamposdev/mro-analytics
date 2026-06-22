#include "mro/classification/syntetos_boylan.hpp"
#include "detail/demand_stats.hpp"

namespace mro::classification {

	namespace {

		const char* label(SbcCategory cat) {
			switch (cat) {
				case SbcCategory::Smooth:       return "Smooth";
				case SbcCategory::Erratic:      return "Erratic";
				case SbcCategory::Intermittent: return "Intermittent";
				case SbcCategory::Lumpy:        return "Lumpy";
				default:                        return "Unknown";
			}
		}

		SbcClassificationResult make_result(SbcCategory cat,
				double adi,
				double cv2,
				double p,
				double mean_nonzero) {
			return SbcClassificationResult{ cat, label(cat), adi, cv2, p, mean_nonzero };
		}

		SbcCategory evaluate_bounds(double adi, double cv2) {
			if (adi <= sbc_thresholds::ADI_INTERMITTENCY) {
				return (cv2 <= sbc_thresholds::CV2_LUMPINESS)
					? SbcCategory::Smooth
					: SbcCategory::Erratic;
			}
			return (cv2 <= sbc_thresholds::CV2_LUMPINESS)
				? SbcCategory::Intermittent
				: SbcCategory::Lumpy;
		}

	} // namespace

	SbcClassificationResult classify_by_syntetos_boylan(const std::vector<double>& time_series) {
		detail::DemandStats s = detail::compute_stats(time_series);

		if (s.nonzero_periods == 0) {
			return make_result(
					SbcCategory::Intermittent, s.adi, s.cv2, s.p, s.mean_nonzero);
		}

		SbcCategory cat = evaluate_bounds(s.adi, s.cv2);
		return make_result(cat, s.adi, s.cv2, s.p, s.mean_nonzero);
	}

} // namespace mro::classification
