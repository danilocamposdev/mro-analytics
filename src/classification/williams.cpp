#include "mro/classification/williams.hpp"
#include "detail/demand_stats.hpp"
#include <cmath>
#include <stdexcept>

namespace mro::classification {

	namespace {

		const char* label(DemandCategory cat) {
			switch (cat) {
				case DemandCategory::Regular:     return "Regular (Smooth)";
				case DemandCategory::Erratic:     return "Erratic";
				case DemandCategory::SlowMovers:  return "Slow Movers";
				case DemandCategory::Lumpy:       return "Lumpy";
				case DemandCategory::HighlyLumpy: return "Highly Lumpy";
				default:                          return "Unknown";
			}
		}

		ClassificationResult make_result(DemandCategory cat,
				double adi,
				double cv2,
				double p,
				double mean_nonzero) {
			return ClassificationResult{ cat, label(cat), adi, cv2, p, mean_nonzero };
		}

		DemandCategory evaluate_bounds(double adi, double cv2) {
			if (adi <= thresholds::ADI_INTERMITTENCY) {
				return (cv2 <= thresholds::CV2_LUMPINESS)
					? DemandCategory::Regular
					: DemandCategory::Erratic;
			}
			if (cv2 <= thresholds::CV2_LUMPINESS) {
				return DemandCategory::SlowMovers;
			}
			return (cv2 <= thresholds::CV2_HIGHLY_LUMPY)
				? DemandCategory::Lumpy
				: DemandCategory::HighlyLumpy;
		}

	} // namespace

	ClassificationResult classify_by_williams(const std::vector<double>& time_series) {
		detail::DemandStats s = detail::compute_stats(time_series);

		if (s.nonzero_periods == 0) {
			return make_result(
					DemandCategory::SlowMovers, s.adi, s.cv2, s.p, s.mean_nonzero);
		}

		DemandCategory cat = evaluate_bounds(s.adi, s.cv2);
		return make_result(cat, s.adi, s.cv2, s.p, s.mean_nonzero);
	}

	ClassificationResult classify_by_williams_lead_time(
			const std::vector<double>& time_series,
			double lt_mean,
			double lt_std_dev)
	{
		if (lt_mean <= 0.0) {
			throw std::invalid_argument("lt_mean must be positive.");
		}
		if (lt_std_dev < 0.0) {
			throw std::invalid_argument("lt_std_dev must be non-negative.");
		}

		detail::DemandStats s = detail::compute_stats(time_series);

		if (s.nonzero_periods == 0) {
			return make_result(
					DemandCategory::SlowMovers, s.adi, s.cv2, s.p, s.mean_nonzero);
		}

		const double n_bar = s.p;
		const double z_bar = s.mean_nonzero;
		const double var_z = s.var_nonzero;
		const double var_n = s.p * (1.0 - s.p);
		const double L_bar = lt_mean;
		const double var_L = lt_std_dev * lt_std_dev;

		const double mu_dlt = n_bar * z_bar * L_bar;

		if (mu_dlt <= 0.0) {
			return make_result(
					DemandCategory::SlowMovers, s.adi, 0.0, s.p, s.mean_nonzero);
		}

		const double term1 = z_bar * z_bar * L_bar * var_n;
		const double term2 = n_bar * L_bar * var_z;
		const double term3 = n_bar * n_bar * z_bar * z_bar * var_L;

		const double var_dlt = term1 + term2 + term3;
		const double cv2_dlt = var_dlt / (mu_dlt * mu_dlt);

		DemandCategory cat = evaluate_bounds(s.adi, cv2_dlt);
		return make_result(cat, s.adi, cv2_dlt, s.p, s.mean_nonzero);
	}

} // namespace mro::classification
