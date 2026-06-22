#include "mro/preprocessing/negative_binomial_overdispersion.hpp"

#include "mro/initialization/parameter_estimators/mle_negative_binomial.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstddef>

namespace mro::preprocessing {

	namespace {

		/// Type-7 (linear-interpolation) quantile, matching R's and
		/// NumPy's default and mro::preprocessing::treat_outliers_adjusted_boxplot()'s
		/// internal convention. Assumes @p sorted_data is already sorted.
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

		void validate_band(const std::pair<double, double>& band, const char* band_name) {
			const auto [low, high] = band;
			if (low < 0.0 || low > 100.0 || high < 0.0 || high > 100.0) {
				throw std::invalid_argument(
						std::string("treat_outliers_negative_binomial: ") + band_name +
						" percentiles must be in [0, 100].");
			}
			if (low >= high) {
				throw std::invalid_argument(
						std::string("treat_outliers_negative_binomial: ") + band_name +
						" low percentile must be < high percentile.");
			}
		}

	} // namespace

	NegativeBinomialOverdispersionResult treat_outliers_negative_binomial(
			const std::vector<double>& time_series,
			double aggressive_dispersion_threshold,
			std::pair<double, double> normal_band,
			std::pair<double, double> aggressive_band,
			bool exclude_zeros)
	{
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		validate_band(normal_band, "normal_band");
		validate_band(aggressive_band, "aggressive_band");

		// Sample used for both the NB fit and the Winsorization
		// percentiles: by default, only the non-zero demands (zeros
		// reflect absence of demand, not an extreme observation, for
		// intermittent MRO series -- same convention as
		// treat_outliers_adjusted_boxplot()).
		std::vector<double> sample;
		for (double v : time_series) {
			if (!exclude_zeros || v != 0.0) {
				sample.push_back(v);
			}
		}

		// `sample` may be empty if every value was zero and excluded;
		// the Winsorization percentiles and the NB fit are both
		// meaningless in that case.
		if (sample.empty()) {
			throw std::invalid_argument(
					"treat_outliers_negative_binomial: every value in the series was "
					"zero (and excluded); there is nothing left to fit or treat.");
		}

		// mu and the sample variance are computed directly here (a
		// trivial O(n) pass) rather than always delegating to
		// calculate_by_mle_negative_binomial(): that function throws
		// std::invalid_argument when the sample is NOT overdispersed
		// (Var <= mu), since the NB shape MLE for r is undefined/diverges
		// there -- which is the expected, common case for roughly the
		// non-overdispersed tail of any real demand portfolio. Treating
		// that as a hard failure here would make outlier treatment
		// unusable for those items, so the overdispersion check is done
		// upfront, and the (comparatively expensive) Newton-Raphson MLE
		// is only invoked when it is actually well-defined.
		double sum_x = 0.0;
		for (double v : sample) {
			sum_x += v;
		}
		const double mu = sum_x / static_cast<double>(sample.size());

		double sum_sq_dev = 0.0;
		for (double v : sample) {
			sum_sq_dev += (v - mu) * (v - mu);
		}
		const double sample_variance = sum_sq_dev / static_cast<double>(sample.size());
		const double ratio = mu > 0.0 ? (sample_variance / mu) : std::numeric_limits<double>::quiet_NaN();
		const bool is_overdispersed = sample_variance > mu + 1e-9; // same tolerance as the MLE's own check

		// r is only meaningful (finite) when the sample is overdispersed;
		// calculate_by_mle_negative_binomial() is only called in that
		// case, both to avoid its exception and to avoid paying for the
		// Newton-Raphson solve when the answer is known to be +infinity.
		double r = std::numeric_limits<double>::infinity();
		if (is_overdispersed) {
			const mro::initialization::parameter_estimators::NegativeBinomialParameters params =
					mro::initialization::parameter_estimators::calculate_by_mle_negative_binomial(
							sample, sample.size());
			r = params.r;
		}

		const bool use_aggressive_band =
				is_overdispersed && ratio >= aggressive_dispersion_threshold;
		const std::pair<double, double>& chosen_band =
				use_aggressive_band ? aggressive_band : normal_band;

		std::vector<double> sorted_sample = sample;
		std::sort(sorted_sample.begin(), sorted_sample.end());
		const double lower_bound = std::max(
				quantile_of_sorted(sorted_sample, chosen_band.first / 100.0), 0.0);
		const double upper_bound = quantile_of_sorted(sorted_sample, chosen_band.second / 100.0);

		NegativeBinomialOverdispersionResult result;
		result.r = r;
		result.mu = mu;
		result.dispersion_ratio = ratio;
		result.is_overdispersed = is_overdispersed;
		result.used_aggressive_band = use_aggressive_band;
		result.treated_series.reserve(time_series.size());
		result.is_outlier.reserve(time_series.size());

		for (double v : time_series) {
			if (exclude_zeros && v == 0.0) {
				result.treated_series.push_back(v);
				result.is_outlier.push_back(false);
				continue;
			}

			if (v < lower_bound) {
				result.treated_series.push_back(lower_bound);
				result.is_outlier.push_back(true);
			} else if (v > upper_bound) {
				result.treated_series.push_back(upper_bound);
				result.is_outlier.push_back(true);
			} else {
				result.treated_series.push_back(v);
				result.is_outlier.push_back(false);
			}
		}

		return result;
	}

} // namespace mro::preprocessing
