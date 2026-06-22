#include "mro/calibration/mase.hpp"

#include <cmath>

#include "mro/calibration/detail/calibration_detail.hpp"

namespace mro::calibration {

	namespace {

		double mase_score(const std::vector<double>& errors, double naive_mad) {
			double sum_abs_error = 0.0;
			for (double err : errors) {
				sum_abs_error += std::abs(err);
			}
			const double mean_abs_error = sum_abs_error / static_cast<double>(errors.size());
			return mean_abs_error / naive_mad;
		}

	} // namespace

	mro::calibration::types::SingleParamCalibResult calibrate_by_mase(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			std::size_t skip_periods,
			const std::function<double(
					const std::vector<double>&,
					const mro::initialization::types::SeedGeneratorsResults&,
					double,
					std::size_t)>& forecast_fn,
			double grid_step) {
		mro::calibration::detail::validate_common_inputs(time_series, skip_periods, grid_step);

		const double naive_mad = mro::calibration::detail::naive_mad_in_sample(time_series, skip_periods);
		if (naive_mad <= 0.0) {
			throw std::invalid_argument(
					"calibrate_by_mase: the initialization window is constant (naive MAD "
					"is zero); MASE is undefined. Use calibrate_by_mad() or "
					"calibrate_by_mse() instead.");
		}

		const auto score_fn = [naive_mad](const std::vector<double>& errors) {
			return mase_score(errors, naive_mad);
		};

		const auto [best_alpha, best_score] = mro::calibration::detail::grid_search_single_param(
				time_series, init, skip_periods, grid_step, forecast_fn, score_fn);

		mro::calibration::types::SingleParamCalibResult result;
		result.alpha = best_alpha;
		result.minimum_error = best_score;
		return result;
	}

	mro::calibration::types::DualParamCalibResult calibrate_by_mase(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			std::size_t skip_periods,
			const std::function<double(
					const std::vector<double>&,
					const mro::initialization::types::SeedGeneratorsResults&,
					double,
					double,
					std::size_t)>& forecast_fn,
			double grid_step) {
		mro::calibration::detail::validate_common_inputs(time_series, skip_periods, grid_step);

		const double naive_mad = mro::calibration::detail::naive_mad_in_sample(time_series, skip_periods);
		if (naive_mad <= 0.0) {
			throw std::invalid_argument(
					"calibrate_by_mase: the initialization window is constant (naive MAD "
					"is zero); MASE is undefined. Use calibrate_by_mad() or "
					"calibrate_by_mse() instead.");
		}

		const auto score_fn = [naive_mad](const std::vector<double>& errors) {
			return mase_score(errors, naive_mad);
		};

		const auto [best_alpha, best_beta, best_score] = mro::calibration::detail::grid_search_dual_param(
				time_series, init, skip_periods, grid_step, forecast_fn, score_fn);

		mro::calibration::types::DualParamCalibResult result;
		result.alpha = best_alpha;
		result.beta = best_beta;
		result.minimum_error = best_score;
		return result;
	}

} // namespace mro::calibration
