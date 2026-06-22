#include "mro/calibration/mse.hpp"

#include "mro/calibration/detail/calibration_detail.hpp"

namespace mro::calibration {

	namespace {

		double mse_score(const std::vector<double>& errors) {
			double sum_sq_error = 0.0;
			for (double err : errors) {
				sum_sq_error += err * err;
			}
			return sum_sq_error / static_cast<double>(errors.size());
		}

	} // namespace

	mro::calibration::types::SingleParamCalibResult calibrate_by_mse(
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

		const auto [best_alpha, best_score] = mro::calibration::detail::grid_search_single_param(
				time_series, init, skip_periods, grid_step, forecast_fn, mse_score);

		mro::calibration::types::SingleParamCalibResult result;
		result.alpha = best_alpha;
		result.minimum_error = best_score;
		return result;
	}

	mro::calibration::types::DualParamCalibResult calibrate_by_mse(
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

		const auto [best_alpha, best_beta, best_score] = mro::calibration::detail::grid_search_dual_param(
				time_series, init, skip_periods, grid_step, forecast_fn, mse_score);

		mro::calibration::types::DualParamCalibResult result;
		result.alpha = best_alpha;
		result.beta = best_beta;
		result.minimum_error = best_score;
		return result;
	}

} // namespace mro::calibration
