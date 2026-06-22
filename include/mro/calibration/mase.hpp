#pragma once

/**
 * @file mase.hpp
 * @brief Calibrates damping constants by minimizing MASE (Mean Absolute
 *        Scaled Error) over a fixed-origin backtest.
 *
 * MASE expresses the forecast error relative to the error of a naive
 * lag-1 ("no-change") forecast measured in-sample, on the
 * initialization window:
 *
 *   MASE = mean(|actual_t - forecast_t|) / MAD_naive_in_sample
 *
 * where MAD_naive_in_sample = mean(|X_i - X_{i-1}|) over the first
 * skip_periods observations. This makes MASE scale-free and comparable
 * across items with different demand magnitudes, which MAD and MSE are
 * not (see Hyndman & Koehler, 2006, on why MASE was proposed for
 * exactly this reason).
 *
 * Backtesting scheme: fixed-origin. The initialization window is the
 * first skip_periods observations (the same skip_periods every
 * forecast_by_*() function already takes). For each candidate damping
 * constant on a grid over (0, 1), the forecasting method is re-run on
 * every prefix time_series[0..t) for t from skip_periods to
 * time_series.size()-1, and the resulting forecast is compared against
 * the actual time_series[t]. The grid value minimizing MASE is
 * returned. No re-estimation of the initialization window itself is
 * performed.
 *
 * Reference:
 *   Hyndman, R.J., Koehler, A.B. (2006). Another look at measures of
 *   forecast accuracy. International Journal of Forecasting, 22(4),
 *   679-688.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "mro/calibration/types/calibration_result.hpp"
#include "mro/initialization/types/initialization_results.hpp"

namespace mro::calibration {

	/**
	 * @brief Calibrates a single damping constant (e.g. for Croston, SBA
	 *        or SBJ) by minimizing MASE over a fixed-origin backtest.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param init Initialization parameters passed through unchanged to
	 * every backtest call to @p forecast_fn (e.g. from a moments/simple-average estimator).
	 * @param skip_periods Number of leading periods treated as the fixed
	 * initialization window (same role as in every forecast_by_*() call).
	 * @param forecast_fn The forecasting method to calibrate; must return
	 * a value whose point_forecast() reflects the input alpha (e.g.
	 * `[](const auto& ts, const auto& i, double a, std::size_t sp) { return
	 * forecast_by_croston(ts, i, a, sp).point_forecast(); }`).
	 * @param grid_step Grid resolution for the alpha search over (0, 1),
	 * exclusive of the endpoints (e.g. 0.01 searches 0.01, 0.02, ..., 0.99).
	 * @return SingleParamCalibResult with the alpha minimizing MASE and
	 * the MASE value attained at it.
	 * @throws std::invalid_argument if @p time_series is empty,
	 * @p skip_periods <= 1, @p skip_periods leaves fewer than 2 trailing
	 * test observations, or @p grid_step is not in (0, 1).
	 */
	mro::calibration::types::SingleParamCalibResult calibrate_by_mase(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			std::size_t skip_periods,
			const std::function<double(
					const std::vector<double>&,
					const mro::initialization::types::SeedGeneratorsResults&,
					double,
					std::size_t)>& forecast_fn,
			double grid_step = 0.01);

	/**
	 * @brief Calibrates two damping constants (e.g. for TSB) by minimizing
	 *        MASE over a fixed-origin backtest.
	 *
	 * Same backtesting scheme as calibrate_by_mase() above, but searching
	 * the (alpha, beta) grid jointly over (0, 1) x (0, 1).
	 *
	 * @param forecast_fn The forecasting method to calibrate; must return
	 * a value whose point_forecast() reflects the input (alpha, beta)
	 * (e.g. `[](const auto& ts, const auto& i, double a, double b, std::size_t sp) {
	 * return forecast_by_tsb(ts, i, a, b, sp).point_forecast(); }`).
	 * @param grid_step Grid resolution for both alpha and beta over
	 * (0, 1). Note the search is O((1/grid_step)^2) backtests, each
	 * itself O(time_series.size()); keep grid_step coarse for long series.
	 * @return DualParamCalibResult with the (alpha, beta) pair minimizing
	 * MASE and the MASE value attained at it.
	 * @throws std::invalid_argument under the same conditions as
	 * calibrate_by_mase() above.
	 */
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
			double grid_step = 0.01);

} // namespace mro::calibration
