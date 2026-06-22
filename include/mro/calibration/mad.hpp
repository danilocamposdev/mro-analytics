#pragma once

/**
 * @file mad.hpp
 * @brief Calibrates damping constants by minimizing MAD (Mean Absolute
 *        Deviation) over a fixed-origin backtest.
 *
 * MAD is the mean absolute forecast error over the test periods:
 *
 *   MAD = mean(|actual_t - forecast_t|)
 *
 * Unlike MASE, MAD is in the same units as the demand series itself, so
 * it is not comparable across items with different demand magnitudes;
 * it is the right choice when calibrating a single item (or items
 * already on a comparable scale) and not aggregating/ranking across a
 * heterogeneous portfolio (use calibrate_by_mase() for that instead).
 *
 * Backtesting scheme: identical to mase_calibration.hpp — fixed-origin,
 * with skip_periods as the initialization window and the grid search
 * run over (0, 1).
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "mro/calibration/types/calibration_result.hpp"
#include "mro/initialization/types/initialization_results.hpp"

namespace mro::calibration {

	/**
	 * @brief Calibrates a single damping constant (e.g. for Croston, SBA
	 *        or SBJ) by minimizing MAD over a fixed-origin backtest.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @param init Initialization parameters passed through unchanged to
	 * every backtest call to @p forecast_fn.
	 * @param skip_periods Number of leading periods treated as the fixed
	 * initialization window.
	 * @param forecast_fn The forecasting method to calibrate; must return
	 * the point forecast for the input alpha (see mase_calibration.hpp
	 * for an example lambda).
	 * @param grid_step Grid resolution for the alpha search over (0, 1),
	 * exclusive of the endpoints.
	 * @return SingleParamCalibResult with the alpha minimizing MAD and
	 * the MAD value attained at it.
	 * @throws std::invalid_argument if @p time_series is empty,
	 * @p skip_periods <= 1, @p skip_periods leaves fewer than 2 trailing
	 * test observations, or @p grid_step is not in (0, 1).
	 */
	mro::calibration::types::SingleParamCalibResult calibrate_by_mad(
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
	 *        MAD over a fixed-origin backtest.
	 *
	 * @param forecast_fn The forecasting method to calibrate; must return
	 * the point forecast for the input (alpha, beta).
	 * @param grid_step Grid resolution for both alpha and beta over
	 * (0, 1). Note the search is O((1/grid_step)^2) backtests.
	 * @return DualParamCalibResult with the (alpha, beta) pair minimizing
	 * MAD and the MAD value attained at it.
	 * @throws std::invalid_argument under the same conditions as
	 * calibrate_by_mad() above.
	 */
	mro::calibration::types::DualParamCalibResult calibrate_by_mad(
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
