#pragma once

/**
 * @file calibration_detail.hpp
 * @brief Shared backtesting/grid-search machinery used by every
 *        damping-constant calibration method (MASE, MAD, MSE).
 *
 * Not part of the public API: this header is included only by the
 * calibration .cpp files. It factors out the one part that MASE, MAD
 * and MSE calibration have in common — running a fixed-origin
 * backtest over a grid of candidate damping constants — so each
 * metric's .cpp only needs to supply how a single period's error is
 * scored.
 *
 * Backtesting scheme (fixed-origin, as opposed to one-step-ahead
 * rolling re-estimation): the initialization window is the first
 * @p skip_periods observations, exactly as in every forecast_by_*()
 * function. For each candidate damping constant(s), and for every test
 * period t from skip_periods to time_series.size()-1, the forecasting
 * method is re-run on the prefix time_series[0..t) with the same
 * skip_periods and candidate constant(s), and its forecast is compared
 * against the actual demand time_series[t]. This re-derives the
 * smoothed state from scratch for every t (the forecasting methods
 * don't expose incremental state), which is the price of treating them
 * as black boxes — see the note on computational cost in the public
 * calibration headers.
 */

#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include "mro/initialization/types/initialization_results.hpp"

namespace mro::calibration::detail {

	/// Signature shared by single-damping-constant forecasting methods
	/// (Croston, SBA, SBJ): given a series prefix and alpha, returns the
	/// point forecast for the next period.
	using SingleParamForecastFn = std::function<double(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			std::size_t skip_periods)>;

	/// Signature shared by dual-damping-constant forecasting methods (TSB).
	using DualParamForecastFn = std::function<double(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			double beta,
			std::size_t skip_periods)>;

	/// Validates the inputs common to every calibration entry point.
	inline void validate_common_inputs(
			const std::vector<double>& time_series,
			std::size_t skip_periods,
			double grid_step) {
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		if (skip_periods <= 1) {
			throw std::invalid_argument("skip_periods must be greater than 1.");
		}
		if (skip_periods >= time_series.size()) {
			throw std::invalid_argument(
					"skip_periods cannot be greater than or equal to time_series size.");
		}
		if (skip_periods + 1 >= time_series.size()) {
			throw std::invalid_argument(
					"calibration requires at least one test period after skip_periods "
					"(skip_periods must leave at least 2 trailing observations, since "
					"the first test forecast itself requires more than skip_periods "
					"observations in its prefix).");
		}
		if (grid_step <= 0.0 || grid_step >= 1.0) {
			throw std::invalid_argument("grid_step must be in (0, 1).");
		}
	}

	/// Mean absolute difference between consecutive observations in the
	/// initialization window, used as the MASE scaling denominator
	/// (the in-sample naive/lag-1 MAD).
	inline double naive_mad_in_sample(
			const std::vector<double>& time_series, std::size_t skip_periods) {
		double sum_abs_diff = 0.0;
		for (std::size_t i = 1; i < skip_periods; ++i) {
			sum_abs_diff += std::abs(time_series[i] - time_series[i - 1]);
		}
		return sum_abs_diff / static_cast<double>(skip_periods - 1);
	}

	/// Runs a fixed-origin backtest of a single-parameter forecasting
	/// method at a given alpha, returning the per-period forecast errors
	/// (actual - forecast) for every test period from skip_periods+1 to
	/// time_series.size()-1.
	///
	/// Test periods start at skip_periods+1, not skip_periods: the
	/// shortest valid prefix passed to forecast_fn must contain *more*
	/// than skip_periods observations (every forecast_by_*() function
	/// rejects skip_periods >= time_series.size()), so the first
	/// available test target is time_series[skip_periods+1], forecast
	/// from the prefix time_series[0..skip_periods+1).
	inline std::vector<double> backtest_errors_single_param(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			std::size_t skip_periods,
			const SingleParamForecastFn& forecast_fn) {
		std::vector<double> errors;
		errors.reserve(time_series.size() - skip_periods - 1);

		for (std::size_t t = skip_periods + 1; t < time_series.size(); ++t) {
			const std::vector<double> prefix(time_series.begin(), time_series.begin() + static_cast<long>(t));
			const double forecast = forecast_fn(prefix, init, alpha, skip_periods);
			errors.push_back(time_series[t] - forecast);
		}

		return errors;
	}

	/// Same as above, for a dual-parameter forecasting method (TSB).
	inline std::vector<double> backtest_errors_dual_param(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			double alpha,
			double beta,
			std::size_t skip_periods,
			const DualParamForecastFn& forecast_fn) {
		std::vector<double> errors;
		errors.reserve(time_series.size() - skip_periods - 1);

		for (std::size_t t = skip_periods + 1; t < time_series.size(); ++t) {
			const std::vector<double> prefix(time_series.begin(), time_series.begin() + static_cast<long>(t));
			const double forecast = forecast_fn(prefix, init, alpha, beta, skip_periods);
			errors.push_back(time_series[t] - forecast);
		}

		return errors;
	}

	/// Grid-searches alpha in (0, 1) (exclusive endpoints, step grid_step)
	/// for the value minimizing score_fn(errors), where errors are the
	/// fixed-origin backtest errors at that alpha.
	template <typename ScoreFn>
	inline std::pair<double, double> grid_search_single_param(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			std::size_t skip_periods,
			double grid_step,
			const SingleParamForecastFn& forecast_fn,
			ScoreFn score_fn) {
		double best_alpha = grid_step;
		double best_score = std::numeric_limits<double>::infinity();

		for (double alpha = grid_step; alpha < 1.0; alpha += grid_step) {
			const std::vector<double> errors =
					backtest_errors_single_param(time_series, init, alpha, skip_periods, forecast_fn);
			const double score = score_fn(errors);

			if (score < best_score) {
				best_score = score;
				best_alpha = alpha;
			}
		}

		return {best_alpha, best_score};
	}

	/// Grid-searches (alpha, beta) jointly in (0, 1) x (0, 1) for the pair
	/// minimizing score_fn(errors), where errors are the fixed-origin
	/// backtest errors at that (alpha, beta).
	template <typename ScoreFn>
	inline std::tuple<double, double, double> grid_search_dual_param(
			const std::vector<double>& time_series,
			const mro::initialization::types::SeedGeneratorsResults& init,
			std::size_t skip_periods,
			double grid_step,
			const DualParamForecastFn& forecast_fn,
			ScoreFn score_fn) {
		double best_alpha = grid_step;
		double best_beta = grid_step;
		double best_score = std::numeric_limits<double>::infinity();

		for (double alpha = grid_step; alpha < 1.0; alpha += grid_step) {
			for (double beta = grid_step; beta < 1.0; beta += grid_step) {
				const std::vector<double> errors = backtest_errors_dual_param(
						time_series, init, alpha, beta, skip_periods, forecast_fn);
				const double score = score_fn(errors);

				if (score < best_score) {
					best_score = score;
					best_alpha = alpha;
					best_beta = beta;
				}
			}
		}

		return {best_alpha, best_beta, best_score};
	}

} // namespace mro::calibration::detail
