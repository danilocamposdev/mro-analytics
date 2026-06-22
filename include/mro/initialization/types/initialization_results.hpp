#pragma once

/**
 * @file initialization_results.hpp
 * @brief Base type for forecast initialization parameters.
 */

namespace mro::initialization::types {

/**
 * @brief Base struct for forecast initialization parameters.
 *
 * Every initialization method (moments, simple average, etc.) returns a
 * struct derived from InitializationParameters. Forecasting methods (e.g.
 * Croston) accept this base type, so the initialization method used is
 * decoupled from the forecasting method applied.
 */
struct SeedGeneratorsResults {
    double z1; ///< Initial demand size estimate.
    double p1; ///< Initial demand interval.
	double prob1; ///< Initial demand probability.

    virtual ~SeedGeneratorsResults() = default;
};

} // namespace mro::initialization::types
