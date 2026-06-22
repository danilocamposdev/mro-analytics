#pragma once

/**
 * @file williams.hpp
 * @brief Williams (1984) demand classification for MRO items.
 *
 * References:
 *   - Williams, T.M. (1984). Stock control with sporadic and slow-moving demand.
 *     Journal of the Operational Research Society, 35(10), 939-948.
 *   - Eaves, A.H.C. & Kingsman, B.G. (2004). Forecasting for the ordering and
 *     stock-holding of spare parts. Journal of the Operational Research Society.
 */

#include <string>
#include <vector>

namespace mro::classification {

/**
 * @brief Demand pattern categories as defined by Williams (1984).
 *
 * @verbatim
 *              CV² <= 0.50       CV² > 0.50
 *  ADI <= 2.80  Regular          Erratic
 *  ADI >  2.80  SlowMovers       Lumpy / HighlyLumpy (CV² > 0.70)
 * @endverbatim
 */
enum class DemandCategory {
    Regular     = 1, ///< Frequent and stable.
    Erratic     = 2, ///< Frequent but volatile in quantity.
    SlowMovers  = 3, ///< Rare but stable in quantity.
    Lumpy       = 4, ///< Sparse and variable.
    HighlyLumpy = 5  ///< Virtually unforecastable.
};

/** @brief Empirical cut-off values for the Williams classification. */
namespace thresholds {
    constexpr double ADI_INTERMITTENCY = 2.80;
    constexpr double CV2_LUMPINESS     = 0.50;
    constexpr double CV2_HIGHLY_LUMPY  = 0.70;
}

/**
 * @brief Output of classify_by_williams functions.
 *
 * Carries the classification decision and the intermediate statistics
 * that produced it.
 */
struct ClassificationResult {
    DemandCategory category;      ///< Demand pattern class (enum).
    std::string    category_name; ///< Human-readable label (e.g. "Lumpy").
    double         adi;           ///< Average Demand Interval = 1/p.
    double         cv2;           ///< CV² used for classification.
    double         p;             ///< Demand probability per period.
    double         mean_nonzero;  ///< Mean demand size on non-zero periods only.
};

/**
 * @brief Classifies a demand time series using Williams (1984).
 *
 * CV² and ADI are computed on non-zero demand periods only.
 * Lead-time variability is not considered; see classify_by_williams_lead_time().
 *
 * @param time_series Chronological sequence of period demands (zeros allowed).
 * @return ClassificationResult with category, label, ADI, CV², p, and mean_nonzero.
 * @throws std::invalid_argument if @p time_series is empty.
 */
ClassificationResult classify_by_williams(const std::vector<double>& time_series);

/**
 * @brief Classifies a demand time series using the full Williams (1984)
 *        variance decomposition with lead-time variability.
 *
 * var_DLT = z² * L * var_n  +  n * L * var_z  +  n² * z² * var_L
 *
 * CV²_DLT is used for classification in place of the plain demand CV².
 *
 * @param time_series Chronological sequence of period demands (zeros allowed).
 * @param lt_mean     Mean lead time, in the same time unit as the series.
 * @param lt_std_dev  Standard deviation of lead time.
 * @return ClassificationResult with category, label, ADI, CV²_DLT, p, and mean_nonzero.
 * @throws std::invalid_argument if @p time_series is empty,
 *         @p lt_mean <= 0, or @p lt_std_dev < 0.
 */
ClassificationResult classify_by_williams_lead_time(
    const std::vector<double>& time_series,
    double lt_mean,
    double lt_std_dev);

} // namespace mro::classification
