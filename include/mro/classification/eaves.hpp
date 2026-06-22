#pragma once

/**
 * @file eaves.hpp
 * @brief Eaves & Kingsman (2004) demand classification for spare parts.
 *
 * Reference:
 *   Eaves, A.H.C. & Kingsman, B.G. (2004). Forecasting for the ordering
 *   and stock-holding of spare parts. Journal of the Operational Research
 *   Society, 55(4), 431-437.
 *
 * Eaves re-estimated Williams' (1984) cut-off values empirically on
 * spare-parts data, proposing ADI = 1.25 and CV² = 0.49 as better
 * discriminators of which forecasting method (SES, Croston, etc.) performs
 * best on a given item. The classification logic is a simple 2x2 matrix
 * (ADI x CV²), without the extra CV² split used by Williams to produce a
 * "Highly Lumpy" category.
 *
 * This method defines its own category enum and result struct; it does
 * NOT reuse types from williams.hpp, so each classification method can
 * evolve independently.
 */

#include <string>
#include <vector>

namespace mro::classification {

/**
 * @brief Demand pattern categories as defined by Eaves & Kingsman (2004).
 *
 * @verbatim
 *               CV² <= 0.49       CV² > 0.49
 *  ADI <= 1.25   Regular           Erratic
 *  ADI >  1.25   SlowMovers        Lumpy
 * @endverbatim
 */
enum class EavesCategory {
    Regular    = 1, ///< Frequent and stable.
    Erratic    = 2, ///< Frequent but volatile in quantity.
    SlowMovers = 3, ///< Rare but stable in quantity.
    Lumpy      = 4  ///< Sparse and variable.
};

/** @brief Empirical cut-off values for the Eaves & Kingsman classification. */
namespace eaves_thresholds {
    constexpr double ADI_INTERMITTENCY = 1.25;
    constexpr double CV2_LUMPINESS     = 0.49;
}

/**
 * @brief Output of classify_by_eaves().
 *
 * Carries the classification decision and the intermediate statistics
 * that produced it.
 */
struct EavesClassificationResult {
    EavesCategory category;      ///< Demand pattern class (enum).
    std::string   category_name; ///< Human-readable label (e.g. "Lumpy").
    double        adi;           ///< Average Demand Interval = 1/p.
    double        cv2;           ///< CV² used for classification.
    double        p;             ///< Demand probability per period.
    double        mean_nonzero;  ///< Mean demand size on non-zero periods only.
};

/**
 * @brief Classifies a demand time series using Eaves & Kingsman (2004).
 *
 * ADI and CV² are computed on non-zero demand periods only. An item with
 * no non-zero demand is classified as SlowMovers.
 *
 * @param time_series Chronological sequence of period demands (zeros allowed).
 * @return EavesClassificationResult with category, label, ADI, CV², p, and mean_nonzero.
 * @throws std::invalid_argument if @p time_series is empty.
 */
EavesClassificationResult classify_by_eaves(const std::vector<double>& time_series);

} // namespace mro::classification
