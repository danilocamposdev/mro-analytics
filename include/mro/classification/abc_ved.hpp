#pragma once

/**
 * @file abc_ved.hpp
 * @brief ABC (Pareto value) and VED (criticality) inventory classification.
 *
 * Unlike Williams/Eaves/Syntetos-Boylan, this is NOT a statistical
 * demand-pattern classification derived from a single item's time series.
 *
 *  - ABC classifies a *portfolio* of items by their share of cumulative
 *    annual consumption value (typically unit_price * annual_demand),
 *    following the Pareto principle ("a small fraction of SKUs account
 *    for most of the value").
 *
 *  - VED classifies a single item by its operational criticality
 *    (Vital / Essential / Desirable). This is a qualitative judgment
 *    normally supplied by engineering/maintenance/operations -- it is
 *    NOT computed from demand data, so this module does not attempt to
 *    derive it; callers must provide it directly.
 *
 * The two are commonly combined into a 3x3 prioritization matrix
 * (e.g. "AV", "CD") used to set differentiated inventory control
 * policies (review frequency, safety stock, criticality of stockouts).
 */

#include <string>
#include <vector>

namespace mro::classification {

/** @brief Annual consumption value category (Pareto-based). */
enum class AbcCategory {
    A = 1, ///< High annual consumption value (top of the Pareto curve).
    B = 2, ///< Medium annual consumption value.
    C = 3  ///< Low annual consumption value.
};

/** @brief Operational criticality category (qualitative, user-supplied). */
enum class VedCategory {
    Vital     = 1, ///< A stockout has severe operational/safety impact.
    Essential = 2, ///< A stockout has moderate impact.
    Desirable = 3  ///< A stockout has minor or no operational impact.
};

/**
 * @brief Cumulative-value cut-off points (as fractions of total annual
 *        consumption value, in (0, 1]) used to separate A, B and C items.
 *
 * Defaults follow the common 80/15/5 Pareto convention:
 *   A: items whose cumulative value share (in descending-value rank
 *      order) falls up to a_cutoff (default 0.80).
 *   B: cumulative share from a_cutoff up to b_cutoff (default 0.95).
 *   C: cumulative share from b_cutoff to 1.0.
 */
struct AbcThresholds {
    double a_cutoff = 0.80;
    double b_cutoff = 0.95;
};

/**
 * @brief Per-item result of the ABC classification.
 */
struct AbcResult {
    AbcCategory category;         ///< A, B or C.
    double      annual_value;     ///< This item's annual consumption value.
    double      value_share;      ///< This item's share of total value (0-1).
    double      cumulative_share; ///< Cumulative share up to and including this item
                                   ///< (in descending-value rank order).
};

/**
 * @brief Classifies a portfolio of items into A/B/C categories based on
 *        their annual consumption value (e.g. unit_price * annual_demand).
 *
 * Items are ranked in descending order of annual_value; an item's
 * category is decided by where its cumulative share of total value
 * (accumulated in that descending order) falls relative to
 * @p thresholds. Ties in value are broken by input order (stable sort).
 *
 * Note: if a single item's value share alone already exceeds a_cutoff,
 * the strict cumulative rule below can place it in B (or even C) instead
 * of A. This is standard behavior for the cumulative-percentage method
 * taught in most inventory-management texts; if you need "largest item
 * is always A" semantics instead, special-case rank 0 at the call site.
 *
 * @param annual_values Annual consumption value of each item (must be >= 0).
 * @param thresholds    Cumulative-value cut-offs for A/B/C (see AbcThresholds).
 * @return One AbcResult per item, in the SAME ORDER as @p annual_values
 *         (i.e. result[i] corresponds to annual_values[i], not to the
 *         internal descending-value rank order).
 * @throws std::invalid_argument if @p annual_values is empty, contains a
 *         negative value, if a_cutoff/b_cutoff are not within (0, 1], or
 *         if a_cutoff >= b_cutoff.
 */
std::vector<AbcResult> classify_by_abc(
    const std::vector<double>& annual_values,
    const AbcThresholds& thresholds = AbcThresholds{});

/**
 * @brief Combines an ABC category with a (user-supplied) VED criticality
 *        rating into a single two-letter prioritization code.
 *
 * Examples: AbcCategory::A + VedCategory::Vital -> "AV";
 *           AbcCategory::C + VedCategory::Desirable -> "CD".
 *
 * @param abc ABC category for the item (from classify_by_abc()).
 * @param ved VED category for the item (supplied by the caller; see
 *            the file-level note on why this is not computed here).
 * @return Two-letter prioritization code, e.g. "AV", "BE", "CD".
 */
std::string combine_abc_ved(AbcCategory abc, VedCategory ved);

} // namespace mro::classification
