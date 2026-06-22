#pragma once

/**
 * @file syntetos_boylan.hpp
 * @brief Syntetos, Boylan & Croston (2005) demand categorization scheme (SBC).
 *
 * Reference:
 *   Syntetos, A.A., Boylan, J.E. & Croston, J.D. (2005). On the
 *   categorization of demand patterns. Journal of the Operational
 *   Research Society, 56(5), 495-503.
 *
 * SBC is probably the most widely cited cut-off scheme in the
 * intermittent-demand literature. It uses a 2x2 ADI x CV² matrix, with
 * cut-offs (ADI = 1.32, CV² = 0.49) derived theoretically by comparing
 * the relative forecast-error performance of SES, Croston and SBA across
 * four resulting demand patterns: Smooth, Erratic, Intermittent, Lumpy.
 *
 * This method defines its own category enum and result struct; it does
 * NOT reuse types from williams.hpp or eaves.hpp, so each classification
 * method can evolve independently.
 */

#include <string>
#include <vector>

namespace mro::classification {

	/**
	 * @brief Demand pattern categories as defined by the Syntetos-Boylan-Croston
	 *        (SBC) scheme.
	 *
	 * @verbatim
	 *               CV² <= 0.49       CV² > 0.49
	 *  ADI <= 1.32   Smooth            Erratic
	 *  ADI >  1.32   Intermittent      Lumpy
	 * @endverbatim
	 */
	enum class SbcCategory {
		Smooth       = 1, ///< Frequent and stable.
		Erratic      = 2, ///< Frequent but volatile in quantity.
		Intermittent = 3, ///< Rare but stable in quantity.
		Lumpy        = 4  ///< Sparse and variable.
	};

	/** @brief Empirical cut-off values for the Syntetos-Boylan-Croston (SBC) scheme. */
	namespace sbc_thresholds {
		constexpr double ADI_INTERMITTENCY = 1.32;
		constexpr double CV2_LUMPINESS     = 0.49;
	}

	/**
	 * @brief Output of classify_by_syntetos_boylan().
	 *
	 * Carries the classification decision and the intermediate statistics
	 * that produced it.
	 */
	struct SbcClassificationResult {
		SbcCategory  category;      ///< Demand pattern class (enum).
		std::string  category_name; ///< Human-readable label (e.g. "Lumpy").
		double       adi;           ///< Average Demand Interval = 1/p.
		double       cv2;           ///< CV² used for classification.
		double       p;             ///< Demand probability per period.
		double       mean_nonzero;  ///< Mean demand size on non-zero periods only.
	};

	/**
	 * @brief Classifies a demand time series using the Syntetos, Boylan &
	 *        Croston (2005) categorization scheme (SBC).
	 *
	 * ADI and CV² are computed on non-zero demand periods only. An item with
	 * no non-zero demand is classified as Intermittent.
	 *
	 * @param time_series Chronological sequence of period demands (zeros allowed).
	 * @return SbcClassificationResult with category, label, ADI, CV², p, and mean_nonzero.
	 * @throws std::invalid_argument if @p time_series is empty.
	 */
	SbcClassificationResult classify_by_syntetos_boylan(const std::vector<double>& time_series);

} // namespace mro::classification
