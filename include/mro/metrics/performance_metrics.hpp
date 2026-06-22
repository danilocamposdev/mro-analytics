#pragma once

/**
 * @file performance_metrics.hpp
 * @brief Inventory performance metrics for MRO / slow-moving items.
 *
 * Implements the standard inventory-control metrics used to evaluate
 * forecasting and replenishment policies under the classical
 * (continuous review, Normally distributed lead-time demand)
 * assumption: Reorder Point (ROP), Safety Stock (SS), Expected
 * Shortage per Cycle (ESC), Fill Rate, Cycle Service Level (CSL) and
 * Total Inventory Cost (TC), the latter including a routine to find
 * the order quantity that minimizes it.
 *
 * Reference:
 *   Silver, E.A., Pyke, D.F., Thomas, D.J. (2017). Inventory and
 *   Production Management in Supply Chains, 4th ed. CRC Press.
 */

#include "mro/forecasting/types/forecast_result.hpp"

namespace mro::metrics {

	/**
	 * @brief Breakdown of the components of the total inventory cost.
	 */
	struct TotalCostResult {
		double ordering_cost;  ///< Annual ordering cost = (D / Q) * K.
		double holding_cost;   ///< Annual holding cost = (Q / 2 + SS) * h.
		double shortage_cost;  ///< Annual shortage cost = (D / Q) * ESC * p.
		double total_cost;     ///< Sum of the three components above.
	};

	// ---------------------------------------------------------------
	// Standard normal helpers (lead-time demand is assumed Normal).
	// ---------------------------------------------------------------

	/**
	 * @brief Standard normal probability density function, phi(z).
	 */
	double standard_normal_pdf(double z);

	/**
	 * @brief Standard normal cumulative distribution function, Phi(z).
	 */
	double standard_normal_cdf(double z);

	/**
	 * @brief Standard normal unit loss function, L(z) = phi(z) - z*(1 - Phi(z)).
	 *
	 * Used to compute the expected shortage per replenishment cycle.
	 */
	double unit_normal_loss(double z);

	/**
	 * @brief Safety factor z associated with a target cycle service level.
	 *
	 * Inverts Phi(z) = csl using StatsLib's qnorm (quantile function of
	 * the standard normal distribution).
	 *
	 * @param csl Target cycle service level, 0 < csl < 1.
	 * @throws std::invalid_argument if csl is not in (0, 1).
	 */
	double z_for_service_level(double csl);

	// ---------------------------------------------------------------
	// Inventory metrics
	// ---------------------------------------------------------------

	/**
	 * @brief Safety stock, SS = z * sigma_L.
	 *
	 * @param sigma_lead_time_demand Standard deviation of demand during lead time.
	 * @param z Safety factor (see z_for_service_level()).
	 * @throws std::invalid_argument if sigma_lead_time_demand < 0.
	 */
	double safety_stock(double sigma_lead_time_demand, double z);

	/**
	 * @brief Reorder point, ROP = mu_L + SS.
	 *
	 * @param mean_lead_time_demand Expected demand during lead time, mu_L.
	 * @param safety_stock_value Safety stock, e.g. from safety_stock().
	 */
	double reorder_point(double mean_lead_time_demand, double safety_stock_value);

	/**
	 * @brief Expected shortage per replenishment cycle, ESC = sigma_L * L(z).
	 *
	 * @param sigma_lead_time_demand Standard deviation of demand during lead time.
	 * @param z Safety factor.
	 * @throws std::invalid_argument if sigma_lead_time_demand < 0.
	 */
	double expected_shortage_per_cycle(double sigma_lead_time_demand, double z);

	/**
	 * @brief Fill rate, fraction of demand met directly from stock on hand:
	 *        FR = 1 - ESC / Q.
	 *
	 * @param order_quantity Replenishment order quantity, Q (> 0).
	 * @param expected_shortage_per_cycle_value ESC, e.g. from expected_shortage_per_cycle().
	 * @throws std::invalid_argument if order_quantity <= 0.
	 */
	double fill_rate(double order_quantity, double expected_shortage_per_cycle_value);

	/**
	 * @brief Cycle service level, CSL = Phi(z): probability of not
	 *        stocking out during a single replenishment cycle.
	 */
	double cycle_service_level(double z);

	/**
	 * @brief Total expected annual inventory cost (EOQ extended with a
	 *        shortage/stockout cost term).
	 *
	 * TC(Q, SS) = (D/Q)*K + (Q/2 + SS)*h + (D/Q)*ESC*p
	 *
	 * @param annual_demand Annual demand, D (> 0).
	 * @param ordering_cost_per_order Fixed cost per order/replenishment, K (>= 0).
	 * @param order_quantity Replenishment order quantity, Q (> 0).
	 * @param holding_cost_per_unit Holding cost per unit per year, h (>= 0).
	 * @param safety_stock_value Safety stock, SS (>= 0).
	 * @param shortage_cost_per_unit Shortage/stockout cost per unit short, p (>= 0).
	 * @param expected_shortage_per_cycle_value Expected shortage per cycle, ESC (>= 0).
	 * @throws std::invalid_argument if annual_demand <= 0 or order_quantity <= 0.
	 */
	TotalCostResult total_inventory_cost(
			double annual_demand,
			double ordering_cost_per_order,
			double order_quantity,
			double holding_cost_per_unit,
			double safety_stock_value,
			double shortage_cost_per_unit,
			double expected_shortage_per_cycle_value);

	/**
	 * @brief Finds the order quantity Q in [q_min, q_max] that minimizes
	 *        total_inventory_cost(), holding all other parameters fixed.
	 *
	 * Total cost is convex in Q for fixed safety stock, so a
	 * dependency-free golden-section search is used.
	 *
	 * @param annual_demand Annual demand, D (> 0).
	 * @param ordering_cost_per_order Fixed cost per order, K (>= 0).
	 * @param holding_cost_per_unit Holding cost per unit per year, h (> 0).
	 * @param safety_stock_value Safety stock, SS (>= 0).
	 * @param shortage_cost_per_unit Shortage cost per unit short, p (>= 0).
	 * @param expected_shortage_per_cycle_value Expected shortage per cycle, ESC (>= 0).
	 * @param q_min Lower bound of the search range for Q (> 0).
	 * @param q_max Upper bound of the search range for Q (> q_min).
	 * @param tolerance Convergence tolerance on Q (> 0).
	 * @return The TotalCostResult evaluated at the minimizing Q.
	 * @throws std::invalid_argument if q_min <= 0, q_max <= q_min, or
	 * holding_cost_per_unit <= 0.
	 */
	TotalCostResult minimize_total_inventory_cost(
			double annual_demand,
			double ordering_cost_per_order,
			double holding_cost_per_unit,
			double safety_stock_value,
			double shortage_cost_per_unit,
			double expected_shortage_per_cycle_value,
			double q_min,
			double q_max,
			double tolerance = 1e-6);

	// ---------------------------------------------------------------
	// Pipeline: ForecastResult -> inventory policy
	// ---------------------------------------------------------------

	/**
	 * @brief Inventory policy parameters derived from a forecasting result.
	 */
	struct InventoryPolicyResult {
		double mu_L;                       ///< Mean demand during the lead time.
		double sigma_L;                    ///< Std. dev. of demand during the lead time.
		double z;                          ///< Safety factor associated with target_csl.
		double safety_stock;               ///< SS = z * sigma_L.
		double reorder_point;              ///< ROP = mu_L + SS.
		double expected_shortage_per_cycle;///< ESC = sigma_L * L(z).
		double cycle_service_level;        ///< CSL = Phi(z) (equals target_csl, returned for convenience).
	};

	/**
	 * @brief Full pipeline from any forecasting method's result to
	 *        inventory policy parameters (ROP, SS, ESC, CSL).
	 *
	 * Works with any model implementing mro::forecasting::ForecastResult
	 * (e.g. CrostonResult, and in the future SbaResult, TsbResult,
	 * SbjResult, ...), so this function never needs to know which
	 * forecasting method produced the result.
	 *
	 * @param result Forecasting result, e.g. from forecast_by_croston().
	 * @param lead_time_periods Lead time, in the same period unit as the
	 * demand series (e.g. quarters), > 0.
	 * @param target_csl Target cycle service level, 0 < target_csl < 1.
	 * @throws std::invalid_argument if lead_time_periods <= 0 or
	 * target_csl is not in (0, 1).
	 */
	InventoryPolicyResult evaluate_policy(
			const mro::forecasting::types::ForecastResult& result,
			double lead_time_periods,
			double target_csl);

} // namespace mro::metrics
