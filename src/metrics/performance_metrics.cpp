#include "mro/metrics/performance_metrics.hpp"

#include <stats.hpp>

#include <stdexcept>

namespace mro::metrics {

	double standard_normal_pdf(double z) {
		return stats::dnorm(z, 0.0, 1.0);
	}

	double standard_normal_cdf(double z) {
		return stats::pnorm(z, 0.0, 1.0);
	}

	double unit_normal_loss(double z) {
		return standard_normal_pdf(z) - z * (1.0 - standard_normal_cdf(z));
	}

	double z_for_service_level(double csl) {
		if (csl <= 0.0 || csl >= 1.0) {
			throw std::invalid_argument(
					"z_for_service_level: csl must be in (0, 1)");
		}
		return stats::qnorm(csl, 0.0, 1.0);
	}

	double safety_stock(double sigma_lead_time_demand, double z) {
		if (sigma_lead_time_demand < 0.0) {
			throw std::invalid_argument(
					"safety_stock: sigma_lead_time_demand must be >= 0");
		}
		return z * sigma_lead_time_demand;
	}

	double reorder_point(double mean_lead_time_demand, double safety_stock_value) {
		return mean_lead_time_demand + safety_stock_value;
	}

	double expected_shortage_per_cycle(double sigma_lead_time_demand, double z) {
		if (sigma_lead_time_demand < 0.0) {
			throw std::invalid_argument(
					"expected_shortage_per_cycle: sigma_lead_time_demand must be >= 0");
		}
		return sigma_lead_time_demand * unit_normal_loss(z);
	}

	double fill_rate(double order_quantity, double expected_shortage_per_cycle_value) {
		if (order_quantity <= 0.0) {
			throw std::invalid_argument("fill_rate: order_quantity must be > 0");
		}
		return 1.0 - expected_shortage_per_cycle_value / order_quantity;
	}

	double cycle_service_level(double z) {
		return standard_normal_cdf(z);
	}

	TotalCostResult total_inventory_cost(
			double annual_demand,
			double ordering_cost_per_order,
			double order_quantity,
			double holding_cost_per_unit,
			double safety_stock_value,
			double shortage_cost_per_unit,
			double expected_shortage_per_cycle_value) {
		if (annual_demand <= 0.0) {
			throw std::invalid_argument(
					"total_inventory_cost: annual_demand must be > 0");
		}
		if (order_quantity <= 0.0) {
			throw std::invalid_argument(
					"total_inventory_cost: order_quantity must be > 0");
		}

		const double orders_per_year = annual_demand / order_quantity;

		TotalCostResult result{};
		result.ordering_cost = orders_per_year * ordering_cost_per_order;
		result.holding_cost =
				(order_quantity / 2.0 + safety_stock_value) * holding_cost_per_unit;
		result.shortage_cost =
				orders_per_year * expected_shortage_per_cycle_value * shortage_cost_per_unit;
		result.total_cost =
				result.ordering_cost + result.holding_cost + result.shortage_cost;

		return result;
	}

	TotalCostResult minimize_total_inventory_cost(
			double annual_demand,
			double ordering_cost_per_order,
			double holding_cost_per_unit,
			double safety_stock_value,
			double shortage_cost_per_unit,
			double expected_shortage_per_cycle_value,
			double q_min,
			double q_max,
			double tolerance) {
		if (q_min <= 0.0) {
			throw std::invalid_argument(
					"minimize_total_inventory_cost: q_min must be > 0");
		}
		if (q_max <= q_min) {
			throw std::invalid_argument(
					"minimize_total_inventory_cost: q_max must be > q_min");
		}
		if (holding_cost_per_unit <= 0.0) {
			throw std::invalid_argument(
					"minimize_total_inventory_cost: holding_cost_per_unit must be > 0");
		}

		auto cost_at = [&](double q) {
			return total_inventory_cost(
					annual_demand, ordering_cost_per_order, q, holding_cost_per_unit,
					safety_stock_value, shortage_cost_per_unit,
					expected_shortage_per_cycle_value).total_cost;
		};

		constexpr double kInvPhi = 0.6180339887498949; // (sqrt(5) - 1) / 2

		double a = q_min;
		double b = q_max;
		double c = b - kInvPhi * (b - a);
		double d = a + kInvPhi * (b - a);
		double fc = cost_at(c);
		double fd = cost_at(d);

		while (b - a > tolerance) {
			if (fc < fd) {
				b = d;
				d = c;
				fd = fc;
				c = b - kInvPhi * (b - a);
				fc = cost_at(c);
			} else {
				a = c;
				c = d;
				fc = fd;
				d = a + kInvPhi * (b - a);
				fd = cost_at(d);
			}
		}

		const double best_q = (a + b) / 2.0;
		return total_inventory_cost(
				annual_demand, ordering_cost_per_order, best_q, holding_cost_per_unit,
				safety_stock_value, shortage_cost_per_unit,
				expected_shortage_per_cycle_value);
	}

	InventoryPolicyResult evaluate_policy(
			const mro::forecasting::types::ForecastResult& result,
			double lead_time_periods,
			double target_csl) {
		if (lead_time_periods <= 0.0) {
			throw std::invalid_argument("evaluate_policy: lead_time_periods must be > 0");
		}

		const auto ltd = result.lead_time_demand(lead_time_periods);
		const double z = z_for_service_level(target_csl);

		InventoryPolicyResult policy{};
		policy.mu_L = ltd.mean;
		policy.sigma_L = ltd.std_dev;
		policy.z = z;
		policy.safety_stock = safety_stock(ltd.std_dev, z);
		policy.reorder_point = reorder_point(ltd.mean, policy.safety_stock);
		policy.expected_shortage_per_cycle = expected_shortage_per_cycle(ltd.std_dev, z);
		policy.cycle_service_level = cycle_service_level(z);

		return policy;
	}

	InventoryPolicyResult evaluate_policy(
			const mro::forecasting::types::LeadTimeDistributionResult& result,
			double target_csl) {
		const auto ltd = result.lead_time_demand();
		const double z = z_for_service_level(target_csl);

		InventoryPolicyResult policy{};
		policy.mu_L = ltd.mean;
		policy.sigma_L = ltd.std_dev;
		policy.z = z;
		policy.safety_stock = safety_stock(ltd.std_dev, z);
		policy.reorder_point = reorder_point(ltd.mean, policy.safety_stock);
		policy.expected_shortage_per_cycle = expected_shortage_per_cycle(ltd.std_dev, z);
		policy.cycle_service_level = cycle_service_level(z);

		return policy;
	}

} // namespace mro::metrics
