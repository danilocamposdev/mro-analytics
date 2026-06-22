#include "mro/forecasting/non_parametrics/willemain_bootstrap.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

namespace mro::forecasting::non_parametrics {

	namespace {

		/// Two-state (zero / non-zero) Markov transition matrix, estimated
		/// by counting observed transitions in the occurrence indicator.
		struct TransitionMatrix {
			double p00;  ///< P(occurrence_{t+1} = 0 | occurrence_t = 0)
			double p01;  ///< P(occurrence_{t+1} = 1 | occurrence_t = 0)
			double p10;  ///< P(occurrence_{t+1} = 0 | occurrence_t = 1)
			double p11;  ///< P(occurrence_{t+1} = 1 | occurrence_t = 1)
		};

		/// Estimates the two-state Markov transition matrix from the
		/// occurrence indicator of time_series, by counting 0->0, 0->1,
		/// 1->0 and 1->1 transitions between consecutive periods.
		TransitionMatrix estimate_transition_matrix(const std::vector<double>& time_series) {
			double count_0_0 = 0.0;
			double count_0_1 = 0.0;
			double count_1_0 = 0.0;
			double count_1_1 = 0.0;

			for (std::size_t i = 0; i + 1 < time_series.size(); ++i) {
				const bool from_nonzero = time_series[i] != 0.0;
				const bool to_nonzero = time_series[i + 1] != 0.0;

				if (!from_nonzero && !to_nonzero) {
					count_0_0 += 1.0;
				} else if (!from_nonzero && to_nonzero) {
					count_0_1 += 1.0;
				} else if (from_nonzero && !to_nonzero) {
					count_1_0 += 1.0;
				} else {
					count_1_1 += 1.0;
				}
			}

			const double total_from_0 = count_0_0 + count_0_1;
			const double total_from_1 = count_1_0 + count_1_1;

			TransitionMatrix matrix{};
			// If a state was never visited (no transitions out of it observed),
			// fall back to staying in the same state with probability 1: there
			// is no data to suggest otherwise, and this keeps the chain
			// well-defined without injecting an arbitrary 50/50 guess.
			matrix.p00 = (total_from_0 > 0.0) ? (count_0_0 / total_from_0) : 1.0;
			matrix.p01 = (total_from_0 > 0.0) ? (count_0_1 / total_from_0) : 0.0;
			matrix.p10 = (total_from_1 > 0.0) ? (count_1_0 / total_from_1) : 0.0;
			matrix.p11 = (total_from_1 > 0.0) ? (count_1_1 / total_from_1) : 1.0;

			return matrix;
		}

		/// Sample standard deviation (denominator n, matching the course
		/// material's sigma_Z computation) of a non-empty vector.
		double sample_std_dev(const std::vector<double>& values) {
			const double n = static_cast<double>(values.size());
			const double mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
			double sum_sq_dev = 0.0;
			for (double val : values) {
				sum_sq_dev += (val - mean) * (val - mean);
			}
			return std::sqrt(sum_sq_dev / n);
		}

		/// Interquartile range using the midpoint-of-halves convention
		/// from the course material (Q1 = median of the lower half, Q3 =
		/// median of the upper half, splitting evenly with no shared
		/// midpoint element), which matches the worked IQR=11 example for
		/// Z=[5,6,7,8,12,15,20,25].
		double interquartile_range(std::vector<double> values) {
			std::sort(values.begin(), values.end());
			const std::size_t n = values.size();
			const std::size_t half = n / 2;

			auto median_of = [](const std::vector<double>& v, std::size_t begin, std::size_t end) {
				const std::size_t len = end - begin;
				if (len % 2 == 1) {
					return v[begin + len / 2];
				}
				return (v[begin + len / 2 - 1] + v[begin + len / 2]) / 2.0;
			};

			const double q1 = median_of(values, 0, half);
			const double q3 = median_of(values, n - half, n);
			return q3 - q1;
		}

		/// Jitter half-width h for the chosen strategy.
		double compute_jitter_half_width(
				const std::vector<double>& demand_sizes,
				JitterStrategy strategy,
				double jitter_alpha) {
			switch (strategy) {
				case JitterStrategy::kAlphaSigma:
					return jitter_alpha * sample_std_dev(demand_sizes);
				case JitterStrategy::kIqr:
					return jitter_alpha * interquartile_range(demand_sizes);
				case JitterStrategy::kSilverman: {
					const double sigma = sample_std_dev(demand_sizes);
					const double iqr_term = interquartile_range(demand_sizes) / 1.34;
					const double n = static_cast<double>(demand_sizes.size());
					return 0.9 * std::min(sigma, iqr_term) * std::pow(n, -0.2);
				}
			}
			return jitter_alpha * sample_std_dev(demand_sizes); // unreachable
		}

	} // namespace

	double WillemainBootstrapResult::point_forecast() const {
		const double n = static_cast<double>(simulated_lead_time_demand.size());
		const double sum = std::accumulate(
				simulated_lead_time_demand.begin(), simulated_lead_time_demand.end(), 0.0);
		return (sum / n) / lead_time_periods;
	}

	double WillemainBootstrapResult::demand_variance() const {
		const double n = static_cast<double>(simulated_lead_time_demand.size());
		const double mean = std::accumulate(
				simulated_lead_time_demand.begin(), simulated_lead_time_demand.end(), 0.0) / n;
		double sum_sq_dev = 0.0;
		for (double val : simulated_lead_time_demand) {
			sum_sq_dev += (val - mean) * (val - mean);
		}
		return (sum_sq_dev / n) / lead_time_periods;
	}

	mro::forecasting::types::LeadTimeDemand WillemainBootstrapResult::lead_time_demand(
			double lead_time_periods_arg) const {
		if (std::abs(lead_time_periods_arg - lead_time_periods) > 1e-9) {
			throw std::invalid_argument(
					"WillemainBootstrapResult::lead_time_demand: lead_time_periods does "
					"not match the horizon this result was simulated for; re-run "
					"forecast_by_willemain_bootstrap() for the new horizon instead of "
					"rescaling, since the Markov-driven occurrence sequence is specific "
					"to the simulated horizon.");
		}

		const double n = static_cast<double>(simulated_lead_time_demand.size());
		const double mean = std::accumulate(
				simulated_lead_time_demand.begin(), simulated_lead_time_demand.end(), 0.0) / n;
		double sum_sq_dev = 0.0;
		for (double val : simulated_lead_time_demand) {
			sum_sq_dev += (val - mean) * (val - mean);
		}
		const double std_dev = std::sqrt(sum_sq_dev / n);

		return mro::forecasting::types::LeadTimeDemand{mean, std_dev};
	}

	double WillemainBootstrapResult::empirical_quantile(double probability) const {
		if (probability < 0.0 || probability > 1.0) {
			throw std::invalid_argument(
					"empirical_quantile: probability must be in [0, 1].");
		}
		if (simulated_lead_time_demand.empty()) {
			throw std::invalid_argument(
					"empirical_quantile: simulated_lead_time_demand is empty.");
		}

		std::vector<double> sorted_values = simulated_lead_time_demand;
		std::sort(sorted_values.begin(), sorted_values.end());

		const double n = static_cast<double>(sorted_values.size());
		// Nearest-rank method: e.g. P95 of 10,000 simulations is the
		// 9,500th smallest value (matching the course material's "500th
		// largest of 10,000" for P95).
		std::size_t rank = static_cast<std::size_t>(std::ceil(probability * n));
		if (rank == 0) {
			rank = 1;
		}
		if (rank > sorted_values.size()) {
			rank = sorted_values.size();
		}
		return sorted_values[rank - 1];
	}

	WillemainBootstrapResult forecast_by_willemain_bootstrap(
			const std::vector<double>& time_series,
			double lead_time_periods,
			std::size_t n_simulations,
			JitterStrategy jitter_strategy,
			double jitter_alpha,
			unsigned int random_seed) {
		if (time_series.empty()) {
			throw std::invalid_argument("Time series cannot be empty.");
		}
		if (time_series.size() < 2) {
			throw std::invalid_argument(
					"forecast_by_willemain_bootstrap: at least 2 observations are "
					"required to estimate the Markov transition matrix.");
		}
		if (lead_time_periods <= 0.0) {
			throw std::invalid_argument(
					"forecast_by_willemain_bootstrap: lead_time_periods must be > 0.");
		}
		if (n_simulations == 0) {
			throw std::invalid_argument(
					"forecast_by_willemain_bootstrap: n_simulations must be > 0.");
		}
		if (jitter_alpha < 0.0) {
			throw std::invalid_argument(
					"forecast_by_willemain_bootstrap: jitter_alpha must be >= 0.");
		}

		std::vector<double> demand_sizes;
		for (double val : time_series) {
			if (val != 0.0) {
				demand_sizes.push_back(val);
			}
		}
		if (demand_sizes.size() < 2) {
			throw std::invalid_argument(
					"forecast_by_willemain_bootstrap: at least 2 non-zero demands are "
					"required to bootstrap demand sizes.");
		}

		const TransitionMatrix transition = estimate_transition_matrix(time_series);
		const double jitter_half_width =
				compute_jitter_half_width(demand_sizes, jitter_strategy, jitter_alpha);

		// Full periods to simulate, plus a fractional final period (if any)
		// whose simulated size, when it occurs, is scaled by the remainder.
		const std::size_t full_periods = static_cast<std::size_t>(std::floor(lead_time_periods));
		const double fractional_period = lead_time_periods - static_cast<double>(full_periods);
		const bool has_fractional_period = fractional_period > 1e-9;

		const bool last_observed_nonzero = time_series.back() != 0.0;

		std::random_device rd;
		std::mt19937 rng(random_seed != 0 ? random_seed : rd());
		std::uniform_real_distribution<double> unit_uniform(0.0, 1.0);
		std::uniform_real_distribution<double> jitter_uniform(
				-jitter_half_width, jitter_half_width);
		std::uniform_int_distribution<std::size_t> size_index(0, demand_sizes.size() - 1);

		std::vector<double> simulated(n_simulations);

		for (std::size_t sim = 0; sim < n_simulations; ++sim) {
			double lead_time_demand_draw = 0.0;
			bool occurrence = last_observed_nonzero;

			auto simulate_period = [&]() -> double {
				const double p_occurrence = occurrence ? transition.p11 : transition.p01;
				occurrence = unit_uniform(rng) < p_occurrence;

				if (!occurrence) {
					return 0.0;
				}

				const double resampled_size = demand_sizes[size_index(rng)];
				return resampled_size + jitter_uniform(rng);
			};

			for (std::size_t period = 0; period < full_periods; ++period) {
				lead_time_demand_draw += simulate_period();
			}

			if (has_fractional_period) {
				lead_time_demand_draw += simulate_period() * fractional_period;
			}

			simulated[sim] = lead_time_demand_draw;
		}

		WillemainBootstrapResult result;
		result.simulated_lead_time_demand = std::move(simulated);
		result.lead_time_periods = lead_time_periods;

		return result;
	}

} // namespace mro::forecasting::non_parametrics
