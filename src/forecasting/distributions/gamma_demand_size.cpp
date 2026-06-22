#include "mro/forecasting/distributions/gamma_demand_size.hpp"

#include <stdexcept>

namespace mro::forecasting::distributions {

	double GammaDemandResult::point_forecast() const {
		return prob1 * (shape * scale);
	}

	double GammaDemandResult::demand_variance() const {
		const double mean_size = shape * scale;
		const double var_size = shape * scale * scale;
		// Variance of a Bernoulli(prob1)-gated random variable X:
		//   Var = prob1 * Var(size) + prob1 * (1 - prob1) * E[size]^2
		return prob1 * var_size + prob1 * (1.0 - prob1) * mean_size * mean_size;
	}

	GammaDemandResult forecast_by_gamma_demand_size(
			const mro::initialization::parameter_estimators::GammaParameters& params,
			double prob1)
	{
		if (prob1 <= 0.0 || prob1 > 1.0) {
			throw std::invalid_argument("prob1 must be in (0, 1].");
		}

		GammaDemandResult result;
		result.shape = params.shape;
		result.scale = params.scale;
		result.prob1 = prob1;

		return result;
	}

} // namespace mro::forecasting::distributions
