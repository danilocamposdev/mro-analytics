#include "mro/forecasting/distributions/negative_binomial.hpp"

namespace mro::forecasting::distributions {

	double NegativeBinomialResult::point_forecast() const {
		return mu;
	}

	double NegativeBinomialResult::demand_variance() const {
		return mu + (mu * mu) / r;
	}

	NegativeBinomialResult forecast_by_negative_binomial(
			const mro::initialization::parameter_estimators::NegativeBinomialParameters& params)
	{
		NegativeBinomialResult result;
		result.r = params.r;
		result.mu = params.mu;
		return result;
	}

} // namespace mro::forecasting::distributions
