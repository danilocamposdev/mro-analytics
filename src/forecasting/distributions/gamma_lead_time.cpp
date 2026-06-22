#include "mro/forecasting/distributions/gamma_lead_time.hpp"

#include <cmath>

namespace mro::forecasting::distributions {

	mro::forecasting::types::LeadTimeDemand GammaLeadTimeResult::lead_time_demand() const {
		return mro::forecasting::types::LeadTimeDemand{
			shape * scale,
			std::sqrt(shape) * scale
		};
	}

	GammaLeadTimeResult forecast_by_gamma_lead_time(
			const mro::initialization::parameter_estimators::GammaParameters& params)
	{
		GammaLeadTimeResult result;
		result.shape = params.shape;
		result.scale = params.scale;

		return result;
	}

} // namespace mro::forecasting::distributions
