#include "mro/forecasting/distributions/zip.hpp"

namespace mro::forecasting::distributions {

	double ZipResult::point_forecast() const {
		return (1.0 - pi) * lambda;
	}

	double ZipResult::demand_variance() const {
		return (1.0 - pi) * lambda * (1.0 + pi * lambda);
	}

	ZipResult forecast_by_zip(
			const mro::initialization::parameter_estimators::ZipParameters& params)
	{
		ZipResult result;
		result.pi = params.pi;
		result.lambda = params.lambda;
		return result;
	}

} // namespace mro::forecasting::distributions
