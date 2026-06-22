cmake_minimum_required(VERSION 3.24)
project(mro_analytics VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# --- GCEM (dependency of StatsLib) ---
FetchContent_Declare(
	gcem
	GIT_REPOSITORY https://github.com/kthohr/gcem.git
	GIT_TAG        v1.18.0
)
FetchContent_MakeAvailable(gcem)

# --- StatsLib (header-only statistical distributions) ---
FetchContent_Declare(
	stats
	GIT_REPOSITORY https://github.com/kthohr/stats.git
	GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(stats)

# --- Biblioteca core C++ ---
add_library(${PROJECT_NAME}
	src/classification/williams.cpp
	src/classification/eaves.cpp
    src/classification/syntetos_boylan.cpp
    src/classification/abc_ved.cpp
	src/initialization/seed_generators/moments.cpp
	src/initialization/seed_generators/classic_ad_hoc.cpp
	src/initialization/seed_generators/backcasting.cpp
	src/initialization/parameter_estimators/mle_zip.cpp
	src/initialization/parameter_estimators/mle_gamma.cpp
	src/initialization/parameter_estimators/mle_negative_binomial.cpp
	src/forecasting/parametrics/croston.cpp
	src/forecasting/parametrics/sba.cpp
	src/forecasting/parametrics/sbj.cpp
	src/forecasting/parametrics/tsb.cpp
	src/forecasting/distributions/zip.cpp
	src/forecasting/distributions/gamma_demand_size.cpp
	src/forecasting/distributions/gamma_lead_time.cpp
	src/forecasting/distributions/negative_binomial.cpp
	src/metrics/performance_metrics.cpp
	src/calibration/mase.cpp
	src/calibration/mad.cpp
	src/calibration/mse.cpp
	src/preprocessing/adjusted_boxplot.cpp
	src/preprocessing/negative_binomial_overdispersion.cpp
)

target_include_directories(${PROJECT_NAME}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${gcem_SOURCE_DIR}/include
	${stats_SOURCE_DIR}/include
)

add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# --- Testes ---
enable_testing()
add_subdirectory(tests)
