#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "mro/forecasting/parametrics/croston.hpp"
#include "mro/initialization/types/initialization_results.hpp"
#include "mro/initialization/seed_generators/moments.hpp"

using mro::forecasting::parametrics::forecast_by_croston;
using mro::forecasting::parametrics::CrostonResult;
using mro::initialization::types::SeedGeneratorsResults;
using mro::initialization::seed_generators::calculate_by_moments;
using Catch::Approx;

namespace {
SeedGeneratorsResults make_init(double z1, double p1) {
    SeedGeneratorsResults init;
    init.z1 = z1;
    init.p1 = p1;
    return init;
}
} // namespace

TEST_CASE("forecast_by_croston - validacao de argumentos", "[croston][validation]") {
    SeedGeneratorsResults valid_init = make_init(2.0, 2.0);
    std::vector<double> ts{1.0, 0.0, 2.0};

    SECTION("serie vazia lanca exception") {
        std::vector<double> empty{};
        REQUIRE_THROWS_AS(forecast_by_croston(empty, valid_init, 0.1), std::invalid_argument);
    }

    SECTION("alpha fora do intervalo (0, 1] lanca exception") {
        REQUIRE_THROWS_AS(forecast_by_croston(ts, valid_init, 0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(forecast_by_croston(ts, valid_init, 1.1), std::invalid_argument);
    }

    SECTION("init.p1 nao positivo lanca exception") {
        SeedGeneratorsResults invalid_init = make_init(2.0, 0.0);
        REQUIRE_THROWS_AS(forecast_by_croston(ts, invalid_init, 0.1), std::invalid_argument);
    }

    // --- NOVA VALIDAÇÃO ADICIONADA ---
    SECTION("skip_periods maior ou igual ao tamanho da serie lanca exception") {
        REQUIRE_THROWS_AS(forecast_by_croston(ts, valid_init, 0.1, 3), std::invalid_argument);
        REQUIRE_THROWS_AS(forecast_by_croston(ts, valid_init, 0.1, 5), std::invalid_argument);
    }
}

TEST_CASE("forecast_by_croston - comportamento basico sem skip_periods", "[croston][core]") {
    std::vector<double> ts{0.0, 0.0, 6.0, 0.0, 0.0, 9.0};
    SeedGeneratorsResults init = make_init(100.0, 100.0);

    // Sem passar skip_periods, assume 0 (comportamento legado mantido)
    CrostonResult result = forecast_by_croston(ts, init, 1.0);

    REQUIRE(result.z_hat == Approx(9.0));
    REQUIRE(result.p_hat == Approx(3.0));
    REQUIRE(result.forecast == Approx(3.0));
}

TEST_CASE("forecast_by_croston - integracao com calculate_by_moments e skip_periods", "[croston][integration]") {
    // Série completa com 8 períodos
    std::vector<double> ts{0.0, 0.0, 3.0, 0.0, 0.0, 2.0, 0.0, 4.0};
    size_t initial_periods = 4; // Usados na inicialização: {0, 0, 3, 0}

    // Inicializacao (moments): nonzero=1, sum=3 -> z1=3
    // prob1 = 1/4 = 0.25 -> p1 = 1/0.25 = 4
    auto init = calculate_by_moments(ts, initial_periods);
    
    REQUIRE(init.z1 == Approx(3.0));
    REQUIRE(init.p1 == Approx(4.0));

    // Agora o Croston pula os primeiros 4 períodos e avalia apenas o restante: {0, 2, 0, 4}
    // O loop começa no índice 4 (valor 0.0).
    // i=4 (val=0): periods_since_demand = 1. z_hat e p_hat continuam {3.0, 4.0}
    // i=5 (val=2): periods_since_demand = 2. Como val > 0:
    //              z_hat = 0.2 * 2 + 0.8 * 3.0 = 2.8
    //              p_hat = 0.2 * 2 + 0.8 * 4.0 = 3.6
    //              periods_since_demand = 0
    // i=6 (val=0): periods_since_demand = 1. Sem alterações.
    // i=7 (val=4): periods_since_demand = 2. Como val > 0:
    //              z_hat = 0.2 * 4 + 0.8 * 2.8 = 3.04
    //              p_hat = 0.2 * 2 + 0.8 * 3.6 = 3.28
    CrostonResult result = forecast_by_croston(ts, init, 0.2, initial_periods);

    REQUIRE(result.z_hat == Approx(3.04).margin(0.001));
    REQUIRE(result.p_hat == Approx(3.28).margin(0.001));
    REQUIRE(result.forecast == Approx(3.04 / 3.28).margin(0.001));
}
