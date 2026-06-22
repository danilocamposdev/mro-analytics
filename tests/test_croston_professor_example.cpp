#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "mro/forecasting/parametrics/croston.hpp"
#include "mro/initialization/types/initialization_results.hpp"
#include "mro/initialization/seed_generators/classic_ad_hoc.hpp"

using mro::forecasting::parametrics::forecast_by_croston;
using mro::forecasting::parametrics::CrostonResult;
using mro::initialization::types::SeedGeneratorsResults;
using mro::initialization::seed_generators::calculate_by_classic_ad_hoc;
using Catch::Approx;

TEST_CASE("forecast_by_croston - Exemplo dos Slides do Professor (Sem skip_periods)", 
          "[croston][professor][gabarito]") {
    
    const std::vector<double> full_series{
        37, 5, 0, 14, 5, 0, 10, 10, 0, 0, 6, 20, 32, 5, 25, 38,
        15, 6, 70, 0, 0, 0, 10, 0
    };

    // No slide o professor usa z1 = 18.667 e p1 = 1.50
    SeedGeneratorsResults init;
    init.z1 = 18.667;
    init.p1 = 1.50;

    const double alpha = 0.1;
    const double tol = 0.005;

    // Cenário idêntico à tabela: processa desde o índice 0 (skip_periods = 0)
    // t=1: val=37 -> z = 0.1*37 + 0.9*18.667 = 20.500; p = 0.1*1 + 0.9*1.5 = 1.450
    {
        auto res = forecast_by_croston({37}, init, alpha, 0);
        REQUIRE(res.z_hat == Approx(20.500).margin(tol));
        REQUIRE(res.p_hat == Approx(1.450).margin(tol));
    }
    // t=2: val=5 -> z = 0.1*5 + 0.9*20.500 = 18.950; p = 0.1*1 + 0.9*1.450 = 1.405
    {
        auto res = forecast_by_croston({37, 5}, init, alpha, 0);
        REQUIRE(res.z_hat == Approx(18.950).margin(tol));
        REQUIRE(res.p_hat == Approx(1.405).margin(tol));
    }
}

TEST_CASE("forecast_by_croston + classic_ad_hoc - Pipeline com skip_periods",
          "[croston][classic_ad_hoc][integration][skip]") {

    const std::vector<double> full_series{
        37, 5, 0, 14, 5, 0, 10, 10, 0, 0, 6, 20, 32, 5, 25, 38,
        15, 6, 70, 0, 0, 0, 10, 0
    };

    size_t initial_periods = 4; // Primeiro ano usado para inicialização

    // Gera as sementes: z1 = 18.667, p1 = 1.50
    auto init = calculate_by_classic_ad_hoc(full_series, initial_periods);
    REQUIRE(init.z1 == Approx(18.667).margin(0.001));
    REQUIRE(init.p1 == Approx(1.50));

    const double alpha = 0.1;

    // Roda o pipeline pulando o primeiro ano. O loop começará no trimestre 5 (índice 4 = 5.0)
    // sementes iniciais prontas para o t=5: z_hat = 18.667, p_hat = 1.50
    // i=4 (val=5): periods_since_demand = 1. Como val > 0:
    //              z_hat = 0.1 * 5 + 0.9 * 18.667 = 17.300
    //              p_hat = 0.1 * 1 + 0.9 * 1.50 = 1.450
    CrostonResult result = forecast_by_croston(full_series, init, alpha, initial_periods);

    REQUIRE(result.z_hat != Approx(18.667)); // Garante que houve evolução
    
    // O teste valida o primeiro passo pós-pulo diretamente na série completa
    auto single_step_check = forecast_by_croston({37, 5, 0, 14, 5}, init, alpha, initial_periods);
    REQUIRE(single_step_check.z_hat == Approx(17.300).margin(0.005));
    REQUIRE(single_step_check.p_hat == Approx(1.450).margin(0.005));
}
