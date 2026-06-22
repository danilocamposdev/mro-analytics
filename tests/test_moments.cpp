#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include "mro/initialization/seed_generators/moments.hpp"

using mro::initialization::seed_generators::calculate_by_moments;
using mro::initialization::seed_generators::MomentsParameters;
using Catch::Approx;

TEST_CASE("calculate_by_moments - validacao de argumentos", "[moments][validation]") {
    SECTION("serie vazia lanca exception") {
        std::vector<double> empty{};
        REQUIRE_THROWS_AS(calculate_by_moments(empty, 2), std::invalid_argument);
    }

    SECTION("initial_periods <= 1 lanca exception") {
        std::vector<double> ts{1.0, 2.0, 3.0};
        REQUIRE_THROWS_AS(calculate_by_moments(ts, 1), std::invalid_argument);
        REQUIRE_THROWS_AS(calculate_by_moments(ts, 0), std::invalid_argument);
    }

    SECTION("initial_periods maior que o tamanho da serie lanca exception") {
        std::vector<double> ts{1.0, 2.0, 3.0};
        REQUIRE_THROWS_AS(calculate_by_moments(ts, 4), std::invalid_argument);
    }

    SECTION("initial_periods igual ao tamanho da serie e valido") {
        std::vector<double> ts{1.0, 2.0, 3.0};
        REQUIRE_NOTHROW(calculate_by_moments(ts, 3));
    }
}

TEST_CASE("calculate_by_moments - calculo correto de z1 e p1", "[moments][correctness]") {
    SECTION("janela com demandas mistas (zeros e nao-zeros)") {
        // Janela de inicializacao: {0, 0, 3, 0}
        // nonzero = 1, sum = 3 -> z1 = 3 / 1 = 3
        // prob1 = nonzero / initial_periods = 1 / 4 = 0.25
        // p1 = 1 / prob1 = 4
        std::vector<double> ts{0.0, 0.0, 3.0, 0.0, 5.0, 7.0};
        MomentsParameters result = calculate_by_moments(ts, 4);

        REQUIRE(result.z1 == Approx(3.0));
        REQUIRE(result.prob1 == Approx(0.25));
        REQUIRE(result.p1 == Approx(4.0));
    }

    SECTION("janela totalmente nao-zero") {
        // Janela: {2, 4, 6} -> nonzero = 3, sum = 12 -> z1 = 4
        // prob1 = 3/3 = 1 -> p1 = 1
        std::vector<double> ts{2.0, 4.0, 6.0, 0.0, 0.0};
        MomentsParameters result = calculate_by_moments(ts, 3);

        REQUIRE(result.z1 == Approx(4.0));
        REQUIRE(result.prob1 == Approx(1.0));
        REQUIRE(result.p1 == Approx(1.0));
    }

    SECTION("apenas a janela inicial e considerada, nao a serie inteira") {
        // Janela: {1, 1} -> z1 = 1, prob1 = 1, p1 = 1
        // O valor 100 fora da janela NAO deve influenciar o resultado.
        std::vector<double> ts{1.0, 1.0, 100.0, 100.0, 100.0};
        MomentsParameters result = calculate_by_moments(ts, 2);

        REQUIRE(result.z1 == Approx(1.0));
        REQUIRE(result.p1 == Approx(1.0));
    }

    SECTION("janela totalmente zero - caso degenerado") {
        // nonzero = 0 -> z1 = 0 (definido explicitamente no codigo)
        // prob1 = 0 -> p1 = 1/0 = +inf
        std::vector<double> ts{0.0, 0.0, 0.0, 5.0};
        MomentsParameters result = calculate_by_moments(ts, 3);

        REQUIRE(result.z1 == Approx(0.0));
        REQUIRE(result.prob1 == Approx(0.0));
        REQUIRE(std::isinf(result.p1));
    }
}
