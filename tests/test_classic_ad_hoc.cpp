#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "mro/initialization/seed_generators/classic_ad_hoc.hpp"

using mro::initialization::seed_generators::calculate_by_classic_ad_hoc;
using mro::initialization::seed_generators::ClassicAdHocParameters;
using Catch::Approx;

TEST_CASE("calculate_by_classic_ad_hoc - validacao de argumentos", "[classic_ad_hoc][validation]") {
	SECTION("serie vazia lanca exception") {
		std::vector<double> empty{};
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(empty, 2), std::invalid_argument);
	}

	SECTION("initial_periods <= 1 lanca exception") {
		std::vector<double> ts{1.0, 2.0, 3.0};
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(ts, 1), std::invalid_argument);
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(ts, 0), std::invalid_argument);
	}

	SECTION("initial_periods maior que o tamanho da serie lanca exception") {
		std::vector<double> ts{1.0, 2.0, 3.0};
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(ts, 4), std::invalid_argument);
	}

	SECTION("janela sem nenhuma demanda nao-zero lanca exception (z1 indefinido)") {
		std::vector<double> ts{0.0, 0.0, 0.0, 5.0};
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(ts, 3), std::invalid_argument);
	}

	SECTION("janela com apenas uma demanda nao-zero lanca exception (sem intervalo para p1)") {
		// Apenas um valor positivo na janela: nao ha gap entre demandas
		// sucessivas para calcular p1.
		std::vector<double> ts{0.0, 7.0, 0.0, 0.0};
		REQUIRE_THROWS_AS(calculate_by_classic_ad_hoc(ts, 4), std::invalid_argument);
	}
}

TEST_CASE("calculate_by_classic_ad_hoc - reproduz o exemplo do professor", "[classic_ad_hoc][golden][professor]") {
	// Slides "Aula_Croston_SBA_SBJ", paginas 30-32.
	// Serie: [37, 5, 0, 14, ...], janela de inicializacao = 4 primeiros trimestres.
	//
	// Demandas nao-zero na janela: 37 (t1), 5 (t2), 14 (t4)
	//   z1 = (37 + 5 + 14) / 3 = 18.667
	//
	// Gaps entre demandas sucessivas: t1->t2 = 1 ; t2->t4 = 2
	//   p1 = (1 + 2) / 2 = 1.50   (media dos GAPS, nao das demandas)
	std::vector<double> ts{37, 5, 0, 14, 5, 0, 10, 10, 0, 0, 6, 20, 32, 5, 25, 38,
		15, 6, 70, 0, 0, 0, 10, 0};

	ClassicAdHocParameters result = calculate_by_classic_ad_hoc(ts, 4);

	REQUIRE(result.z1 == Approx(18.667).margin(0.001));
	REQUIRE(result.p1 == Approx(1.50));
}

TEST_CASE("calculate_by_classic_ad_hoc - calculo correto em outros cenarios", "[classic_ad_hoc][correctness]") {
	SECTION("demandas consecutivas sem zeros (gap sempre 1)") {
		// Janela: {2, 4, 6} -> nonzero=3, z1=4
		// Gaps: t1->t2=1, t2->t3=1 -> p1 = (1+1)/2 = 1
		std::vector<double> ts{2.0, 4.0, 6.0, 0.0};
		ClassicAdHocParameters result = calculate_by_classic_ad_hoc(ts, 3);

		REQUIRE(result.z1 == Approx(4.0));
		REQUIRE(result.p1 == Approx(1.0));
	}

	SECTION("apenas a janela inicial e considerada, nao a serie inteira") {
		// Janela: {0, 5, 0, 10} -> z1 = 7.5, gap unico = 2 -> p1 = 2
		// Os valores fora da janela (100, 100) nao devem influenciar o resultado.
		std::vector<double> ts{0.0, 5.0, 0.0, 10.0, 100.0, 100.0};
		ClassicAdHocParameters result = calculate_by_classic_ad_hoc(ts, 4);

		REQUIRE(result.z1 == Approx(7.5));
		REQUIRE(result.p1 == Approx(2.0));
	}

	SECTION("multiplos gaps distintos sao corretamente calculados") {
		// Janela: {3, 0, 0, 6, 0, 9} -> nonzero=3 (idx 0,3,5), z1=6
		// Gaps: idx0->idx3 = 3 ; idx3->idx5 = 2 -> p1 = (3+2)/2 = 2.5
		std::vector<double> ts{3.0, 0.0, 0.0, 6.0, 0.0, 9.0};
		ClassicAdHocParameters result = calculate_by_classic_ad_hoc(ts, 6);

		REQUIRE(result.z1 == Approx(6.0));
		REQUIRE(result.p1 == Approx(2.5));
	}
}
