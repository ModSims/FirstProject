#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include "maths.h"

TEST_CASE("Spectral Radius - Matrix 1", "[spectral]")
{
    Eigen::MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.2, 0.5;

    SECTION("getSpectralRadiusViaEigen")
    {
        double spectral_radius = Maths::Spectral::getSpectralRadiusViaEigen(A);
        REQUIRE(spectral_radius == 0.9);
    }
    SECTION("getGelfandsSpectralApproximation")
    {
        int max_iterations = 100;
        double spectral_radius = Maths::Spectral::getGelfandsSpectralApproximation(A, max_iterations);
        REQUIRE(0.88 < spectral_radius);
        REQUIRE(spectral_radius < 0.92);
    }
}
