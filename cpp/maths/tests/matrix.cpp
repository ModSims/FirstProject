#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include "maths.h"

TEST_CASE("Functional analysis of matrix - Matrix 1", "[functional_analysis]")
{
    Eigen::MatrixXd A(2, 2);
    A << 0.7, -0.4,
        -0.2, 0.5;

    SECTION("getConditionNumber")
    {
        double condition_number = Maths::Matrix::getConditionNumber(A);
        REQUIRE(condition_number <= 3.666);
    }

    SECTION("isSymmetric")
    {
        bool is_symmetric = Maths::Matrix::isSymmetric(A);
        REQUIRE(is_symmetric == false);
    }

    SECTION("isDiagonallyDominant")
    {
        bool is_diagonally_dominant = Maths::Matrix::isDiagonallyDominant(A);
        REQUIRE(is_diagonally_dominant == true);
    }

    SECTION("isSPD")
    {
        bool is_spd = Maths::Matrix::isSPD(A);
        REQUIRE(is_spd == false);
    }
}

TEST_CASE("Functional analysis of matrix - Matrix 2", "[functional_analysis]")
{
    Eigen::MatrixXd A(2, 2);
    A << 0.7, -0.2,
        -0.2, 0.5;

    SECTION("getConditionNumber")
    {
        double condition_number = Maths::Matrix::getConditionNumber(A);
        REQUIRE(condition_number <= 3.666);
    }

    SECTION("isSymmetric")
    {
        bool is_symmetric = Maths::Matrix::isSymmetric(A);
        REQUIRE(is_symmetric == true);
    }

    SECTION("isDiagonallyDominant")
    {
        bool is_diagonally_dominant = Maths::Matrix::isDiagonallyDominant(A);
        REQUIRE(is_diagonally_dominant == true);
    }

    SECTION("isSPD")
    {
        bool is_spd = Maths::Matrix::isSPD(A);
        REQUIRE(is_spd == true);
    }
}
