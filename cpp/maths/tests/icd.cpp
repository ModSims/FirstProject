#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include "maths.h"
#include <iostream>

TEST_CASE("Incomplete Cholesky Decomposition", "[functional_analysis]")
{
    Eigen::MatrixXd A(2, 2);
    Eigen::MatrixXd P(2, 2);
    Eigen::MatrixXd L(2, 2);
    Eigen::MatrixXd F(2, 2);

    A << 0.7, -0.4,
        -0.2, 0.5;
    P << 1.6129, 0.645161,
         0.645161, 2.25806;
    L << 0.83666, 0,
        -0.23905, 0.665475;
    F << 1.19523, 0,
        -4.1833, 1.50269;

    Maths::Matrix::IncompleteCholeskyDecomposition icd(A);

    REQUIRE(icd.getL().isApprox(L, 0.0001));
    REQUIRE(icd.getF().isApprox(F, 0.0001));
    REQUIRE(icd.getP().isApprox(P, 0.0001));
}