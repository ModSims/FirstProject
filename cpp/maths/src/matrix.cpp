#include <iostream>
#include <Eigen/Dense>
#include "maths.h"

double Maths::Matrix::getConditionNumber(Eigen::MatrixXd& A) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    return cond;
}

bool Maths::Matrix::isSymmetric(Eigen::MatrixXd& A) {
    return A.transpose() == A;
}

bool Maths::Matrix::isDiagonallyDominant(Eigen::MatrixXd& A) {
    for (int i = 0; i < A.rows(); ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < A.cols(); ++j) {
            if (i != j) {
                row_sum += abs(A(i, j));
            }
        }

        if (abs(A(i, i)) < row_sum) {
            return false;
        }
    }

    return true;
}

bool Maths::Matrix::isSPD(Eigen::MatrixXd& A) {
    return Maths::Matrix::isSymmetric(A) && Maths::Matrix::isDiagonallyDominant(A);
}

Maths::Matrix::IncompleteCholeskyDecomposition::IncompleteCholeskyDecomposition(Eigen::MatrixXd& A) {
    double eps = 1e-10;
    this->L = Eigen::MatrixXd::Zero(A.rows(), A.cols());
    this->F = Eigen::MatrixXd::Zero(A.rows(), A.cols());

    for (int k = 0; k < A.rows(); ++k) {
        double row_sum = 0.0;
        for (int j = 0; j < k; ++j) {
            row_sum += pow(this->L(k, j), 2);
        }
        
        double lkk = sqrt(A(k, k) - row_sum);
        if (abs(lkk) > eps) {
            this->L(k, k) = lkk;
            this->F(k, k) = 1/lkk;
        }

        for (int i = k + 1; i < A.rows(); ++i) {
            double col_sum = 0.0;
            for (int j = 0; j < k; ++j) {
                col_sum += pow(this->L(i, j), 2);
            }

            double lik = (1/lkk)*(A(i, k) - col_sum);

            if (abs(lik) > eps) {
                this->L(i, k) = lik;
                this->F(i, k) = 1/lik;
            }
        }
    }
}

Eigen::MatrixXd Maths::Matrix::IncompleteCholeskyDecomposition::getL() {
    return this->L;
}

Eigen::MatrixXd Maths::Matrix::IncompleteCholeskyDecomposition::getF() {
    return this->F;
}


Eigen::MatrixXd Maths::Matrix::IncompleteCholeskyDecomposition::getP() {
    this->P = L.transpose().inverse() * L.inverse();
    return this->P;
}
