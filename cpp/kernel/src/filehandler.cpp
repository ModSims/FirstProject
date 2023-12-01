#include <iostream>
#include <fstream>
#include <iomanip>

#include "kernel.h"

using namespace Eigen;

void Kernel::saveMatrix(const char* filename, MatrixXd* matrix) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Loop through matrix and write to file
        for (int i = 0; i < matrix->rows(); i++) {
            for (int j = 0; j < matrix->cols(); j++) {
                file << std::fixed << std::setprecision(16) << matrix->coeffRef(i, j) << " ";
            }
            file << std::endl;
        }

        // Close the file
        file.close();
        std::cout << "Matrix saved to " << filename << std::endl;
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
    }
}

void Kernel::saveVector(const char* filename, VectorXd* vector) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Loop through vector and write to file
        for (int i = 0; i < vector->rows(); i++) {
            file << std::fixed << std::setprecision(16) << vector->coeffRef(i) << " ";
            file << std::endl;
        }

        // Close the file
        file.close();
        std::cout << "Vector saved to " << filename << std::endl;
    } else {
        std::cerr << "Error opening file " << filename << std::endl;
    }
}