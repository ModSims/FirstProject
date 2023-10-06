#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "maths.h"

using namespace Maths;

TEST_CASE( "Radius Calculator", "[Radius]" ) {
    int iterations = 100;
    std::vector<double> spectral_radius_1(iterations);
    std::vector<double> spectral_radius_2(iterations);
    Maths::Radius radiusCalculator;

    for (int i = 0; i < iterations; ++i) {
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(100, 100);
        
        auto start_time = std::chrono::high_resolution_clock::now();
        double result_1  = radiusCalculator.calculateSpectralRadius(A);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        spectral_radius_1[i] = result_1 ;
        double time_function_1 = duration.count();

        start_time = std::chrono::high_resolution_clock::now();
        double result_2 = radiusCalculator.gelfands_spectral_approximation(A);
        end_time = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        spectral_radius_2[i] = result_2;
        double time_function_2 = duration.count();

        if (i % 100 == 0) {
            std::cout << "Iteration " << i << " of " << iterations << std::endl;
        }        
        // Output the execution times if needed
        std::cout << "Time for function 1 (microseconds): " << time_function_1 << std::endl;
        std::cout << "Time for function 2 (microseconds): " << time_function_2 << std::endl;
    }

}