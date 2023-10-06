#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "math.h"

using namespace Math;

TEST_CASE( "Radian Calculator", "[Radians]" ) {
    int iterations = 10000;
    std::vector<double> spectral_radius_1(iterations);
    std::vector<double> spectral_radius_2(iterations);
    Math::Radians radiansCalculator;

    for (int i = 0; i < iterations; ++i) {
        // make a random matrix
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(100, 100);
        // start timer and time function 1
        // Start timer and time function 1
        auto start_time = std::chrono::high_resolution_clock::now();
        double result_1  = radiansCalculator.calculateSpectralRadius(A);
        auto end_time = std::chrono::high_resolution_clock::now();

        // Calculate elapsed time
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Save the result and time
        spectral_radius_1[i] = result_1 ;
        double time_function_1 = duration.count();

        // Start timer and time function 2
        start_time = std::chrono::high_resolution_clock::now();
        double result_2 = radiansCalculator.gelfands_spectral_approximation(A);
        end_time = std::chrono::high_resolution_clock::now();

        // Calculate elapsed time
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Save the result and time
        spectral_radius_2[i] = result_2;
        double time_function_2 = duration.count();

        if (i % 100 == 0) {
            std::cout << "Iteration " << i << " of " << iterations << std::endl;
        }

        REQUIRE(spectral_radius_1[i] == spectral_radius_2[i]);
        
        // Output the execution times if needed
        std::cout << "Time for function 1 (microseconds): " << time_function_1 << std::endl;
        std::cout << "Time for function 2 (microseconds): " << time_function_2 << std::endl;
    }

}