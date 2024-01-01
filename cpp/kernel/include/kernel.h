#pragma once
#include <chrono>
#include <Eigen/Dense>

namespace Kernel {
    using namespace Eigen;
    class Timer {
    public:
        Timer() {};
        Timer(double dt) : m_dt(dt) {};
        void start();
        void update();
        void stop();
        void reset();
        bool isStarted() const;
        bool isStopped() const;
        int getCurrentTimeStep() const;
        int getDurationInSeconds();        
    private:
        std::chrono::steady_clock::time_point m_startTime;
        std::chrono::steady_clock::time_point m_endTime;
        int m_timeSteps = 0;
        double m_dt = 0.0001;
    };

    void saveMatrix(const char* filename, MatrixXd* matrix);
    void saveVector(const char* filename, VectorXd* vector);
}
