#pragma once
#include <chrono>

namespace Kernel {
    class Timer {
    public:
        Timer(double dt) : m_dt(dt) {};
        void start();
        void update();
        void stop();
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
}
