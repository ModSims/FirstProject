#include "kernel.h"

using namespace Kernel;

bool Timer::isStarted() const {
    if (m_startTime == std::chrono::steady_clock::time_point()) {
        return false;
    }
    return true;
}

bool Timer::isStopped() const {
    if (m_endTime == std::chrono::steady_clock::time_point()) {
        return false;
    }
    return true;
}

void Timer::start() {
    m_startTime = std::chrono::steady_clock::now();
}

void Timer::stop() {
    m_endTime = std::chrono::steady_clock::now();
}

void Timer::update() {
    m_timeSteps++;
}

int Timer::getDurationInSeconds() {
    const long long count = std::chrono::duration_cast<std::chrono::microseconds>(m_endTime - m_startTime).count();
    return static_cast<int>(count) / 1000000;
}

int Timer::getCurrentTimeStep() const {
    return m_timeSteps;
}