#pragma once

#include <iostream>
#include <chrono>

class ScopeTimer {
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;
public:
    ScopeTimer() : start(Clock::now()), msg("Duration") {}

    explicit ScopeTimer(const std::string& msg) : start(Clock::now()), msg(msg) {
        std::cout << msg << "..." << std::endl;
    }

    ~ScopeTimer() {
        TimePoint end = Clock::now();
        int64_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << msg << ": " << duration << " ms" << std::endl;
    }

private:
    const TimePoint start;
    const std::string msg;
};
