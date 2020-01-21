#include <chrono>
#include <iostream>

static std::chrono::time_point<std::chrono::steady_clock> start____;

inline void MeasureTimeFrom() {
  start____ = std::chrono::steady_clock::now();
}

inline int MeasureNanoSecondTime() {
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start____).count();
} 

inline int MeasureMicroSecondTime() {
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - start____).count();
} 

void DisplayNanoSecondTime(std::string prefix="") {
  std::cout << prefix << " " << MeasureNanoSecondTime() << " ns\n";
}

void DisplayMicroSecondTime(std::string prefix="") {
  std::cout << prefix << " " << MeasureMicroSecondTime() << " ns\n";
}


