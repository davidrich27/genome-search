#ifndef _CLOCK_HPP
#define _CLOCK_HPP

#include <iostream>
#include <iomanip>
#include <chrono>

class Clock {

	protected:

	std::chrono::steady_clock::time_point start_time;

	public:

	Clock() {
		tick();
	}

	void tick() {
		start_time = std::chrono::steady_clock::now();
	}

	float tock() {
		std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
		return float(std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()) / 1e6f;
	}

	// Print elapsed time with newline
	void ptock() {
		float elapsed = tock();
		std::cout << "Took " << elapsed << " seconds" << std::endl;
	}
	
};

#endif