#pragma once
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <random>
#include <iterator>
#include <chrono>
#include <unordered_map>


extern std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
//uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

extern uint64_t timeSeed;
extern std::seed_seq ss;
// initialize a uniform distribution between 0 and 1
extern std::uniform_real_distribution<double> unif;
