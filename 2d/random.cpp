#include "libs.h"
std::mt19937_64 rng;
uint64_t timeSeed = 24101;
std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
// initialize a uniform distribution between 0 and 1
std::uniform_real_distribution<double> unif(0, 1);

