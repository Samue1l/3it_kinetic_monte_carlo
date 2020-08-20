#pragma once
#include "libs.h"
#include "energy.h"
#include "initialization.h"

void add_rate(float rate, std::vector<double> &rates, std::vector<double> &c_rates) ;
void decrease_rate(float rate,  std::vector<double> &rates, std::vector<double> &c_rates);
void recompute_rates(unsigned char* lat, float* arr, int x, int y,int z, std::vector<double> &rates, std::vector<double> &c_rates, std::unordered_map<float, std::vector<int>> &dict, const int mx, const int my, const int mz, const float beta, const int n);
float get_choice(std::vector<double> &l, std::vector<double> &weights) ;
std::vector<int> get_transitions(float* arr, float rate, int size) ;
void element_wise_product(std::vector<double> &v1,std::vector<double> &v2,std::vector<double> &out) ;
double update_lattice(unsigned char* lat,unsigned char* copied_lat, int* epi, float* arr,std::vector<double> &rates, std::vector<double> &c_rates,std::unordered_map<float, std::vector<int>> &dict,const int mx, const int my, const int mz, const int n, const float beta);
