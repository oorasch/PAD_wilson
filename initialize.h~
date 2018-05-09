#ifndef FILE_INITIALIZE
#define FILE_INITIALIZE
#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <random>
#include "constants.h"

//const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::ranlux48 generator(seed);

/*
 *
 * name: neib_init
 * @param N1, N2
 * @return filled neighbour field
 *
 */
void neib_init (const int &N1, const int &N2, int neib[][4]);

void sig_init (std::vector<std::vector<int>>& sig_link);

void calc_M_array(const int neib [][4], const std::vector<std::vector<int>>& sig_link, const std::vector<int>& s_site, std::vector<double>& M_prime);

void bessel_init(std::vector<double>& vec);

#endif
