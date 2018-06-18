//Created by Joseph Dunbar -- June 6, 2018

#ifndef VALUENOISE

#define VALUENOISE
#include "vectors.h"
//Functions return a value between 0-1, and the gradient at the given coordinate  
struct double5 ValueNoise4D (double x, double y, double z, double w, double frequency, int octaves, double lacunarity, double persistence);
struct double4 ValueNoise3D (double x, double y, double z, double frequency, int octaves, double lacunarity, double persistence);
struct double3 ValueNoise2D (double x, double y, double frequency, int octaves, double lacunarity, double persistence);
struct double2 ValueNoise1D (double x, double frequency, int octaves, double lacunarity, double persistence);
#endif