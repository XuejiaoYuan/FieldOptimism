#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <chrono>
#include <omp.h>
#include <cstdarg>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
using namespace std;

#include "../3rdparty/Eigen/Core"
#include "../3rdparty/Eigen/LU"
#include "../3rdparty/Eigen/src/Geometry/OrthoMethods.h"
using namespace Eigen;

#include "../Common/vector_arithmetic.cuh"
#include <device_launch_parameters.h>
#include <device_atomic_functions.h>
#include <vector_types.h>

#include <jsoncons\json.hpp>
using namespace jsoncons;



#define Epsilon		1e-6
#define VERTEXSCALE 1000000000
#define PI acos(double(-1))
#define HELIOSTAT_REFLECTIVITY 0.88 
#define RECEIVER_SLICE 0.05		// *ATTENTION*: w与l应被RECEIVER_SLICE整除
#define HELIO_SLICE 0.04		// 以该长度分割定日镜

#define DEVICE_LIST_SIZE 10		
#define DEVICE_BACKUP_SIZE 5

#define		MAX_ALL_THREADS			0x40000000000
#define		MAX_BLOCK_SINGLE_DIM	0x0ffff
#define		MAX_THREADS				1024

#define SIGMA_SUN	2.24e-3 //2.41e-3
#define SIGMA_S		1e-3


//#define DEBUG
//#define OUTPUTRES
//#define READFILE
//#define CLIPPER
//#define CALC_TIME

typedef enum {
	RectLayoutType, CrossRectLayoutType, RadialStaggerLayoutType, SpiralLayoutType
}LayoutType;


typedef enum {
	Initial, GroundMode, ReceiverMode, LayoutMode, HeliostatMode
}InputMode;

