#include "../Common/CommonFunc.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../GaussLegendre/GaussLegendre.cuh"

__global__ void fluxIntegral(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n);

__device__ float calcSigma(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, int helioIndex, int recvIndex);