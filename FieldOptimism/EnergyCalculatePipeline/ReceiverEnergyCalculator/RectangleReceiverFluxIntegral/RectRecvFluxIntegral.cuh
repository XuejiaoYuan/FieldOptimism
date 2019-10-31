#include "../../../Common/CommonFunc.h"
#include "../../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../../GaussLegendre/GaussLegendre.cuh"

__global__ void calcRectRecvFluxSum(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n);

__global__ void calcHelioRectRecvFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_helio_energy, const int m, const int n);

__device__ float calcRectRecvFluxIntegralCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, const int m, const int n);

__device__ float calcRecvFluxIntegralCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, int helioIndex, int recvIndex, int i, int j, int m, int n);