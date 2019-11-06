#include "../../../Common/CommonFunc.h"
#include "../../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../../GaussLegendre/GaussLegendre.cuh"
#include "../RectangleReceiverFluxIntegral/RectRecvFluxIntegral.cuh"

void calcCylinderRecvEnergySum(int m, int n, int helioNum, IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl_handler, float* d_helio_energy);

__global__ void calcCylinderRecvFluxSum(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n);

__global__ void calcHelioCylinderRecvFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_helio_energy, const int m, const int n);

__device__ float calcCylinderRecvFluxIntegralCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, const int m, const int n);
