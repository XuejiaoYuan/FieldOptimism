#include "../../../Common/CommonFunc.h"
#include "../../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../../GaussLegendre/GaussLegendre.cuh"

void calcRectRecvFluxDistribution(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, float* d_recv_flux, int row, int col);

__global__ void calcRectRecvDiscreteFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_recv_flux, int row, int col);

__device__ float calcRectRecvDiscreteFluxCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, int row, int col);