#include "../../../Common/CommonFunc.h"
#include "../../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../../GaussLegendre/GaussLegendre.cuh"

// 
// [CUDA] 计算全镜场下接收器平面的辐射能密度分布
//
void calcRectRecvFluxDistribution(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, float* d_recv_flux, int row, int col);

__global__ void calcRectRecvDiscreteFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_recv_flux, int row, int col);

__device__ float calcRectRecvDiscreteFluxCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, int row, int col);