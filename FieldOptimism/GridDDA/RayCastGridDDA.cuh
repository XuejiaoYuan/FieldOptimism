#pragma once
#include "cuda_runtime.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../DataStructure/DeviceList.cuh"
#include "../Common/vector_arithmetic.cuh"

__global__ void rayCastGridDDACore(float3 dir, HeliostatDeviceArgument helio_device_args, LayoutDeviceArgument l_args, DeviceList<int2>* d_related_grid_list, bool shadowDir);

void rayCastGridDDA(Vector3d dir, HeliostatDeviceArgument helio_device_args, LayoutDeviceArgument l_args, DeviceList<int2>* d_related_grid_list, bool shadowDir);