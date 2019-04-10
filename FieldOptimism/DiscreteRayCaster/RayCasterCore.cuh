#include "../RayCastingArgument/RayCastingArgument.h"


vector<double> rayCastCore(Vector3d sunray_dir, RayCastHelioDeviceArgument& h_args);

__global__ void rayCollisionCalc(float3 sunray_dir, RayCastHelioDeviceArgument h_args);

__device__ bool collision(int helioIndex, float3 originPos, float3 sunray_dir, RayCastHelioDeviceArgument& h_args, bool shadowDir);