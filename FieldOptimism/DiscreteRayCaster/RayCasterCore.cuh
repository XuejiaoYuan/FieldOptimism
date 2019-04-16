#include "../RayCastingArgument/RayCastingArgument.h"


vector<double> rayCastCore(Vector3d sunray_dir, RayCastHelioDeviceArgument& h_args, string save_path, CalcMode calc_mode);

__global__ void rayCollisionCalc(float3 sunray_dir, RayCastHelioDeviceArgument h_args, int* d_hit_index);

__device__ bool collision(int helioIndex, float3 originPos, float3 sunray_dir, RayCastHelioDeviceArgument& h_args, bool shadowDir, int* d_hit_index);

