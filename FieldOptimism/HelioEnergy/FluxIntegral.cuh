#include "../Common/CommonFunc.h"
#include "../RayCastingArgument/RayCastingArgument.h"

__global__ void fluxIntegral(HeliostatDeviceArgument h_args, ReceiverDeviceArgument r_args, int m, int n);