#include "FluxIntegral.cuh"

__global__ void fluxIntegral(HeliostatDeviceArgument h_args, ReceiverDeviceArgument r_args, int m, int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats*r_args.numberOfReceivers) return;

	int helioIndex = myId / (m*n*r_args.numberOfReceivers);
	int recvIndex = (myId % (m*n*r_args.numberOfReceivers)) / (m*n);
	int row_col = (myId % (m*n*r_args.numberOfReceivers)) % (m*n);
	int curRow = row_col / n;
	int curCol = row_col%n;



}