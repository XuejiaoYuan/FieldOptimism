#include "CylinderRecvFluxIntegral.cuh"

void calcCylinderRecvEnergySum(int m, int n, int helioNum, IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl_handler, float* d_helio_energy)
{
	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, helioNum*m*n);

	calcHelioCylinderRecvFlux << <nBlocks, nThreads >> > (h_args, r_args, gl_handler, d_helio_energy, m, n);
	cudaDeviceSynchronize();

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
}

__global__ void calcCylinderRecvFluxSum(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n) {
	float res = calcCylinderRecvFluxIntegralCore(h_args, r_args, gl, m, n);
	if (res < Epsilon) return;
	atomicAdd(d_total_energy, res);
}

__global__ void calcHelioCylinderRecvFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_helio_energy, const int m, const int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats) return;

	float res = calcCylinderRecvFluxIntegralCore(h_args, r_args, gl, m, n);

	int helioIndex = myId / (m*n*r_args.numberOfReceivers);
	atomicAdd(d_helio_energy + helioIndex, res);
}


__device__ float calcCylinderRecvFluxIntegralCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, const int m, const int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats) return -1;

	int helioIndex = myId / (m*n);
	int recvIndex = helioIndex;
	int row_col = myId % (m*n);
	int i = row_col / n;
	int j = row_col % n;

	return calcRecvFluxIntegralCore(h_args, r_args, gl, helioIndex, recvIndex, i, j, m, n);
}