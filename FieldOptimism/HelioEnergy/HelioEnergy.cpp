#include "HelioEnergy.h"

void HelioEnergy::calcHelioEnergy()
{
	int helioNum = solar_scene->helios.size();
	float* h_flux_sum = new float[helioNum];
	for (int i = 0; i < helioNum; ++i) h_flux_sum[i] = 0;

	float* d_flux_sum = nullptr;
	cudaMalloc((void**)&d_flux_sum, sizeof(float)*helioNum);
	cudaMemcpy(d_flux_sum, h_flux_sum, sizeof(float)*helioNum, cudaMemcpyHostToDevice);

	int nThreads = 1024;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, solar_scene->recvs.size()*m*n*helioNum);

	IntegralHelioDeviceArgumet h_args;
	h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));

	fluxIntegral << <nBlocks, nThreads >> > (h_args, r_args, gl_handler, d_flux_sum, m, n);

	cudaMemcpy(h_flux_sum, d_flux_sum, sizeof(float)*helioNum, cudaMemcpyDeviceToHost);

	cudaFree(d_flux_sum);
	delete[] h_flux_sum;
}
