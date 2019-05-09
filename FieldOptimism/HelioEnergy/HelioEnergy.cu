#include "HelioEnergy.cuh"


vector<float> HelioEnergy::calcHelioEnergy(FieldUpdateMode mode)
{
	int helioNum = solar_scene->helios.size();
	float* h_helio_energy = new float[helioNum];
	for (int i = 0; i < helioNum; ++i) h_helio_energy[i] = 0;

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args.numberOfReceivers*m*n*helioNum);

	switch (mode)
	{
	case HelioUpdateMode:
		h_args.setHelioDevicePos(solar_scene->helios);
	case SunUpdateMode:
		h_args.setHelioDeviceArguments(solar_scene->helios);
		h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
		if(calcCenterMode) h_args.setHelioCenterBias(solar_scene->helios);
		break;
	default:
		break;
	}

	float* d_helio_energy = nullptr;
	cudaMalloc((void**)&d_helio_energy, sizeof(float)*helioNum);
	cudaMemcpy(d_helio_energy, h_helio_energy, sizeof(float)*helioNum, cudaMemcpyHostToDevice);
	calcHelioFluxIntegral << <nBlocks, nThreads >> > (h_args, r_args, gl_handler, d_helio_energy, m, n);
	cudaDeviceSynchronize();
	cudaMemcpy(h_helio_energy, d_helio_energy, sizeof(float)*helioNum, cudaMemcpyDeviceToHost);
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}

	vector<float> res;
	for (int i = 0; i < helioNum; ++i) res.push_back(h_helio_energy[i]);

	delete[] h_helio_energy;
	h_helio_energy = nullptr;
	cudaFree(d_helio_energy);
	return res;
}

float HelioEnergy::calcTotalEnergy(FieldUpdateMode mode)
{
	int helioNum = solar_scene->helios.size();
	float* h_total_energy = new float;
	*h_total_energy = 0;

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args.numberOfReceivers*m*n*helioNum);

	switch (mode)
	{
	case HelioUpdateMode:
		h_args.setHelioDevicePos(solar_scene->helios);
	case SunUpdateMode:
		h_args.setHelioDeviceArguments(solar_scene->helios);
		h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
		if(calcCenterMode) h_args.setHelioCenterBias(solar_scene->helios);
		break;
	default:
		break;
	}

	cudaMemcpy(d_total_energy, h_total_energy, sizeof(float), cudaMemcpyHostToDevice);
	calcFieldFluxIntegral << <nBlocks, nThreads >> > (h_args, r_args, gl_handler, d_total_energy, m, n);
	cudaDeviceSynchronize();
	cudaMemcpy(h_total_energy, d_total_energy, sizeof(float), cudaMemcpyDeviceToHost);
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}

	float res = *h_total_energy;
	delete h_total_energy;
	h_total_energy = nullptr;
	return res;
}
