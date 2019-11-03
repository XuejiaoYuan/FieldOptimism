#include "ReceiverEnergyCalculator.cuh"

float ReceiverEnergyCalculator::calcRecvEnergySum()
{
	ReceiverType recv_type = solar_scene->recvs[0]->recv_type;

	// 1. Set heliostat device arguments
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);
	h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
	if (calcCenterMode) h_args.setHelioCenterBias(solar_scene->helios);

	// 2. Calculate field total energy
	// 2.1 Set receivert arguments
	// 2.2 Calculate flux density integral
	float res = 0;
	switch (recv_type)
	{
	case RectangularRecvType:
	case PolyhedronRecvType:
		res = calcRectRecvEnergySum();
		break;
	case CylinderRecvType:
		res = calcCylinderRecvEnergySum();
		break;
	default:
		throw runtime_error("[Error RecvEnergyCalculator] Wrong receiver type!!!\n");
	}

	return res;
}

void ReceiverEnergyCalculator::setDeviceHelioEnergy() {
	int helioNum = solar_scene->helios.size();
	h_helio_energy = new float[helioNum];
	for (int i = 0; i < helioNum; ++i)
		h_helio_energy[i] = 0;
	cudaMalloc((void**)&d_helio_energy, sizeof(float)*helioNum);
	cudaMemcpy(d_helio_energy, h_helio_energy, sizeof(float)*helioNum, cudaMemcpyHostToDevice);
}

float ReceiverEnergyCalculator::calcHelioEnergySum() {
	int helioNum = solar_scene->helios.size();
	cudaMemcpy(h_helio_energy, d_helio_energy, sizeof(float)*helioNum, cudaMemcpyDeviceToHost);
	float sum = 0;
	fstream out("cuda.txt", ios_base::out);
	for (int i = 0; i < helioNum; ++i) {
		sum += h_helio_energy[i];
		out << h_helio_energy[i] << endl;
	}
	out.close();
	delete[] h_helio_energy;
	return sum;
}

float ReceiverEnergyCalculator::calcRectRecvEnergySum()
{
	int helioNum = solar_scene->helios.size();
	r_args = new ReceiverDeviceArgument();
	r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args->numberOfReceivers*m*n*helioNum);

	setDeviceHelioEnergy();
	calcHelioRectRecvFlux << <nBlocks, nThreads >> > (h_args, *r_args, *gl_handler, d_helio_energy, m, n);
	cudaDeviceSynchronize();
	float sum = calcHelioEnergySum();

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
	return sum;
}

float ReceiverEnergyCalculator::calcCylinderRecvEnergySum()
{
	int helioNum = solar_scene->helios.size();
	r_args = new CylinderReceiverDeviceArgument();
	r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, m*n*helioNum);
	
	setDeviceHelioEnergy();
	//calcCylinderRecvFluxSum << <nBlocks, nThreads >> > (h_args, *r_args, *gl_handler, d_total_energy, m, n);
	calcHelioCylinderRecvFlux << <nBlocks, nThreads >> > (h_args, *r_args, *gl_handler, d_helio_energy, m, n);
	cudaDeviceSynchronize();
	float sum = calcHelioEnergySum();

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
	return sum;
}
