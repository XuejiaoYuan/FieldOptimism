#include "RecvEnergyCalc.cuh"

void RecvEnergyCalc::calcRecvEnergySum(float * d_total_energy)
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
	switch (recv_type)
	{
	case RectangularRecvType:
	case PolyhedronRecvType:
		calcRectRecvEnergySum(d_total_energy);
		break;
	case CylinderRecvType:
		calcCylinderRecvEnergySum(d_total_energy);
		break;
	case CircularTruncatedConeRecvType:
		throw runtime_error("[Error RecvEnergyCalc] Receiver type not implement!!!\n");
	default:
		throw runtime_error("[Error RecvEnergyCalc] Wrong receiver type!!!\n");
	}

	// 3. Clear 
	h_args.clear();
	r_args->clear();
	delete r_args;
}

void RecvEnergyCalc::calcRectRecvEnergySum(float * d_total_energy)
{
	int helioNum = solar_scene->helios.size();
	r_args = new ReceiverDeviceArgument();
	r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args->numberOfReceivers*m*n*helioNum);

	calcRectRecvFluxSum << <nBlocks, nThreads >> > (h_args, *r_args, *gl_handler, d_total_energy, m, n);
	cudaDeviceSynchronize();
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
}

void RecvEnergyCalc::calcCylinderRecvEnergySum(float * d_total_energy)
{
	int helioNum = solar_scene->helios.size();
	r_args = new CylinderReceiverDeviceArgument();
	r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);

	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, m*n*helioNum);

	calcCylinderRecvFluxSum << <nBlocks, nThreads >> > (h_args, *r_args, *gl_handler, d_total_energy, m, n);
	cudaDeviceSynchronize();
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
}
