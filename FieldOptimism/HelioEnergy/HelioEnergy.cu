#include "HelioEnergy.cuh"


void HelioEnergy::calcHelioEnergy(float sigma, FieldUpdateMode mode)
{
	int helioNum = solar_scene->helios.size();
	float* h_total_energy = new float;
	*h_total_energy = 0;

	cudaMemcpy(d_total_energy, h_total_energy, sizeof(float), cudaMemcpyHostToDevice);

	int nThreads = 1024;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args.numberOfReceivers*m*n*helioNum);

	h_args.sigma = sigma;
	switch (mode)
	{
	case HelioUpdateMode:
		h_args.setHelioDevicePos(solar_scene->helios);
	case SunUpdateMode:
		h_args.setHelioDeviceArguments(solar_scene->helios);
		h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
		break;
	default:
		break;
	}

	fluxIntegral << <nBlocks, nThreads >> > (h_args, r_args, gl_handler, d_total_energy, m, n);
	cudaDeviceSynchronize();

	cudaMemcpy(h_total_energy, d_total_energy, sizeof(float), cudaMemcpyDeviceToHost);

	cout << *h_total_energy << endl;
	delete h_total_energy;
}
