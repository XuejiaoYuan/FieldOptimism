#include "EnergyCalculator.h"
#include "../Tool/Timer/Timer.h"
#include "ReceiverEnergyCalc/RecvEnergyCalc.cuh"

float* EnergyCalculator::d_energy_sum;

float EnergyCalculator::calcEnergySum(SolarScene * solar_scene, int M, int N, int m, int n)
{
	// 1. Get guass legendre calculation handler
	GaussLegendre* gl_hander = GaussLegendre::getInstance(M, N);

	// 2. Calculate receiver flux integral
	bool calcCenterMode = false;
	if (solar_scene->getModelType() == bHFLCAL)
		calcCenterMode = true;
	RecvEnergyCalc recv_energy_calc(solar_scene, gl_hander, m, n, calcCenterMode);
	recv_energy_calc.calcRecvEnergySum(getDeviceTotalEnergy());

	// 3. Return filed energy
	float h_energy_sum;
	cudaMemcpy(&h_energy_sum, d_energy_sum, sizeof(float), cudaMemcpyDeviceToHost);
	return h_energy_sum;
}

EnergyCalculator::~EnergyCalculator()
{
	if (d_energy_sum) {
		cudaFree(d_energy_sum);
		d_energy_sum = nullptr;
	}
}

float * EnergyCalculator::getDeviceTotalEnergy()
{
	if (d_energy_sum == nullptr) {
		cudaMalloc((void**)&d_energy_sum, sizeof(float));
	}
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
	float* h_energy_sum = new float;
	*h_energy_sum = 0;
	cudaMemcpy(d_energy_sum, h_energy_sum, sizeof(float), cudaMemcpyHostToDevice);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}

	delete h_energy_sum;
	return d_energy_sum;
}
