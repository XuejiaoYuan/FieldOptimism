#include "EnergyCalcPipeline.h"
#include "../DataStructure/Timer.h"

float* EnergyCalcPipeline::d_total_energy;

float EnergyCalcPipeline::calcFieldEnergySum(SolarScene * solar_scene, int M, int N, int m, int n, bool calcCenterMode)
{
	// 1. Get guass legendre calculation handler
	GaussLegendre* gl_hander = GaussLegendre::getInstance(M, N);

	// 2. Calculate receiver flux integral
	RecvEnergyCalc recv_energy_calc(solar_scene, gl_hander, m, n, calcCenterMode);
	recv_energy_calc.calcRecvEnergySum(getDeviceTotalEnergy());

	// 3. Return filed energy
	float h_total_energy;
	cudaMemcpy(&h_total_energy, d_total_energy, sizeof(float), cudaMemcpyDeviceToHost);
	return h_total_energy;
}

EnergyCalcPipeline::~EnergyCalcPipeline()
{
	if (d_total_energy) {
		cudaFree(d_total_energy);
		d_total_energy = nullptr;
	}
}

float * EnergyCalcPipeline::getDeviceTotalEnergy()
{
	if (d_total_energy == nullptr) {
		cudaMalloc((void**)&d_total_energy, sizeof(float));
	}
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}
	float* h_total_energy = new float;
	*h_total_energy = 0;
	cudaMemcpy(d_total_energy, h_total_energy, sizeof(float), cudaMemcpyHostToDevice);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
	}

	delete h_total_energy;
	return d_total_energy;
}
