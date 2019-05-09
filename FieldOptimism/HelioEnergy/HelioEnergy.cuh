#pragma once
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../GaussLegendre/GaussLegendre.cuh"
#include "FluxIntegral.cuh"


class HelioEnergy {
public:
	HelioEnergy(SolarScene* solar_scene, int M, int N, int m, int n, bool calcCenterMode = false) 
		:M(M), N(N), m(m), n(n), solar_scene(solar_scene), d_total_energy(nullptr), calcCenterMode(calcCenterMode) {
		r_args.setRecvDeviceArguments(*(solar_scene->recvs[0]));
		gl_handler.initNodeWeight(M, N);
		cudaMalloc((void**)&d_total_energy, sizeof(float));
		h_args.setHelioDevicePos(solar_scene->helios);
	}
	vector<float> calcHelioEnergy(FieldUpdateMode mode);
	float calcTotalEnergy(FieldUpdateMode mode);
	void setHFCALMode(bool calcCenterMode) { this->calcCenterMode = calcCenterMode; }
	~HelioEnergy() {
		h_args.clear();
		gl_handler.clear();
		cudaFree(d_total_energy);
		d_total_energy = nullptr;
	}
	
	bool calcCenterMode;

private:
	int M;	// 积分节点
	int N;	// 积分节点
	int m;	// 区域分割
	int n;	// 区域分割
	IntegralHelioDeviceArgumet h_args;
	ReceiverDeviceArgument r_args;
	GaussLegendre gl_handler;
	SolarScene* solar_scene;
	float* d_total_energy;
};