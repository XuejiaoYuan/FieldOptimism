#pragma once
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../GaussLegendre/GaussLegendre.cuh"
#include "FluxIntegral.cuh"

class HelioEnergy {
public:
	HelioEnergy(SolarScene* _solar_scene, int M, int N, int m, int n) :M(M), N(N), m(m), n(n), solar_scene(_solar_scene) {
		r_args.setRecvDeviceArguments(*(solar_scene->recvs[0]));
		gl_handler.initNodeWeight(M, N);
	}
	void calcHelioEnergy();

private:
	int M;	// 积分节点
	int N;	// 积分节点
	int m;	// 区域分割
	int n;	// 区域分割
	ReceiverDeviceArgument r_args;
	GaussLegendre gl_handler;
	SolarScene* solar_scene;
};