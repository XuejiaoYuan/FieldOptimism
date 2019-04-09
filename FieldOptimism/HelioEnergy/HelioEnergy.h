#pragma once
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../GaussLegendre/GaussLegendre.cuh"

class HelioEnergy {
public:
	HelioEnergy(Receiver& recv, int M, int N, int m, int n) :M(M), N(N), m(m), n(n) {
		r_args.setRecvDeviceArguments(recv);
		gl_handler.initNodeWeight(M, N);
	}
	void calcHelioEnergy(SolarScene* solar_scene);

private:
	int M;	// ���ֽڵ�
	int N;	// ���ֽڵ�
	int m;	// ����ָ�
	int n;	// ����ָ�
	ReceiverDeviceArgument r_args;
	GaussLegendre gl_handler;
};