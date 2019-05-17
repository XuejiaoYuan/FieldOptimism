#pragma once
#include"../DataStructure/SolarScene.h"
#include "../GaussLegendre/GaussLegendre.cuh"

class EnergyCalcPipeline {
public:
	static float calcFieldEnergySum(SolarScene* solar_scene, int M, int N, int m, int n, bool calcCenterMode = false);
	~EnergyCalcPipeline();

private:
	static float* d_total_energy;
	static float* getDeviceTotalEnergy();
};
