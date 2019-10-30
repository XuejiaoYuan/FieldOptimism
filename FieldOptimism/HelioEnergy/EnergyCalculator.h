#pragma once
#include"../DataStructure/SolarScene.h"
#include "../GaussLegendre/GaussLegendre.cuh"

class EnergyCalculator {
public:
	static float calcEnergySum(SolarScene* solar_scene, int M, int N, int m, int n);
	~EnergyCalculator();

private:
	static float* d_energy_sum;
	static float* getDeviceTotalEnergy();
};
