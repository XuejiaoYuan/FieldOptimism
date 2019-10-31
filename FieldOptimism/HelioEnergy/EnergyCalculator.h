#pragma once
#include"../DataStructure/SolarScene.h"
#include "../GaussLegendre/GaussLegendre.cuh"

class EnergyCalculator {
public:
	static float calcEnergySum(SolarScene* solar_scene,json& gaussian_params);
	~EnergyCalculator();

private:
	static float* d_energy_sum;
	static float* getDeviceTotalEnergy();
};
