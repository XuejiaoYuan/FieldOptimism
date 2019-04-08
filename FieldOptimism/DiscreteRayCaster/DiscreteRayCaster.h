#pragma once
#include "../DataStructure/SolarScene.h"
#include "../GridDDA/GridDDA.h"
#include "RayCasterCore.cuh"

class DiscreteRayCaster {
public:
	void clearArguments();
	void rayCasting(SolarScene* solar_sene);

private:
};