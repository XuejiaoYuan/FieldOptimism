#pragma once
#include "../DataStructure/SolarScene.h"


class DiscreteRayCaster {
public:
	DiscreteRayCaster(SolarScene* _solar_scene) :solar_scene(_solar_scene) {}
	void clearArguments();
	void rayCasting();

private:
	SolarScene* solar_scene;
};