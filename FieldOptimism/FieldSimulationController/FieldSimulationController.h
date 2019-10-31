#pragma once
#include "../EnergyCalculatePipeline/EnergyCalculatePipeline.h"
#include "../DifferentialEvolution/DifferentialEvolution.h"

enum ControllerType
{
	FieldOptimismType, FluxSimulationType
};

class BaseController
{
public:
	virtual void handler(int argc, char** argv) = 0;
};


class FieldOptimismController :public BaseController{
public:
	void handler(int argc, char** argv);
};

class FluxSimulationController : public BaseController{
public:
	void handler(int argc, char** argv);
};


class ControllerCreator {
public:
	static BaseController* getController(ControllerType type) {
		switch (type)
		{
		case FieldOptimismType:
			return new FieldOptimismController();
		case FluxSimulationType:
			return new FluxSimulationController();
		default:
			return NULL;
		}
	}
};