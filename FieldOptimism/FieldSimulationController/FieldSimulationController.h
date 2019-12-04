#pragma once
#include "../EnergyCalculatePipeline/EnergyCalculatePipeline.h"
#include "../DifferentialEvolution/DifferentialEvolution.h"

enum ControllerMode
{
	FieldOptimismMode, HelioFluxSimulationMode, FieldFluxSimulationMode, EnergyCalculateMode
};

class BaseController
{
public:
	virtual void handler(int argc, char** argv) = 0;
};

class EnergyCalculateController: public BaseController {
public:
	void handler(int argc, char** argv);
};

class FieldOptimismController :public BaseController{
public:
	void handler(int argc, char** argv);
};

class HelioFluxSimulationController : public BaseController{
public:
	void handler(int argc, char** argv);
};

class FieldFluxSimulationController : public BaseController {
public:
	void handler(int argc, char** argv);
};


class ControllerCreator {
public:
	static BaseController* getController(ControllerMode type) {
		switch (type)
		{
		case FieldOptimismMode:
			return new FieldOptimismController();
		case HelioFluxSimulationMode:
			return new HelioFluxSimulationController();
		case FieldFluxSimulationMode:
			return new FieldFluxSimulationController();
		case EnergyCalculateMode:
			return new EnergyCalculateController();
		default:
			return NULL;
		}
	}
};