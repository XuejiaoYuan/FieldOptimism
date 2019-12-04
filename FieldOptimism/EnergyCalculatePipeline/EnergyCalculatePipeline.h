#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include "../ShadowBlockCalculator/ShadowBlockCalculator.h"
#include "../SigmaFitting/SigmaFitting.h"
#include "ReceiverEnergyCalculator/ReceiverEnergyCalculator.h"

enum EnergyCalculatePipelineMode
{
	SceneEnergyMode, HelioFluxDensityMode, FieldFluxDensityMode
};

class EnergyCalculatePipeline
{
public:
	EnergyCalculatePipeline(ArgumentParser& argumentParser) : argumentParser(&argumentParser) {
		solar_scene = new SolarScene();
	};
	double handler(json& field_args);
	double handler();
	void saveSolarScene();
	virtual ~EnergyCalculatePipeline();

protected:
	SolarScene* solar_scene;
	ArgumentParser* argumentParser;
	void initSolarScene(json& field_args);
	double traverseTimeHandler();
	void oneTimeHandlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler);
	virtual void handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
	virtual double removeHeliostat();
};


class HelioFluxCalculatePipeline :public EnergyCalculatePipeline {
public:
	HelioFluxCalculatePipeline(ArgumentParser& argumentParser) :EnergyCalculatePipeline(argumentParser) {}
	void handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
	double removeHeliostat() { return 0; }
};

class FieldFluxCalculatePipeline : public EnergyCalculatePipeline {
public:
	FieldFluxCalculatePipeline(ArgumentParser& argumentParser) : EnergyCalculatePipeline(argumentParser) {}
	void handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
	double removeHeliostat() { return 0; }
};

class EnergyCalculateCreator {
public:
	static EnergyCalculatePipeline* getPipeline(EnergyCalculatePipelineMode mode, ArgumentParser& argumentParser) {
		switch (mode)
		{
		case SceneEnergyMode:
			return new EnergyCalculatePipeline(argumentParser);
		case HelioFluxDensityMode:
			return new HelioFluxCalculatePipeline(argumentParser);
		case FieldFluxDensityMode:
			return new FieldFluxCalculatePipeline(argumentParser);
		default:
			cerr << "Wrong mode!" << endl;
			return NULL;
		}
	}
};