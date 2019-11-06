#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include "../ShadowBlockCalculator/ShadowBlockCalculator.h"
#include "../SigmaFitting/SigmaFitting.h"
#include "ReceiverEnergyCalculator/ReceiverEnergyCalculator.h"

enum EnergyCalcMode
{
	SceneEnergyMode, FluxDensityMode
};

class EnergyCalculatePipeline
{
public:
	EnergyCalculatePipeline(ArgumentParser& argumentParser) : argumentParser(&argumentParser) {
		solar_scene = new SolarScene();
	};
	double handler(json& field_args);
	virtual ~EnergyCalculatePipeline();

protected:
	SolarScene* solar_scene;
	ArgumentParser* argumentParser;
	void initSolarScene(json& field_args);
	void handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler);
	virtual void handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
	virtual double removeHeliostat();
};


class FluxCalculatePipeline :public EnergyCalculatePipeline {
public:
	FluxCalculatePipeline(ArgumentParser& argumentParser) :EnergyCalculatePipeline(argumentParser) {}
	void handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
	double removeHeliostat() { return 0; }
};

class EnergyCalculateCreator {
public:
	static EnergyCalculatePipeline* getPipeline(EnergyCalcMode mode, ArgumentParser& argumentParser) {
		switch (mode)
		{
		case SceneEnergyMode:
			return new EnergyCalculatePipeline(argumentParser);
		case FluxDensityMode:
			return new FluxCalculatePipeline(argumentParser);
		default:
			cerr << "Wrong mode!" << endl;
			return NULL;
		}
	}
};