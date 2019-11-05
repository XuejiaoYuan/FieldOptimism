#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include "../ShadowBlockCalculator/ShadowBlockCalculator.h"
#include "../SigmaFitting/SigmaFitting.h"
#include "ReceiverEnergyCalculator/ReceiverEnergyCalculator.cuh"

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
	double handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler);
	virtual double handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
};


class FluxCalculatePipeline :public EnergyCalculatePipeline {
public:
	FluxCalculatePipeline(ArgumentParser& argumentParser) :EnergyCalculatePipeline(argumentParser) {}
	double handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
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