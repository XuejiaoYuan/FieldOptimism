#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include "../ShadowBlockCalculator/ShadowBlockCalculator.h"
#include "../EnergyCalculatePipeline/EnergyCalculator.h"
#include "../SigmaFitting/SigmaFitting.h"


enum EnergyCalcMode
{
	SceneEnergyMode, FluxDensityMode
};

class EnergyCalculatePipeline
{
public:
	double handler(ArgumentParser& argumentParser, json& field_args);
	EnergyCalculatePipeline();
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
	double handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbkHandler);
};

class EnergyCalculateCreator {
public:
	static EnergyCalculatePipeline* getPipeline(EnergyCalcMode mode) {
		switch (mode)
		{
		case SceneEnergyMode:
			return new EnergyCalculatePipeline();
		case FluxDensityMode:
			return new FluxCalculatePipeline();
		default:
			cerr << "Wrong mode!" << endl;
			return NULL;
		}
	}
};