#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include "../ShadowBlockCalculator/ShadowBlockCalculator.h"
#include "../HelioEnergy/EnergyCalculator.h"
#include "../SigmaFitting/SigmaFitting.h"


enum EnergyCalcMode
{
	SceneEnergyMode, FluxDensityMode
};

class EnergyCalcPipeline
{
public:
	double handler(ArgumentParser& argumentParser, json& field_args);
	EnergyCalcPipeline();
	virtual ~EnergyCalcPipeline();

protected:
	SolarScene* solar_scene;
	ArgumentParser* argumentParser;
	void initSolarScene(json& field_args);
	virtual double handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler);
};


class FluxCalcPipeline :public EnergyCalcPipeline{
public:
	double handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler);
};

class EnergyCalcCreator {
public:
	static EnergyCalcPipeline* getPipeline(EnergyCalcMode mode) {
		switch (mode)
		{
		case SceneEnergyMode:
			return new EnergyCalcPipeline();
		case FluxDensityMode:
			return new FluxCalcPipeline();
		default:
			cerr << "Wrong mode!" << endl;
			return NULL;
		}
	}
};