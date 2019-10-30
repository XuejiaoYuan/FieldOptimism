//
// Created by Amber on 2018/4/3.
//
#pragma once

#ifndef HELIOSHADOW_SOLARSCENE_H
#define HELIOSHADOW_SOLARSCENE_H


#include "Layout.h"
#include "Receiver.h"
#include "../Tool/ArgumentParser/ArgumentParser.h"


class SolarScene {
public:
    SolarScene():model_type(HFLCAL), calc_sigma(true){}
	~SolarScene();
	bool changeHeliosNormal(const Vector3d&sunray_dir);		//change the surface normal of heiliostats
	void saveSolarScene(string scene_savepath);
	void saveHeSolarScene(string scene_savepath);
	void setModelStatus(const ModelType& model_type, const bool& calc_sigma);
	ModelType getModelType() { return model_type; }

    vector<Layout*>  layouts;
    vector<Heliostat*> helios;
	vector<Receiver*> recvs;
	Vector3d sunray_dir;

private:
	ModelType model_type;
	bool calc_sigma;
};

#endif //HELIOSHADOW_SOLARSCENE_H
