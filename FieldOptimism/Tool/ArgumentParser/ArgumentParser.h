#pragma once
#include <string>
#include <iostream>
#include <fstream>
using namespace std;


#include "../../DataStructure/Heliostat.h"
#include "../../DataStructure/Receiver.h"
#include "../../DataStructure/SunRay.h"


enum TimeType
{
	MonthType, DayType, HourType
};

class ArgumentParser
{
public:
	ArgumentParser() {};
	~ArgumentParser() {
		delete heliostat;
		for (auto& recv : recvs)
			delete recv;
	}
	bool parser(int argc, char** argv);
	string getInputPath() { return input_path; }
	string getOutputPath() { return output_path; }
	string getOutputFn() { return output_path + '/' + output_fn; }
	json getConfig() { return config; }
	int getNumOfHelio() { return num_of_helio; }
	ModelType getModelType() { return model_type; }
	bool getCalcSigma() { return calc_sigma; }
	LayoutType getLayoutType() { return layout_type;}
	Heliostat getHelio() { return *heliostat; }
	vector<Receiver*> getReceivers() { return recvs; }
	vector<int> getTimeParams(TimeType type);
	int getMinuteGap() { return minute_gap; }
	SunRay getSunray() { return sunray; }

private:
	string config_path;

	string input_path;
	string sunray_fn;
	string output_path;
	string output_fn;

	// config
	json config;

	// Scene
	LayoutType layout_type;
	int num_of_helio;
	ModelType model_type;
	bool calc_sigma;

	// Sunray
	SunRay sunray;

	// Heliostat
	Heliostat* heliostat;

	// Receiver
	vector<Receiver*> recvs;

	// Time
	vector<int> months;
	vector<int> days;
	vector<int> hours;
	int minute_gap;

	void getTimeParams();
};
