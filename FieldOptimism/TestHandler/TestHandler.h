#pragma once
#include "../sdbkCalculator/sdbkCalculator.h"
#include "../DataStructure/SunRay.h"
#include "../DiscreteRayCaster/DiscreteRayCaster.h"
#include "../DataStructure/Timer.h"
#include "../HelioEnergy/HelioEnergy.cuh"
#include "../SigmaFitting/SigmaFitting.h"
#include "../CalcPolygonCenter/CalcPolygonCenter.h"
#include <direct.h>

class TestHandler {
public:
	TestHandler(int m_inter) :minute_inter(m_inter) {
		sample_day = { make_int2(3, 21), make_int2(6, 22), make_int2(9, 23), make_int2(12, 21) };
		sample_hour = { 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	}
	void initSolarScene(string scene_file_path, string sunray_filepath);
	void createSolarScene(string scene_filepath);
	void testSunRay();
	void testHelioPredict();
	void testRayCasting();
	void testPolygonClipping();
	void testTotalEnergyCalc(int M, int N, int m, int n);
	void testHelioEnergyCalc(int M, int N, int m, int n);
	void testOneTimeHelioEnergy(int M, int N, int m, int n);
	void testFitSigme(int M, int N, int m, int n);
	void testCalcPolygonCenter();

	SolarScene *solar_scene;
	SunRay sunray;
	vector<int2> sample_day;
	vector<int> sample_hour;
	int minute_inter;
};