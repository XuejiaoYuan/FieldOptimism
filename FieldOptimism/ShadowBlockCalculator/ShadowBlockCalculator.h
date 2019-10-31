#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../GaussLegendre/GaussLegendre.h"
#include "../GridDDA/GridDDA.h"
#include "../PolygonCenterCalculator/PolygonCenterCalculator.h"


class SdBkCalc
{
public:
	SdBkCalc(SolarScene* _solar_scene) :calcCenterMode(false), output_path("") {
		this->solar_scene = _solar_scene;
		if (solar_scene->getModelType() == bHFLCAL) this->calcCenterMode = true;
	}
	void calcSceneShadowBlock();
	void calcSceneFluxDistribution(vector<int>& test_helio_index, const double DNI, json& config);
	void setOutputPath(string& _output_path) { output_path = _output_path; }
	
protected:
	GaussLegendreCPU* gl;
	map<int, set<vector<int>>> rela_block_grid_index;
	SolarScene* solar_scene;
	bool calcCenterMode;
	string output_path;
	double calcHelioShadowBlock(int helio_index);
	double helioClipper(Heliostat * helio, const vector<Vector3d>& dir, const vector<unordered_set<int>>& estimate_grids);
	float calcHelio2RecvEnergy(vector<Vector3d>& recv_v, Vector3d& recv_n, const int rows, const int cols, Heliostat* helio, Vector3d& fc_center, double DNI, double cos_phi);
};


class SdBkCalcTest:public SdBkCalc{
public:
	SdBkCalcTest(SolarScene* solar_scene):SdBkCalc(solar_scene) {}
	void readRayTracingRes();
	void readRayTracingRes(int index);
	void rayTracingSdBk();
	void normalSdBk();
	void boundingSphereSdBk();
	void neighRowSdBk();
	void improvedNeighSdBk();
	void use3dddaSdBk();
	void totalHeliosTest(const string& _save_path);
	void singleHelioTest(const int _helio_index);
	void setDir(Vector3d _sd_dir, Vector3d _bk_dir) { sd_dir = _sd_dir; bk_dir = _bk_dir; }
	void setTestIndex(const int _helio_index) { helio_index = _helio_index; }

private:
	void readRayTracingCore(string file_name, vector<unordered_set<int>>& sd_index_set, vector<unordered_set<int>>& bk_index_set);
	bool calcIntersect(Vector3d& ori_v, Vector3d& dir, unordered_set<int>& index_set);
	bool checkHelioDis(Vector3d& dir, Heliostat* helio, Vector3d& Hloc, Vector3d& HIloc);
	bool checkBoundingBox(Vector3d& Hloc, Vector3d& HIloc, Vector3d& dir, double diameter);
	void getStartEndIndex(Heliostat* helio, int& start, int& end);
	bool checkEffectRegion(Vector3d dir, Vector3d& HIloc, Vector3d& Hloc, double diameter);
	bool checkCenterDist(Heliostat* H, Heliostat* HI, Vector3d& dir);
	void checkEstimateHelio(Vector3d& dir, unordered_set<int>& helio_set, int& ac, unordered_set<int>& gt_helio_set, bool shadowDir);
	void saveTestRes(const string& file_name, const int sd_ac, const int bk_ac, const unordered_set<int>& sd_set, const unordered_set<int>& bk_set);
	Vector3d sd_dir, bk_dir;
	int helio_index;
	string save_path;
	unordered_set<int> gt_sd_helio_index;		// 通过ray tracing计算得到，实际造成阴影的其他定日镜
	unordered_set<int> gt_bk_helio_index;		// 通过ray tracing计算得到，实际造成遮挡的其他定日镜

	vector<unordered_set<int>> v_gt_sd_helio_index;
	vector<unordered_set<int>> v_gt_bk_helio_index;
};