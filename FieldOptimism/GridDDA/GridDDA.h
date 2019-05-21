#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../DataStructure/Timer.h"


class GridDDA {
public:
	void predictRelatedHelio(SolarScene* solar_scene, RayCastHelioDeviceArgument& h_args, bool shadowDir = true);
	bool checkBoundingBox(const Vector3d& Hloc, const Vector3d& Hnormal, const Vector3d& HIloc, const Vector3d& dir, double diameter);
	void testHandler(SolarScene* solar_scene);
	void rayCastGridDDA(SolarScene *solar_scene, Heliostat * helio, Vector3d dir, unordered_set<int>& rela_helio_index);	// 计算shadow方向相关定日镜
	void rayCastGridDDA(SolarScene* solar_scene, Heliostat* helio, set<vector<int>>& rela_grid_index);	// 计算block方向相关网格
	void getBlockGrid2HelioIndex(SolarScene* solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Heliostat* helio);

private:
	void saveTestRes(const string& file_name, int helioNum, int* d_related_index, int list_size);
	double GridDDACore(Vector3d& dir, Heliostat* helio, Layout* layout, set<vector<int>>& relative_grid_label);
	void getGrid2HelioIndex(SolarScene* solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Vector3d& dir, Heliostat* helio);
};



