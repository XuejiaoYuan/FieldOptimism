#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../Tool/Timer/Timer.h"

class SolarScene;

class GridDDA {
public:
	void rayCastForShadow(SolarScene *solar_scene, Heliostat * helio, Vector3d dir, unordered_set<int>& rela_helio_index);	// ����shadow������ض��վ�
	void rayCastForBlock(SolarScene* solar_scene, Heliostat* helio, set<vector<int>>& rela_grid_index);	// ����block�����������
	void getBlockHelioFromGrid(SolarScene* solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Heliostat* helio);

private:
	double GridDDACore(Vector3d& dir, Heliostat* helio, Layout* layout, set<vector<int>>& relative_grid_label);
	void getHelioFromGridCore(SolarScene* solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Vector3d& dir, Heliostat* helio);
	bool checkBoundingBox(const Vector3d& Hloc, const Vector3d& Hnormal, const Vector3d& HIloc, const Vector3d& dir, double diameter);

};



