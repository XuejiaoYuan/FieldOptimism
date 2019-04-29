#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../DataStructure/DeviceList.cuh"
#include "../DataStructure/Timer.h"


class GridDDA {
public:
	void predictRelatedHelio(SolarScene* solar_scene, RayCastHelioDeviceArgument& h_args, bool shadowDir = true);
	bool checkBoundingBox(const Vector3d& Hloc, const Vector3d& Hnormal, const Vector3d& HIloc, const Vector3d& dir, double diameter);
	void testHandler(SolarScene* solar_scene);
	void rayCastGridDDA(SolarScene *solar_scene, Heliostat * helio, Vector3d dir, unordered_set<int>& rela_helio_index, bool shadowDir);


private:
	void saveTestRes(const string& file_name, int helioNum, int* d_related_index, int list_size);
};



