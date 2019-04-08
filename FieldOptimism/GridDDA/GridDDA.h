#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../RayCastingArgument/RayCastingArgument.h"
#include "../DataStructure/DeviceList.cuh"
#include "../DataStructure/Timer.h"


class GridDDA {
public:
	~GridDDA() {
		delete[] d_shadow_related_index;
		delete[] d_block_related_index;
		d_shadow_related_index = nullptr;
		d_block_related_index = nullptr;
	}
	void predictRelatedHelio(SolarScene* solar_scene, HeliostatDeviceArgument& h_args, LayoutDeviceArgument l_args, bool shadowDir = true);
	bool checkBoundingBox(const Vector3d& Hloc, const Vector3d& Hnormal, const Vector3d& HIloc, const Vector3d& dir, double diameter, bool shadowDir= true);
	void testHandler(SolarScene* solar_scene);

	DeviceList<int>* d_shadow_related_index = nullptr;
	DeviceList<int>* d_block_related_index = nullptr;

private:
	//void setRelatedHelioIndex(vector<vector<int>>& rela_index, int*& d_rela_index, int*& d_start_index, int indexSum);
	void saveTestRes(const string& file_name, int helioNum, DeviceList<int>* d_related_index);
};



