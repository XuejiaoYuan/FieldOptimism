#include "DiscreteRayCaster.h"


//
// [离散光线跟踪] 
//
void DiscreteRayCaster::rayCasting(SolarScene* solar_scene, string save_path, CalcMode calc_mode)
{
	// 1. 预计算相关定日镜
	RayCastHelioDeviceArgument h_args;
	auto helio_size = solar_scene->helios[0]->helio_size;
	h_args.setHelioDeviceOrigins(HELIO_SLICE, helio_size.x(), helio_size.z());
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);

	GridDDA ddaHandler;
	ddaHandler.predictRelatedHelio(solar_scene, h_args, true);
	ddaHandler.predictRelatedHelio(solar_scene, h_args, false);


	// 2. 离散光线跟踪
	Timer t;
	t.resetStart();
	vector<double> sdbk_res = rayCastCore(solar_scene->sunray_dir, h_args, save_path, calc_mode);
	t.printDuration("ray casting");

	fstream outFile(save_path + "rayCasting_" + to_string(h_args.helio_slice_length) + "_" + to_string(h_args.helio_slice_width) + ".txt", ios_base::out);
	for (int i = 0; i < sdbk_res.size(); ++i)
		outFile << setprecision(12) <<  sdbk_res[i] << endl;
	outFile.close();

	h_args.clear();
}
