#include "DiscreteRayCaster.h"

//
// [清除历史参数]
//
void DiscreteRayCaster::clearArguments()
{
	
}

//
// [离散光线跟踪] 
//
void DiscreteRayCaster::rayCasting(SolarScene* solar_scene)
{
	// 1. 预计算相关定日镜
	HeliostatDeviceArgument h_args;
	h_args.setHelioDeviceOrigins(100, 100);
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);

	GridDDA ddaHandler;
	ddaHandler.predictRelatedHelio(solar_scene, h_args, true);
	ddaHandler.predictRelatedHelio(solar_scene, h_args, false);


	// 2. 离散光线跟踪
	vector<double> sdbk_res = rayCastCore(solar_scene->sunray_dir, h_args);
}
