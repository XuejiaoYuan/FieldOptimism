#include "DiscreteRayCaster.h"

//
// [�����ʷ����]
//
void DiscreteRayCaster::clearArguments()
{
	
}

//
// [��ɢ���߸���] 
//
void DiscreteRayCaster::rayCasting(SolarScene* solar_scene)
{
	// 1. Ԥ������ض��վ�
	HeliostatDeviceArgument h_args;
	h_args.setHelioDeviceOrigins(100, 100);
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);

	GridDDA ddaHandler;
	ddaHandler.predictRelatedHelio(solar_scene, h_args, true);
	ddaHandler.predictRelatedHelio(solar_scene, h_args, false);


	// 2. ��ɢ���߸���
	vector<double> sdbk_res = rayCastCore(solar_scene->sunray_dir, h_args);
}
