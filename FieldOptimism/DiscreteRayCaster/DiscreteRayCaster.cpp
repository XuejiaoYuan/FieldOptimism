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
	h_args.setHelioDeviceOrigins(150, 150);
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);

	GridDDA ddaHandler;
	ddaHandler.predictRelatedHelio(solar_scene, h_args, true);
	ddaHandler.predictRelatedHelio(solar_scene, h_args, false);


	// 2. ��ɢ���߸���
	Timer t;
	t.resetStart();
	vector<double> sdbk_res = rayCastCore(solar_scene->sunray_dir, h_args);
	t.printDuration("ray casting");

	fstream outFile("sdBkRes/rayCasting_" + to_string(h_args.helio_slice_length) + "_" + to_string(h_args.helio_slice_width) + ".txt", ios_base::out);
	for (int i = 0; i < sdbk_res.size(); ++i)
		outFile << setprecision(12) <<  sdbk_res[i] << endl;
	outFile.close();
}
