#include "EnergyCalcPipeline.h"

double EnergyCalcPipeline::handler(ArgumentParser& _argumentParser, json& field_args) {
	// 1. 初始化镜场
	argumentParser = &_argumentParser;
	initSolarScene(field_args);

	// 2. 设置采样时间
	vector<int> months = argumentParser->getTimeParams(MonthType);
	vector<int> days = argumentParser->getTimeParams(DayType);
	vector<int> hours = argumentParser->getTimeParams(HourType);
	int minute_gap = argumentParser->getMinuteGap();

	// 3. 计算采样时刻镜场能量
	vector<int> time_params;
	SdBkCalc* sdbk_handler = new SdBkCalc(solar_scene);
	double total_sum = 0;
	for (auto& m : months)
		for (auto& d : days)
			for (auto& h : hours)
				for (int t = 0; t < 60; t += minute_gap) {
					time_params = { m, d, h, t };
					total_sum += handlerCore(time_params, argumentParser->getSunray(), sdbk_handler);
				}
	return total_sum;
}


void EnergyCalcPipeline::initSolarScene(json& field_args) {
	// 1. 设置唯一Receiver
	solar_scene->recvs = argumentParser->getReceivers();

	// 2. 设置Layout
	Layout* layout = LayoutCreator::getLayout(argumentParser->getLayoutType());
	solar_scene->layouts.push_back(layout);

	// 3. 设置Heliostat在Layout的布局
	solar_scene->layouts[0]->createHelioAndLayout(*argumentParser, field_args, solar_scene->helios);

	// 4. 设置镜面仿真模型
	solar_scene->setModelStatus(argumentParser->getModelType(), argumentParser->getCalcSigma());
}

double EnergyCalcPipeline::handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	// 1. 调整镜场定日镜朝向
	Vector3d sunray_dir = sunray.changeSunRay(time_param);
	double DNI = sunray.calcDNI(time_param);
	solar_scene->changeHeliosNormal(sunray_dir);

	// 2. 计算阴影遮挡率
	sdbk_handler->calcSceneShadowBlock();

	// 3. 计算镜场能量
	double gpu_res = EnergyCalculator::calcEnergySum(solar_scene, 8, 8, 2, 2);

	return gpu_res;
}


EnergyCalcPipeline::EnergyCalcPipeline()
{
	solar_scene = new SolarScene();
}

EnergyCalcPipeline::~EnergyCalcPipeline()
{
	delete solar_scene;
}

double FluxCalcPipeline::handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	json flux = argumentParser->getConfig();
	vector<int> rt_helio_index(flux["RayTracingTestHelioIndex"].as<vector<int>>());
	vector<int> test_helio_index(flux["TestHelioIndex"].as<vector<int>>());
	vector<Heliostat*> rt_helio;
	for (auto& h_idx : rt_helio_index)
		rt_helio.push_back(solar_scene->helios[h_idx]);
	
	SigmaFitting sg_handler;
	double DNI = sunray.calcDNI(time_param);
	string time_str = "M" + to_string(time_param[0]) + "D" + to_string(time_param[1])
		+ "H" + to_string(time_param[2]) + "m" + to_string(time_param[3]) + "/";
	sg_handler.fromFluxPeak2Sigma(argumentParser->getInputPath() + time_str + "max_flux.txt", rt_helio, solar_scene->recvs[0], DNI);
	
	sdbk_handler->setOutputPath(argumentParser->getOutputPath() + time_str);
	sdbk_handler->calcSceneFluxDistribution(test_helio_index, DNI, flux.get_with_default("GaussianParams").as<json>());
	
	return 0;
}
