#include "EnergyCalculatePipeline.h"

double EnergyCalculatePipeline::handler(json& field_args) {
	// 1. 初始化镜场
	//cout << "2. Initialize solar scene" << endl;
	initSolarScene(field_args);

	// 2. 设置采样时间
	vector<int> months = argumentParser->getTimeParams(MonthType);
	vector<int> days = argumentParser->getTimeParams(DayType);
	vector<int> hours = argumentParser->getTimeParams(HourType);
	int minute_gap = argumentParser->getMinuteGap();

	// 3. 计算采样时刻镜场能量
	//cout << "3. Start simulate flux density distribution" << endl;
	vector<int> time_params;
	SdBkCalc sdbk_handler(solar_scene);
	double total_sum = 0;
	for (auto& m : months)
		for (auto& d : days)
			for (auto& h : hours)
				for (int t = 0; t < 60; t += minute_gap) {
					//cout << "[Time: " << m << "." << d << ' ' << h << ":" << t << "]" << endl;
					time_params = { m, d, h, t };
					handlerCore(time_params, argumentParser->getSunray(), &sdbk_handler);
				}

	// 4. 剔除无效定日镜
	return removeHeliostat();
}


void EnergyCalculatePipeline::initSolarScene(json& field_args) {
	// 1. 设置唯一Receiver
	//cout << "\t2.1 Create receiver" << endl;
	solar_scene->recvs = argumentParser->getReceivers();

	// 2. 设置Layout
	//cout << "\t2.2 Create layout" << endl;
	Layout* layout = LayoutCreator::getLayout(argumentParser->getLayoutType());
	solar_scene->layouts.push_back(layout);

	// 3. 设置Heliostat在Layout的布局
	//cout << "\t2.3 Create heliostats" << endl;
	solar_scene->layouts[0]->createHelioAndLayout(*argumentParser, field_args, solar_scene->helios);
	//solar_scene->saveSolarScene(argumentParser->getOutputPath());

	// 4. 设置镜面仿真模型
	//cout << "\t2.4 Initialize simulation model" << endl;
	solar_scene->setModelStatus(argumentParser->getModelType(), argumentParser->getCalcSigma());
}

void EnergyCalculatePipeline::handlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	// 1. 调整镜场定日镜朝向
	//cout << "\t3.1 Adjust heliostats' normal" << endl;
	Vector3d sunray_dir = sunray.changeSunRay(time_param);
	double DNI = sunray.calcDNI(time_param);
	solar_scene->changeHeliosNormal(sunray_dir);
	solar_scene->DNI = DNI;

	// 2. 计算阴影遮挡率
	//cout << "\t3.2 Calculate heliostats' shading and blocking factor" << endl;
	sdbk_handler->calcSceneShadowBlock();

	// 3. 计算镜场能量
	//cout << "\t3.3 Calculate field energy" << endl;
	handlerFunc(solar_scene, time_param, sunray, sdbk_handler);
}

void EnergyCalculatePipeline::handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	// 1. Get guass legendre calculation handler
	json gaussian_params = argumentParser->getConfig()["FluxParams"]["GaussianParams"].as<json>();
	int M = gaussian_params.get_with_default("M").as<int>();
	int N = gaussian_params.get_with_default("N").as<int>();
	int m = gaussian_params.get_with_default("m").as<int>();
	int n = gaussian_params.get_with_default("n").as<int>();
	GaussLegendre* gl_hander = GaussLegendre::getInstance(M, N);

	// 2. Calculate receiver flux integral
	bool calcCenterMode = false;
	if (solar_scene->getModelType() == bHFLCAL)
		calcCenterMode = true;
	ReceiverEnergyCalculator recv_energy_calc(solar_scene, gl_hander, m, n, calcCenterMode);
	recv_energy_calc.calcRecvEnergySum();
}

double EnergyCalculatePipeline::removeHeliostat() {
	sort(solar_scene->helios.begin(), solar_scene->helios.end(), [](Heliostat* a, Heliostat* b) {return a->energy > b->energy; });
	double total_sum = 0;
	for (int i = 0; i < solar_scene->layouts[0]->real_helio_num; ++i) {
		total_sum += solar_scene->helios[i]->energy;
	}
	return total_sum;
}


EnergyCalculatePipeline::~EnergyCalculatePipeline()
{
	delete solar_scene;
}

void FluxCalculatePipeline::handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	json flux = argumentParser->getConfig()["FluxParams"].as<json>();
	//vector<int> test_helio_index(flux["TestHelioIndex"].as<vector<int>>());
	vector<int> test_helio_index;
	for (int i = 0; i < solar_scene->helios.size(); ++i)
		test_helio_index.push_back(i);

	double DNI = sunray.calcDNI(time_param);
	string time_str = "M" + to_string(time_param[0]) + "D" + to_string(time_param[1])
		+ "H" + to_string(time_param[2]) + "m" + to_string(time_param[3]) + "/";

	if (!solar_scene->isCalcSigma()) {
		SigmaFitting sg_handler;
		vector<Heliostat*> rt_helio;
		vector<int> rt_helio_index(flux["RayTracingTestHelioIndex"].as<vector<int>>());
		for (auto& h_idx : rt_helio_index)
			rt_helio.push_back(solar_scene->helios[h_idx]);
		sg_handler.fromFluxPeak2Sigma(argumentParser->getInputPath() + flux["MaxFluxPath"].as_string() + time_str + "max_flux.txt", rt_helio, solar_scene->recvs[0], DNI);
	}

	sdbk_handler->setOutputPath(argumentParser->getOutputPath() + time_str);
	sdbk_handler->calcSceneFluxDistribution(test_helio_index, DNI, flux.get_with_default("GaussianParams").as<json>());
}
