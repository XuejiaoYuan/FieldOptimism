#include "EnergyCalculatePipeline.h"
#include <direct.h>


//
// [能量计算接口] 镜场优化调用
//
double EnergyCalculatePipeline::handler(json& field_args) {
	// 1. 初始化镜场
	initSolarScene(field_args);
	
	double res = traverseTimeHandler();

	return res;
}

//
// [能量计算接口] 固定镜场能量计算调用
//
double EnergyCalculatePipeline::handler()
{
	// 1. 初始化镜场
	solar_scene->initSolarScene(argumentParser->getConfig()["Path"]["ScnPath"].as<string>());
	solar_scene->setModelStatus(argumentParser->getModelType(), argumentParser->getCalcSigma());

	double res = traverseTimeHandler();
	
	saveHeliostatEnergy();
	return res;
}

//
// [能量计算接口] 计算全镜场各定日镜的总能量
//
void EnergyCalculatePipeline::saveHeliostatEnergy() {
	fstream out(argumentParser->getOutputFn(), ios_base::out);
	vector<Heliostat*> tmp = solar_scene->helios;
	sort(tmp.begin(), tmp.end(), [](Heliostat* a, Heliostat* b) {return a->energy > b->energy; });
	for (int i = 0; i <argumentParser->getNumOfHelio(); ++i) {
		Heliostat* h = tmp[i];
		out << h->helio_pos.x() << ' ' << h->helio_pos.y() << ' ' << h->helio_pos.z() << ' ' << setprecision(16) << h->energy << endl;
	}
	out.close();
}

//
// [能量计算接口] 保存优化镜场结果
//
void EnergyCalculatePipeline::saveSolarScene()
{
	sort(solar_scene->helios.begin(), solar_scene->helios.end(), [](Heliostat* a, Heliostat* b) { return a->energy > b->energy; });
	solar_scene->helios.resize(argumentParser->getNumOfHelio());
	solar_scene->saveSolarScene(argumentParser->getOutputFn());
}

//
// [能量计算接口] 按采样时间计算镜场能量
//
double EnergyCalculatePipeline::traverseTimeHandler() {
	// 2. 设置采样时间
	vector<int> months = argumentParser->getTimeParams(MonthType);
	vector<int> days = argumentParser->getTimeParams(DayType);
	vector<int> hours = argumentParser->getTimeParams(HourType);
	int minute_gap = argumentParser->getMinuteGap();

	// 3. 计算采样时刻镜场能量
	//cout << "3. Start simulate flux density distribution" << endl;
	vector<int> time_params;
	SdBkCalc sdbk_handler(solar_scene);
	bool status = true;
	for (auto& m : months)
		for (auto& d : days)
			for (auto& h : hours)
				for (int t = 0; t < 60; t += minute_gap) {
					//cout << "[Time: " << m << "." << d << ' ' << h << ":" << t << "]" << endl;
					time_params = { m, d, h, t };
					status = oneTimeHandlerCore(time_params, argumentParser->getSunray(), &sdbk_handler);
					if (!status) break;
				}

	// 4. 剔除无效定日镜
	if (!status)
		return 0;
	return removeHeliostat();
}

//
// [能量计算接口] 创建新镜场
//
void EnergyCalculatePipeline::initSolarScene(json& field_args) {
	// 1. 设置唯一Receiver
	solar_scene->recvs = argumentParser->getReceivers();

	// 2. 设置Layout
	Layout* layout = LayoutCreator::getLayout(argumentParser->getLayoutType());
	solar_scene->layouts.push_back(layout);

	// 3. 设置镜面仿真模型
	solar_scene->setModelStatus(argumentParser->getModelType(), argumentParser->getCalcSigma());

	// 4. 设置Heliostat在Layout的布局
	solar_scene->layouts[0]->createHelioAndLayout(*argumentParser, field_args, solar_scene->helios);
	//solar_scene->saveSolarScene(argumentParser->getOutputPath());

}

// 
// [能量计算接口] 计算某采样时刻下的镜场能量
//
bool EnergyCalculatePipeline::oneTimeHandlerCore(vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	// 1. 调整镜场定日镜朝向
	Vector3d sunray_dir = sunray.changeSunRay(time_param);
	double DNI = sunray.calcDNI(time_param);
	solar_scene->changeHeliosNormal(sunray_dir);
	solar_scene->DNI = DNI;

	// 2. 计算阴影遮挡率
	bool status = sdbk_handler->calcSceneShadowBlock();

	// 3. 计算镜场能量
	if(status)
		handlerFunc(solar_scene, time_param, sunray, sdbk_handler);

	return status;
}

//
// [能量计算接口] 计算镜场采样时刻下的镜场能量核心函数
//
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

//
//[能量计算接口] 剔除多余定日镜
//
double EnergyCalculatePipeline::removeHeliostat() {
	vector<Heliostat*> tmp = solar_scene->helios;
	sort(tmp.begin(), tmp.end(), [](Heliostat* a, Heliostat* b) {return a->energy > b->energy; });
	double total_sum = 0;
	for (int i = 0; i < solar_scene->layouts[0]->real_helio_num; ++i) {
		total_sum += tmp[i]->energy;
	}
	return total_sum;
}


EnergyCalculatePipeline::~EnergyCalculatePipeline()
{
	delete solar_scene;
}

//
// [定日镜辐射能计算接口] 计算测试定日镜的辐射能密度分布结果
//
void HelioFluxCalculatePipeline::handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	json flux = argumentParser->getConfig()["FluxParams"].as<json>();
	vector<int> test_helio_index(flux["TestHelioIndex"].as<vector<int>>());

	double DNI = sunray.calcDNI(time_param);
	string time_str = "M" + to_string(time_param[0]) + "D" + to_string(time_param[1])
		+ "H" + to_string(time_param[2]) + "m" + to_string(time_param[3]) + "/";

	// 若有光线跟踪峰值，则使用单点拟合计算sigma
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


//
// [镜场辐射能计算接口] 计算全镜场在接收器平面的辐射能密度分布
//
void FieldFluxCalculatePipeline::handlerFunc(SolarScene* solar_scene, vector<int>& time_param, SunRay& sunray, SdBkCalc* sdbk_handler) {
	// 1. Get guass legendre calculation handler
	json flux = argumentParser->getConfig()["FluxParams"].as<json>();
	json gaussian_params = flux["GaussianParams"].as<json>();
	int M = gaussian_params.get_with_default("M").as<int>();
	int N = gaussian_params.get_with_default("N").as<int>();
	int m = gaussian_params.get_with_default("m").as<int>();
	int n = gaussian_params.get_with_default("n").as<int>();
	GaussLegendre* gl_hander = GaussLegendre::getInstance(M, N);

	// 2. Check if calculate sigma
	double DNI = sunray.calcDNI(time_param);
	string time_str = "M" + to_string(time_param[0]) + "D" + to_string(time_param[1])
		+ "H" + to_string(time_param[2]) + "m" + to_string(time_param[3]) + "/";

	if (!solar_scene->isCalcSigma()) {
		SigmaFitting sg_handler;
		sg_handler.fromFluxPeak2Sigma(argumentParser->getInputPath() + flux["MaxFluxPath"].as_string() + time_str + "max_flux.txt", solar_scene->helios, solar_scene->recvs[0], DNI);
	}

	// 2. Calculate receiver discrete flux 
	bool calcCenterMode = false;
	if (solar_scene->getModelType() == bHFLCAL)
		calcCenterMode = true;
	ReceiverEnergyCalculator recv_energy_calc(solar_scene, gl_hander, m, n, calcCenterMode);
	vector<float> res = recv_energy_calc.calcRecvFluxDistribution();	

	// 3. Save receiver discrete flux 
	int recvFaceNum = solar_scene->recvs[0]->recv_face_num;
	Vector2i row_col = solar_scene->recvs[0]->rows_cols;
	int row = row_col.x(), col = row_col.y();
	for (int k = 0; k < recvFaceNum; ++k) {
		string path = argumentParser->getOutputPath() + time_str;
		_mkdir(path.c_str());
		fstream out(path + "R" + to_string(k) + ".txt", ios_base::out);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j)
				out << res[k*row*col + i*col + j] * DNI << ' ';
			out << endl;
		}
		out.close();
	}
}
