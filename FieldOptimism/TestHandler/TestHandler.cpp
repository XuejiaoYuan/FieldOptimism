#include "TestHandler.h"

void TestHandler::initSolarScene(string scene_filepath, string sunray_filepath) {
	createSolarScene(scene_filepath);
	sunray.calcSunRay(sunray_filepath);
}


void TestHandler::createSolarScene(string scene_filepath)
{
	solar_scene = new SolarScene();
	solar_scene->initFieldParam(scene_filepath);

	vector<vector<double>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<double>{ 10, 10, 10 });			// 定日镜间隔
		field_args.push_back(new vector<double>{ -80 });					// 第一行定日镜与接收器之间的距离
		field_args.push_back(new vector<double>{ 100 });					// 定日镜行数
		field_args.push_back(new vector<double>{ 100 });					// 定日镜列数
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<double>{ double(0.1) });			// 定日镜包围盒安全距离
		field_args.push_back(new vector<double>{ double(114.9417) });	// 第一个同心圆与接收器之间的距离(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });				// 第一个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 1.4f });				// 第二个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 2.0f });					// 第三个同心环中定日镜的分布间隔
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);
}

void TestHandler::testHelioPredict()
{
	Vector3d sunray_dir;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeSolarScene(sunray_dir);
				SdBkCalcTest sdbk_test_handler(solar_scene);
				string path = "SdBkRes/MM" + to_string(day.x) + "DD" + to_string(day.y) + "HH" + to_string(hour) + "mm" + to_string(minute);
				_mkdir(path.c_str());
				sdbk_test_handler.totalHeliosTest(path + "/");
			}
		}
	}
}

void TestHandler::testRayCasting()
{
	Vector3d sunray_dir;
	DiscreteRayCaster rayCaster;

	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeSolarScene(sunray_dir);
				SdBkCalcTest sdbk_test_handler(solar_scene);
				string path = "SdBkRes/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute) + "/";
				_mkdir(path.c_str());
				rayCaster.rayCasting(solar_scene, path, TestMode);
			}
		}
	}
}

void TestHandler::testEnergyCalc(int M, int N, int m, int n)
{
	Vector3d sunray_dir;
	
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, M, N, m, n);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	
	fstream outFile("SdBkRes/energy_calc_t_g" + to_string(M) + "x" + to_string(N) + "_n" + to_string(m) + "x" + to_string(n) + ".txt", ios_base::out);
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeSolarScene(sunray_dir);
				
				// 1. CPU计算能量
				Timer::resetStart();
				sdbk_calc->calcTotalEnergy(1);
				outFile << Timer::getDuration() << ' ';

				// 2. GPU计算能量
				Timer::resetStart();
				sdbk_calc->calcTotalShadowBlock();
				e_handler.calcHelioEnergy(1.31, SunUpdateMode);
				outFile << Timer::getDuration() << endl;
			}
		}
	}
	outFile.close();
}


