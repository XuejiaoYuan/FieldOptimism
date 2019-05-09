#include "TestHandler.h"
#include <cmath>
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
		//field_args.push_back(new vector<double>{ 10, 10, 10 });			// 定日镜间隔
		//field_args.push_back(new vector<double>{ -80 });					// 第一行定日镜与接收器之间的距离
		//field_args.push_back(new vector<double>{ 100 });					// 定日镜行数
		//field_args.push_back(new vector<double>{ 100 });					// 定日镜列数

		field_args.push_back(new vector<double>{ 5, 5, 5 });			// 定日镜间隔
		field_args.push_back(new vector<double>{ -75 });					// 第一行定日镜与接收器之间的距离
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
	solar_scene->saveSolarScene("SdBkRes/");
}

void TestHandler::testSunRay()
{
	fstream outFile("SdBkRes/sunray_dir.txt", ios_base::out);
	Vector3d sunray_dir;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				outFile << sunray_dir.x() << ' ' << sunray_dir.y() << ' ' << sunray_dir.z() << endl;
			}
		}
	}
	outFile.close();
}

void TestHandler::testHelioPredict()
{
	Vector3d sunray_dir;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);
				SdBkCalcTest sdbk_test_handler(solar_scene);
				string path = "SdBkRes/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
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
				solar_scene->changeHeliosNormal(sunray_dir);
				SdBkCalcTest sdbk_test_handler(solar_scene);
				string path = "SdBkRes/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());
				rayCaster.rayCasting(solar_scene, path + "/", TestMode);
			}
		}
	}
}

void TestHandler::testPolygonClipping()
{
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	Vector3d sunray_dir;
	sdbk_calc->initBlockRelaIndex(sunray_dir);

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);

				string path = "SdBkRes/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());

				Timer::resetStart();
				sdbk_calc->calcTotalShadowBlock();
				auto t = Timer::getDuration();
				outFile.open("SdBkRes/PC_sdbk_t.txt", ios_base::app);
				outFile << t << endl;
				outFile.close();

				outFile.open(path + "/PC_sdbk.txt", ios_base::out);
				for (int i = 0; i < solar_scene->helios.size(); ++i)
					outFile << solar_scene->helios[i]->sd_bk << endl;
				outFile.close();
			}
		}
	}
}

void TestHandler::testTotalEnergyCalc(int M, int N, int m, int n)
{
	Vector3d sunray_dir;
	
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, M, N, m, n);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	sdbk_calc->initBlockRelaIndex(sunray_dir);
	
	string path = "SdBkRes/";
	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);
				

				// 1. CPU计算能量
				Timer::resetStart();
				double cpu_res = sdbk_calc->calcTotalEnergy(1);
				double cpu_min_t = Timer::getDuration();
					

				// 2. GPU计算能量
				Timer::resetStart();
				sdbk_calc->calcTotalShadowBlock();
				double gpu_res = e_handler.calcTotalEnergy(SunUpdateMode);
				double gpu_min_t = Timer::getDuration();

				
				outFile.open(path + "calc_energy_t.txt", ios_base::app);
				outFile << cpu_min_t << ' ' << gpu_min_t << endl;
				outFile.close();
				
				outFile.open(path + "calc_energy.txt", ios_base::app);
				outFile << cpu_res << ' ' << gpu_res << endl;
				outFile.close();
			}
		}
	}
}

void TestHandler::testHelioEnergyCalc(int M, int N, int m, int n) {
	Vector3d sunray_dir;

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, M, N, m, n);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	sdbk_calc->initBlockRelaIndex(sunray_dir);

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);

				sdbk_calc->calcTotalShadowBlock();
				
				//sdbk_calc->calcSingleFluxSum(0, 1);
				
				vector<float> gpu_res = e_handler.calcHelioEnergy(SunUpdateMode);
				double sum = e_handler.calcTotalEnergy(SunUpdateMode);

				string path = "SdBkRes/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());

				outFile.open(path + "/calc_helio_energy.txt", ios_base::out);
				outFile << sum << endl;
				for (auto&r : gpu_res)
					outFile << r << endl;
				outFile.close();

				outFile.open(path + "/sigma.txt", ios_base::out);
				for (auto&h : solar_scene->helios)
					outFile << h->sigma << endl;
				outFile.close();
			}
		}
	}

}

void TestHandler::testOneTimeHelioEnergy(int M, int N, int m, int n)
{
	Vector3d sunray_dir;
	float DNI = 1000;
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, M, N, m, n, true);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);

	fstream outFile;
	
	sunray_dir = sunray.changeSunRay(80.2, 185);		// mid day: 80.2, 185
	solar_scene->changeHeliosNormal(sunray_dir);
	sdbk_calc->setHFCALMode(true);
	sdbk_calc->initBlockRelaIndex(sunray_dir);
	sdbk_calc->calcTotalShadowBlock();

	SigmaFitting sg_handler;
	//sg_handler.fromFluxPeak2Sigma("SdBkRes/he_M3D21H8m0/max_flux.txt", solar_scene->helios, solar_scene->recvs[0], DNI);
	vector<int> test_helio_index = { 0, 980, 1960, 2113 };
	vector<Heliostat*> test_helio_list;
	test_helio_list = solar_scene->helios;
	//for (auto i : test_helio_index)
	//	test_helio_list.push_back(solar_scene->helios[i]);
	//sg_handler.fromFluxPeak2Sigma("Inputfiles/QMCRT/mid_day/iter_500_tri/max_flux.txt", test_helio_list, solar_scene->recvs[0], DNI);
	
	//DNI = 681;
	//sunray_dir = sunray.changeSunRay(9.5, 120);		// mid day: 80.2, 185
	//solar_scene->changeHeliosNormal(sunray_dir);
	//sdbk_calc->calcTotalShadowBlock();

	string path = "Outputfiles/QMCRT/mid_day/iter_500_tri_1/total_res/";
	_mkdir(path.c_str());

	//sdbk_calc->save_path = path;
	//for (auto h : test_helio_list)
	//	sdbk_calc->calcSingleFluxSum(h->helio_index, DNI);

	vector<float> gpu_res = e_handler.calcHelioEnergy(SunUpdateMode);
	double sum = e_handler.calcTotalEnergy(SunUpdateMode);

	outFile.open(path + "calc_helio_energy.txt", ios_base::out);
	//outFile << sum * DNI << endl;
	for (auto&r : gpu_res)
		outFile << r * DNI << endl;
	outFile.close();

	outFile.open(path + "sigma.txt", ios_base::out);
	for (auto&h : test_helio_list)
		outFile << h->sigma << endl;
	outFile.close();
}

void TestHandler::testFitSigme(int M, int N, int m, int n)
{
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, M, N, m, n);
	sunray.changeSunRay({ 3, 21, 8, 0 });
	Vector3d sunray_dir = sunray.changeSunRay(80.2, 185);
	solar_scene->changeHeliosNormal(sunray_dir);
	sdbk_calc->setHFCALMode(true);
	sdbk_calc->initBlockRelaIndex(sunray_dir);
	sdbk_calc->calcTotalShadowBlock();

	SigmaFitting sg_handler;
	sg_handler.fromFluxPeak2Sigma("SdBkRes/he_M3D21H8m0/max_flux.txt", solar_scene->helios, solar_scene->recvs[0], 1000);

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				//sunray_dir = sunray.changeSunRay(time_param);
				
				solar_scene->changeHeliosNormal(sunray_dir);
				sdbk_calc->calcTotalShadowBlock();

				string path = "SdBkRes/he_M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());

				vector<float> gpu_res = e_handler.calcHelioEnergy(SunUpdateMode);
				
				outFile.open(path + "/new_fit_helio_energy.txt", ios_base::out);
				for (auto&r : gpu_res)
					outFile << r * 1000 << endl;
				outFile.close();

				outFile.open(path + "/new_fit_sigma.txt", ios_base::out);
				for (auto&h : solar_scene->helios)
					outFile << h->sigma << endl;
				outFile.close();
			}
		}
	}

}

void TestHandler::testCalcPolygonCenter()
{
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	sunray.changeSunRay({ 3, 21, 8, 0 });
	Vector3d sunray_dir = sunray.changeSunRay(80.2, 185);
	solar_scene->changeHeliosNormal(sunray_dir);
	sdbk_calc->setHFCALMode(true);
	sdbk_calc->initBlockRelaIndex(sunray_dir);
	sdbk_calc->calcTotalShadowBlock();

	Paths subj(1), clips;
	vector<float2> v_list = { { 0, 0 },{ 0, 2 },{ 2, 2 },{ 2, 0 } };
	for (int i = 0; i < v_list.size(); ++i)
		subj[0] << IntPoint(v_list[i].x*VERTEXSCALE, v_list[i].y*VERTEXSCALE);

	Path clip1, clip2;
	clip1 << IntPoint(0, 0) << IntPoint(0, VERTEXSCALE) << IntPoint(0.5*VERTEXSCALE, VERTEXSCALE) << IntPoint(0.5*VERTEXSCALE, 0);
	clip2 << IntPoint(1.5*VERTEXSCALE, 0) << IntPoint(1.5*VERTEXSCALE, VERTEXSCALE) << IntPoint(2 * VERTEXSCALE, VERTEXSCALE) << IntPoint(2 * VERTEXSCALE,0);
	clips.push_back(clip1);
	clips.push_back(clip2);

	Clipper c;
	Paths solution;
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clips, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	c.Clear();
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(solution, ptClip, true);
	c.Execute(ctDifference, solution, pftNonZero, pftNonZero);
	
	Vector3d center;
	PolygonCenterCalculator::calcPolygonCenter(solution, center);
}
