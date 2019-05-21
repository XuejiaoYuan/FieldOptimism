#include "TestHandler.h"
#include <cmath>

void TestHandler::initSolarScene(string scene_filepath, string sunray_filepath, string output_filepath) {
	createSolarScene(scene_filepath, output_filepath);
	sunray.calcSunRay(sunray_filepath);
}

void TestHandler::createSolarScene(string scene_filepath, string output_filepath)
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
		field_args.push_back(new vector<double>{ double(114.9417) });		// 第一个同心圆与接收器之间的距离(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });					// 第一个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 1.4f });					// 第二个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 2.0f });					// 第三个同心环中定日镜的分布间隔
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);
	
	vector<string> field_path = { "crossRectField/", "rectField/", "fermatField/" };
	this->output_filepath = output_filepath + field_path[solar_scene->layouts[0]->layout_type];
	string save_path = output_filepath + "FieldFiles/" + field_path[solar_scene->layouts[0]->layout_type];
	_mkdir(this->output_filepath.c_str());
	_mkdir(save_path.c_str());
	
	solar_scene->saveSolarScene(save_path);
}

void TestHandler::testSunRay()
{
	string save_path = output_filepath + "SunRayDir/";
	_mkdir(save_path.c_str());
	
	Vector3d sunray_dir;
	float DNI;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				DNI = sunray.calcDNI(time_param);

				string f_name = "he_scene_T" + to_string(solar_scene->layouts[0]->layout_type) + "_H" + to_string(solar_scene->helios.size());
				string time_str = "M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				fstream outFile(save_path + f_name + "_"  + time_str + ".json", ios_base::out);
				outFile << "sun_dir " << sunray_dir.x() << ' ' << sunray_dir.y() << ' ' << sunray_dir.z() << endl;
				outFile << "csr 0.1" << endl;
				outFile << "dni " << DNI << endl;
				outFile << "num_of_sunshape_groups 128\n"
					<< "num_per_sunshape_group 2\n"
					<< "inverse_transform_sampling_groups 64\n"
					<< "receiver_pixel_length 0.05\n"
					<< "helio_disturb_std 0.001\n"
					<< "helio_reflected_rate 0.88\n"
					<< "helio_pixel_length 0.01" << endl;
				outFile.close();

				outFile.open(save_path + f_name + "_run.bat", ios_base::out | ios_base::app);
				outFile << "SolarEnergyRayTracing.exe "
					<< "-c InputFiles/he_scene_T1_H9950/DiffSunDir/he_scene_T1_H9950_" + time_str + ".json " 
					<< "-s InputFiles/he_scene_T1_H9950/he_scene_T1_H9950.scn -h "
					<< "InputFiles/he_scene_T1_H9950/he_scene_T1_H9950.txt "
					<< "-o OutputFiles/he_scene_T1_H9950/" + time_str + "/\n" << endl;
				outFile.close();
			}
		}
	}
	
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
				string path = output_filepath + "HelioPredict/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
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
				string path = output_filepath + "RayCasting/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
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
	sdbk_calc->initBlockRelaIndex();

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);

				string path = output_filepath + "PolygonClipping/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());

				Timer::resetStart();
				sdbk_calc->calcTotalShadowBlock();
				auto t = Timer::getDuration();
				outFile.open(path + "/PC_sdbk_t.txt", ios_base::app);
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
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	sdbk_calc->initBlockRelaIndex();
	
	fstream outFile;
	string save_path = output_filepath + "TotalEnergy/";
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

				outFile.open(output_filepath + "TotalEnergy/calc_energy_cpu_t.txt", ios_base::app);
				outFile << cpu_min_t << endl;
				outFile.close();
				
				outFile.open(output_filepath + "TotalEnergy/calc_energy_cpu.txt", ios_base::app);
				outFile << cpu_res << endl;
				outFile.close();
			}
		}
	}
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);

				// 2. GPU计算能量
				Timer::resetStart();
				sdbk_calc->calcTotalShadowBlock();
				double gpu_res = EnergyCalcPipeline::calcFieldEnergySum(solar_scene, M, N, m, n, false);
				double gpu_min_t = Timer::getDuration();

				outFile.open(output_filepath + "TotalEnergy/calc_energy_gpu_t.txt", ios_base::app);
				outFile << gpu_min_t << endl;
				outFile.close();
				
				outFile.open(output_filepath + "TotalEnergy/calc_energy_gpu.txt", ios_base::app);
				outFile << gpu_res << endl;
				outFile.close();
			}
		}
	}


}

void TestHandler::testHelioEnergyCalc(int M, int N, int m, int n) {
	Vector3d sunray_dir;

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	sdbk_calc->initBlockRelaIndex();

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				sunray_dir = sunray.changeSunRay(time_param);
				solar_scene->changeHeliosNormal(sunray_dir);

				sdbk_calc->calcTotalShadowBlock();
				
				//sdbk_calc->calcSingleFluxSum(0, 1);
				
				//vector<float> gpu_res = e_handler.calcHelioEnergy(SunUpdateMode);
				//double sum = e_handler.calcTotalEnergy(SunUpdateMode);
				double sum = EnergyCalcPipeline::calcFieldEnergySum(solar_scene, M, N, m, n, false);

				string path = output_filepath + "HelioEnergy/M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute);
				_mkdir(path.c_str());

				outFile.open(path + "/calc_helio_energy.txt", ios_base::out);
				outFile << sum << endl;
				//for (auto&r : gpu_res)
				//	outFile << r << endl;
				//outFile.close();

				outFile.open(path + "/sigma.txt", ios_base::out);
				for (auto&h : solar_scene->helios)
					outFile << h->sigma << endl;
				outFile.close();
			}
		}
	}

}

void TestHandler::testHelioFluxDensityModel(int M, int N, int m, int n)
{
	vector<int> test_helio_index = { 75, 1699, 3060, 4960, 6000};
	vector<Heliostat*> helio_list4calc;
	for (auto i : test_helio_index)
		helio_list4calc.push_back(solar_scene->helios[i]);
	
	string path = output_filepath + "ModelTest/LWRatioCompare/";
	_mkdir(path.c_str());

	vector<Heliostat*> test_helio_list;
	for (auto i : test_helio_index)
		test_helio_list.push_back(solar_scene->helios[i]);
	
	Vector3d sunray_dir;
	float DNI = 1;
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	sdbk_calc->gl = new GaussLegendreCPU(M, N, m, n);
	sdbk_calc->initBlockRelaIndex();
	sdbk_calc->setHFCALMode(true);
	
	bool calcLWRatio = true;
	bool calcSigma = false;	

	fstream outFile;
	for (auto& day : sample_day) {
		for (auto& hour : sample_hour) {
			for (int minute = 0; minute < 60; minute += minute_inter) {
				vector<int> time_param = { day.x, day.y, hour, minute };
				//vector<int> time_param = { 6, 21, 8, 15 };
				string time_str = "M" + to_string(day.x) + "D" + to_string(day.y) + "H" + to_string(hour) + "m" + to_string(minute) + "/";
				sunray_dir = sunray.changeSunRay(time_param);
				DNI = sunray.calcDNI(time_param);
				solar_scene->changeHeliosNormal(sunray_dir, calcLWRatio, calcSigma);
				
				outFile.open(path + "lw_ratio_m11.txt", ios_base::out | ios_base::app);
				for (auto& h : test_helio_list) {	
					outFile << h->l_w_ratio << ' ';
				}
				outFile << endl;
				outFile.close();
				
				sdbk_calc->calcTotalShadowBlock();

				SigmaFitting sg_handler;
				sg_handler.fromFluxPeak2Sigma("Inputfiles/QMCRT/ModelTest/he_scene_T1_H9950/" + time_str + "max_flux.txt", helio_list4calc, solar_scene->recvs[0], DNI);

				sdbk_calc->save_path = path + time_str;
				_mkdir(sdbk_calc->save_path.c_str());
				for (auto h : test_helio_list)
					sdbk_calc->calcSingleFluxSum(h->helio_index, DNI);
			}
		}
	}

}

void TestHandler::testOneTimeHelioEnergy(int M, int N, int m, int n)
{
	bool calcLWRatio = true;
	bool calcSigma = true;

	// 1. Initialize solar scene
	vector<int> time_params = { 3, 21, 10, 0 };
	Vector3d sunray_dir = sunray.changeSunRay(time_params);
	double DNI = sunray.calcDNI(time_params);
	solar_scene->changeHeliosNormal(sunray_dir, calcLWRatio, calcSigma);

	// 2. Calculate shadow and block ratio
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	sdbk_calc->setHFCALMode(calcLWRatio);
	sdbk_calc->initBlockRelaIndex();
	sdbk_calc->calcTotalShadowBlock();

	// 3. Calculate total energy
	double gpu_res = EnergyCalcPipeline::calcFieldEnergySum(solar_scene, M, N, m, n, calcLWRatio);
	cout << "Total Energy: " << gpu_res * DNI << endl;
 }


void TestHandler::testCalcPolygonCenter()
{
	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	sunray.changeSunRay({ 3, 21, 8, 0 });
	Vector3d sunray_dir = sunray.changeSunRay(80.2, 185);
	solar_scene->changeHeliosNormal(sunray_dir);
	sdbk_calc->setHFCALMode(true);
	sdbk_calc->initBlockRelaIndex();
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
