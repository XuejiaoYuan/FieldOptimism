
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSF.h"
#include "./DiscreteRayCaster/DiscreteRayCaster.h"
#include "./DataStructure/Timer.h"
#include "./HelioEnergy/HelioEnergy.cuh"


int main(int argc, char** argv) {

	// 计算阴影与遮挡参数
	// scene_file 镜场布置文件
	// sunray_file 太阳参数文件
	// options: -a_c 计算每个时刻下所有镜子的阴影遮挡结果(clipper);
	//			-a_r 计算每个时刻下所有镜子的阴影遮挡结果（raytracing）;
	//			-s_l 计算每个时刻下所有镜子的阴影遮挡结果（LSPIA）
	// outfile_options: -f 存储计算阴影遮挡的结果(clipper)
	//					-o 不存储结果
	// save_path 输出文件的存储地址
	if (argc != 6) {
		cout << argc << endl;
		cout << "Usage: [scene_file] [sunray_file] [options] [outfile_options] [save_path]" << endl;
		return -1;
	}


	string scene_filepath = string(argv[1]);
	string sunray_filepath = string(argv[2]);
	string options = string(argv[3]);
	string outfile_options = string(argv[4]);
	string save_path = string(argv[5]);

	SunRay sunray;
	Vector3d sunray_dir = sunray.calcSunRay(sunray_filepath);
	MatrixXd a;

	SolarScene* solar_scene = new SolarScene();
	solar_scene->initFieldParam(scene_filepath);

	// TODO 镜场参数优化
	vector<vector<double>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<double>{10, 10, 10});			// 定日镜间隔
		field_args.push_back(new vector<double>{ -80 });					// 第一行定日镜与接收器之间的距离
		field_args.push_back(new vector<double>{ 100 });					// 定日镜行数
		field_args.push_back(new vector<double>{ 100 });					// 定日镜列数
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<double>{ double(0.1) });			// 定日镜包围盒安全距离
		field_args.push_back(new vector<double>{ double( 114.9417 )});	// 第一个同心圆与接收器之间的距离(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });				// 第一个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 1.4f });				// 第二个同心环中定日镜的分布间隔
		field_args.push_back(new vector<double>{ 2.0f});					// 第三个同心环中定日镜的分布间隔
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);

	// test gpu dda
	sunray_dir = sunray.changeSunRay({ 3, 21, 8, 0 });
	Timer::resetStart();
	solar_scene->changeSolarScene(sunray_dir);
	Timer::printDuration("change solar scene");

	//DiscreteRayCaster rayCaster;
	//rayCaster.rayCasting(solar_scene);

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene);
	HelioEnergy e_handler(solar_scene, 8, 8, 3, 3);
	sdbk_calc->gl = new GaussLegendreCPU(8, 8);

	//vector<vector<int>> daily = {
	//	{ 3,21 },		// 春分
	//	{ 6,22 },		// 夏至
	//	{ 9,23 },		// 秋分
	//	{ 12,22 }		// 冬至
	//};
	fstream outFile("energy_calc_t_g3x3_n8x8_3.txt", ios_base::out);
	vector<int> sample_hour = { 8, 10, 12, 14, 16 };
	vector<int> time_param(4);
	for (int month = 1; month <= 12; month++) {
		for (int day = 1; day < 29; day += 3) {
			for (int hour = 8; hour < 17; hour++) {
				for (int min = 0; min < 30; min += 10) {
					cout << "Calc start" << endl;
					time_param[0] = month;
					time_param[1] = day;
					time_param[2] = hour;
					time_param[3] = min;
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
	}
	outFile.close();


	// 多种阴影遮挡预处理测试
	//vector<int> time_param = { 1,22, 8, 0 };
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);
	SdBkCalcTest sdbk_test_handler(solar_scene);
	sdbk_test_handler.totalHeliosTest("SdBkRes");
	

	return 1;
}