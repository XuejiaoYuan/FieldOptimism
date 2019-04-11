
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSF.h"
#include "./DiscreteRayCaster/DiscreteRayCaster.h"
#include "./DataStructure/Timer.h"
#include "./HelioEnergy/HelioEnergy.cuh"


int main(int argc, char** argv) {

	// ������Ӱ���ڵ�����
	// scene_file ���������ļ�
	// sunray_file ̫�������ļ�
	// options: -a_c ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ����(clipper);
	//			-a_r ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ������raytracing��;
	//			-s_l ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ������LSPIA��
	// outfile_options: -f �洢������Ӱ�ڵ��Ľ��(clipper)
	//					-o ���洢���
	// save_path ����ļ��Ĵ洢��ַ
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

	// TODO ���������Ż�
	vector<vector<double>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<double>{10, 10, 10});			// ���վ����
		field_args.push_back(new vector<double>{ -80 });					// ��һ�ж��վ��������֮��ľ���
		field_args.push_back(new vector<double>{ 100 });					// ���վ�����
		field_args.push_back(new vector<double>{ 100 });					// ���վ�����
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<double>{ double(0.1) });			// ���վ���Χ�а�ȫ����
		field_args.push_back(new vector<double>{ double( 114.9417 )});	// ��һ��ͬ��Բ�������֮��ľ���(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });				// ��һ��ͬ�Ļ��ж��վ��ķֲ����
		field_args.push_back(new vector<double>{ 1.4f });				// �ڶ���ͬ�Ļ��ж��վ��ķֲ����
		field_args.push_back(new vector<double>{ 2.0f});					// ������ͬ�Ļ��ж��վ��ķֲ����
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
	//	{ 3,21 },		// ����
	//	{ 6,22 },		// ����
	//	{ 9,23 },		// ���
	//	{ 12,22 }		// ����
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
					
					// 1. CPU��������
					Timer::resetStart();
					sdbk_calc->calcTotalEnergy(1);
					outFile << Timer::getDuration() << ' ';
					
					// 2. GPU��������
					Timer::resetStart();
					sdbk_calc->calcTotalShadowBlock();
					e_handler.calcHelioEnergy(1.31, SunUpdateMode);
					outFile << Timer::getDuration() << endl;
				}
			}
		}
	}
	outFile.close();


	// ������Ӱ�ڵ�Ԥ�������
	//vector<int> time_param = { 1,22, 8, 0 };
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);
	SdBkCalcTest sdbk_test_handler(solar_scene);
	sdbk_test_handler.totalHeliosTest("SdBkRes");
	

	return 1;
}