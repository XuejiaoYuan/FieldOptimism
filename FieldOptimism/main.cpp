#include "TestHandler/TestHandler.h"

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


	TestHandler t_handler(30);
	t_handler.initSolarScene(scene_filepath, sunray_filepath);

	// 1. ��Ӱ�ڵ�����
	//t_handler.testRayCasting();
	//t_handler.testHelioPredict();
	//t_handler.testPolygonClipping();

	// 2. ��������
	t_handler.testEnergyCalc(8, 8, 2, 2);

 	return 1;
}