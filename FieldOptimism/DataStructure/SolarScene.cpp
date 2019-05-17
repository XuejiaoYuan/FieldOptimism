//
// Created by Amber on 2018/4/3.
//

#include "SolarScene.h"


bool SolarScene::changeHeliosNormal(const Vector3d & sunray_dir, bool calcLWRatio, bool calcSimga)
{
	this->sunray_dir = Vector3d(sunray_dir.x(), sunray_dir.y(), sunray_dir.z() );

#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++)
		helios[i]->changeSurfaceNormal(sunray_dir, calcLWRatio, calcSimga);

	layouts[0]->setHelioLayout(helios);
	
	return true;
}

//
// [镜场优化预处理] 镜场优化中固定参数设置
//
bool SolarScene::initFieldParam(const string& file_name)
{
	fstream inFile(file_name, ios_base::in);
	if (inFile.fail()) {
		cerr << "Can't open the filed parameter file!" << endl;
		return false;
	}

	string line, word;
	stringstream line_stream;
	InputMode input_mode = Initial;
	int row = 0;
	int col = 0;
	Layout* layout;
	int helio_type;
	int layout_type;
	Vector2d helio_gap;
	Vector2i helio_matrix;
	Vector3d helio_size;
	Vector3d helio_pos;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;

		if (word == "#") {
			line_stream >> word;
			if (word == "Receiver") {
				input_mode = ReceiverMode;
				continue;
			}
			else if (word == "Heliostats") {
				input_mode = HeliostatMode;
				continue;
			}
			else if (word == "Grid") {
				input_mode = LayoutMode;
				continue;
			}
		}

		switch (input_mode)
		{
		case ReceiverMode: {
			int recv_type;
			ReceiverCreator recv_creator;
			line_stream >> recv_type;
			Receiver* recv = recv_creator.getReceiver((ReceiverType)recv_type);
			recv->init_recv(inFile, input_mode);
			recvs.push_back(recv);
			break;
		}
		case HeliostatMode: {
			if (word == "gap")
				line_stream >> helio_gap.x() >> helio_gap.y();
			else if (word == "matrix")
				line_stream >> helio_matrix.x() >> helio_matrix.y();
			else if (word == "pos")
				line_stream >> helio_pos.x() >> helio_pos.y() >> helio_pos.z();
			else if (word == "size")
				line_stream >> helio_size.x() >> helio_size.y() >> helio_size.z();
			else if (word == "end")
				input_mode = Initial;
			break;
		}
		case LayoutMode: {
			if (word == "Scene") {
				line_stream >> layout_type;
				LayoutCreator layout_creator;
				layout = layout_creator.getLayout((LayoutType)layout_type);
			}
			else if (word == "row")
				line_stream >> row;
			else if (word == "col")
				line_stream >> col;
			else if (word == "type")
				line_stream >> helio_type;
			else if (word == "end") {
				input_mode = Initial;
				layouts.push_back(layout);
				layout = nullptr;
			}
			break;
		}
		case Initial:
			break;
		default:
			break;
		}
	}
	inFile.close();
	layouts[0]->helio_type = (HelioType)helio_type;
	layouts[0]->helio_pos = helio_pos;
	layouts[0]->helio_gap = helio_gap;
	layouts[0]->helio_matrix = helio_matrix;
	layouts[0]->helio_size = helio_size;

	return true;
}

//
// [镜场优化] 镜场优化时改变优化参数以调整镜场位置
//
bool SolarScene::adjustFieldParam(const vector<vector<double>*>& field_args)
{
	if (!helios.empty()) {
		for (auto&h : helios)
			delete h;
	}
	helios.clear();
	layouts[0]->adjustHelioLayout(helios, field_args, recvs);
	return true;
}


void SolarScene::saveSolarScene(string scene_savepath)
{
	fstream outFile(scene_savepath +"scene_T" + to_string(layouts[0]->layout_type) + "_H" + to_string(helios.size()) + ".scn", ios::out);
	if (outFile.fail()) {
		cerr << "Can't write to this file!" << endl;
	}

	outFile << "# Ground Boundary" << endl;
	outFile << "ground " << scene_length << ' ' << scene_width << endl;
	outFile << "ngrid 1"  << endl;

	outFile << "\n# Receiver attributes" << endl;
	for (auto&recv : recvs) {
		outFile << "Recv " << recv->recv_type << endl;
		outFile << "pos " << recv->recv_pos.x() << ' ' << recv->recv_pos.y() << ' ' << recv->recv_pos.z() << endl;
		outFile << "size " << recv->recv_size.x() << ' ' << recv->recv_size.y() << ' ' << recv->recv_size.z() << endl;
		outFile << "norm " << recv->recv_normal.x() << ' ' << recv->recv_normal.y() << ' ' << recv->recv_normal.z() << endl;
		outFile << "face 0" << endl;
	}
	outFile << "end" << endl;

	for (int i = 1; i <= layouts.size(); i++) {
		outFile << "\n# Grid" << i << " attributes" << endl;
		Layout* layout = layouts[i - 1];
		outFile << "Grid 0" << endl;
		outFile << "pos " << layout->layout_bound_pos.x() << ' ' << layout->layout_bound_pos.y() << ' ' << layout->layout_bound_pos.z() << endl;
		outFile << "size " << layout->layout_size.x() << ' ' << layout->layout_size.y() << ' ' << layout->layout_size.z() << endl;
		outFile << "inter " << layout->helio_interval.x() << ' ' << layout->helio_interval.y() << ' ' << layout->helio_interval.z() << endl;
		outFile << "n " << helios.size() << endl;
		outFile << "type 0" << endl;
		outFile << "end" << endl;
		

		outFile << "\n# Heliostats" << endl;
		float2 gap;
		int2 matrix;
		outFile << "gap " << helios[0]->helio_gap.x() << ' ' << helios[0]->helio_gap.y() << endl;
		outFile << "matrix " << helios[0]->helio_matrix.x() << ' ' << helios[0]->helio_matrix.y() << endl;


		for (auto&helio : helios) {
			outFile << "helio " << helio->helio_pos.x() << ' ' << helio->helio_pos.y() << ' ' << helio->helio_pos.z() << endl;
			outFile << helio->helio_size.x() << ' ' << helio->helio_size.y() << ' ' << helio->helio_size.z() << endl;

		}
	}
	outFile.close();

	outFile.open(scene_savepath + "scene_T" + to_string(layouts[0]->layout_type) + "_H" + to_string(helios.size()) + ".txt", ios_base::out);
	outFile << helios.size() << endl;
	for (int i = 0; i < helios.size(); ++i)
		outFile << i << ' ';
	outFile << endl;
	outFile.close();
}
