//
// Created by Amber on 2018/4/3.
//

#include "SolarScene.h"


bool SolarScene::changeHeliosNormal(const Vector3d & sunray_dir)
{
	this->sunray_dir = Vector3d(sunray_dir.x(), sunray_dir.y(), sunray_dir.z());

#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++)
		helios[i]->changeSurfaceNormal(sunray_dir, model_type, calc_sigma);

	layouts[0]->storeHelioToLayout(helios);
	
	return true;
}

void SolarScene::saveSolarScene(string scene_savepath)
{
	fstream outFile(scene_savepath +"scene_T" + to_string(layouts[0]->layout_type) + "_H" + to_string(helios.size()) + ".scn", ios::out);
	if (outFile.fail()) {
		cerr << "Can't write to this file!" << endl;
	}

	outFile << "# Ground Boundary" << endl;
	outFile << "ground 0 0" << endl;
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
		outFile << "pos " << layout->layout_bound_pos.x() << ' ' << 0 << ' ' << layout->layout_bound_pos.y() << endl;
		// TODO: There may some error for QMCRT
		outFile << "size " << layout->layout_size.x() << ' ' << layout->helio_interval.x() << ' ' << layout->layout_size.y() << endl;
		outFile << "inter " << layout->helio_interval.x() << ' ' << layout->helio_interval.x() << ' ' << layout->helio_interval.y() << endl;
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

void SolarScene::saveHeSolarScene(string scene_savepath) {
	fstream outFile(scene_savepath + "he_scene_T" + to_string(layouts[0]->layout_type) + "_H" + to_string(helios.size()) + ".scn", ios::out);
	if (outFile.fail()) {
		cerr << "Can't write to this file!" << endl;
	}

	outFile << "// // Scene area" << endl;
	outFile << "-431602080.000000 -431602080.000000\n" << endl;

	outFile << "// Receiver position" << endl;
	outFile << "0.000000 110.000000 0.000000 12.000000 12.000000 1.000000 0\n" << endl;

	outFile << "// Field layout type" << endl;
	outFile << "0 rectangular\n" << endl;

	outFile << "// heliostat intervals" << endl;
	outFile << "5.000000 5.000000\n" << endl;

	outFile << "// row  col\n" << endl;

	outFile << "// Heliostat slice matrix" << endl;
	outFile << "1 1\n" << endl;

	outFile << "// Heliostat slice gap" << endl;
	outFile << "0.200000 0.200000\n" << endl;

	outFile << "// Reflector numbers" << endl;
	outFile << helios.size() << "\n" << endl;

	outFile << "// Reflector position" << endl;

	for (auto&helio : helios) {
		outFile << helio->helio_pos.x() << ' ' << helio->helio_pos.y() << ' ' << -helio->helio_pos.z() << ' '
			<< helio->helio_size.x() << ' ' << helio->helio_size.z() << ' ' << helio->helio_size.y() << endl;
	}

	outFile.close();
}

void SolarScene::setModelStatus(const ModelType & model_type, const bool & calc_sigma)
{
	this->model_type = model_type;
	this->calc_sigma = calc_sigma;
}

SolarScene::~SolarScene() {
	for (auto& layout : layouts) {
		delete layout;
		layout = NULL;
	}
	layouts.clear();

	for (auto& helio : helios) {
		delete helio;
		helio = NULL;
	}
	helios.clear();

	// Receiver will be deleted later
}


void SolarScene::initSolarScene(string scn_path)
{
	fstream in(scn_path, ios_base::in);
	if (in.fail()) {
		cerr << "[ERROR] Can't open scn file!!!" << endl;
		return;
	}

	string line, word;
	stringstream line_stream;
	InputMode input_mode = Initial;
	int row = 0;
	int col = 0;
	Layout* layout;
	int layout_type;
	Vector2d helio_gap;
	Vector2i helio_matrix;
	while (getline(in, line)) {
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
			recv->initRecv(in, input_mode);
			recvs.push_back(recv);
			break;
		}
		case HeliostatMode: {
			if (word == "gap")
				line_stream >> helio_gap.x() >> helio_gap.y();
			else if (word == "matrix")
				line_stream >> helio_matrix.x() >> helio_matrix.y();
			else {
				Heliostat* helio = new Heliostat(RectangularHelioType);
				helio->helio_index = helios.size();
				helio->helio_gap = helio_gap;
				helio->helio_matrix = helio_matrix;
				line_stream >> helio->helio_pos.x() >> helio->helio_pos.y() >> helio->helio_pos.z();

				getline(in, line);
				line_stream.clear();
				line_stream.str(line);
				line_stream >> helio->helio_size.x() >> helio->helio_size.y() >> helio->helio_size.z();
				helio->initFluxParam(recvs);
				helios.push_back(helio);
				helio = nullptr;
			}
			break;
		}
		case LayoutMode: {
			line_stream >> layout_type;
			LayoutCreator layout_creator;
			layout = layout_creator.getLayout((LayoutType)layout_type);
			layout->initLayout(in, input_mode);
			layouts.push_back(layout);
			break;
		}
		case Initial:
			break;
		default:
			break;
		}
	}

	in.close();
}
