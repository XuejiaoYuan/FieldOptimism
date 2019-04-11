//
// Created by Amber on 2018/4/3.
//

#include "Layout.h"

void Layout::initLayout(fstream& inFile, InputMode& input_mode, int& helio_type) {
	stringstream line_stream;
	string line, word;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> layout_bound_pos.x() >> layout_bound_pos.y() >> layout_bound_pos.z();
		else if (word == "size")
			line_stream >> layout_size.x() >> layout_size.y() >> layout_size.z();
		else if (word == "inter")
			line_stream >> helio_interval.x() >> helio_interval.y() >> helio_interval.z();
		else if (word == "n")
			line_stream >> helio_num;
		else
			line_stream >> helio_type;
	}
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
}

inline void Layout::setHelioLayout(Heliostat* helio)
{
	vector<int> row_col(2, 0);
	set<vector<int>> pos;
	for (auto&v : helio->vertex) {
		row_col[0] = (v.z() - layout_bound_pos.z()) / helio_interval.z();	// smaller z is, smaller row is
		row_col[1] = (v.x() - layout_bound_pos.x()) / helio_interval.x();	// smaller x is, smaller col is
		if (pos.count(row_col) == 0) {
			pos.insert(row_col);
			helio_layout[row_col[0]][row_col[1]].push_back(helio);
		}
	}
	
}

void Layout::setHelioLayout(vector<Heliostat*>& helios) 
{
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.clear();
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));

//#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++)
		setHelioLayout(helios[i]);
}

void Layout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	// 构建矩形镜场排布
	// 1. 调整Layout的各个参数
	helio_interval = Vector3d((*field_args[0])[0], (*field_args[0])[1], (*field_args[0])[2]);			// 定日镜间隔
	double z_start = (*field_args[1])[0];					// 定日镜第一行距离与接收器之间距离
	int row = int((*field_args[2])[0]);						// 镜场行数
	int col = int((*field_args[3])[0]);						// 镜场列数
	layout_bound_pos = Vector3d(
		-(col / 2.0) * helio_interval.x(),
		1,
		z_start - (row + 0.5)*helio_interval.z()
	);
	layout_size = Vector3d(
		helio_interval.x()*col,
		helio_interval.y(),
		helio_interval.z()*row
	);

	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));



	int cnt = 0;
	vector<int> row_col(2, 0);
	for (int i = 0; i < row; i++) {
		z_start -= helio_interval.z();
		double x_start = -(col / 2.0 - 0.5)*helio_interval.x();
		for (int j = 0; j < col; j++) {
			helio = new Heliostat((HelioType)helio_type, helio_gap, helio_matrix, helio_size, helios.size(), Vector3d(x_start, helio_pos.y(), z_start));
			helio->calcFluxParam(recvs);
			setHelioLayout(helio);
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}

}

void CrossRectLayout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	helio_interval = Vector3d((*field_args[0])[0], (*field_args[0])[1], (*field_args[0])[2]);		// 定日镜间隔
	cout << helio_interval.x() << ' ' << helio_interval.y() << ' ' << helio_interval.z() << endl;
	double z_start = (*field_args[1])[0];					// 定日镜第一行距离与接收器之间距离
	int row = int((*field_args[2])[0]);						// 镜场行数
	int col = int((*field_args[3])[0]);						// 镜场列数
	layout_bound_pos = Vector3d(
		-(col / 2.0) * helio_interval.x(),
		1,
		z_start - (row + 0.5)*helio_interval.z()
	);
	cout << layout_bound_pos.x() << ' ' << layout_bound_pos.y() << ' ' << layout_bound_pos.z() << endl;
	layout_size = Vector3d(
		helio_interval.x()*col,
		helio_interval.y(),
		helio_interval.z()*row
	);
	cout << layout_size.x() << ' ' << layout_size.y() << ' ' << layout_size.z() << endl;
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));


	int cnt = 0;
	for (int i = 0; i < row; i++) {
		int tmp_col = col - i % 2;
		z_start -= helio_interval.z();
		double x_start = -(tmp_col / 2.0 - 0.5)*helio_interval.x();
		for (int j = 0; j < tmp_col; j++) {
			helio = new Heliostat((HelioType)helio_type, helio_gap, helio_matrix, helio_size, helios.size(), Vector3d(x_start, helio_pos.y(), z_start));
			helio->calcFluxParam(recvs);
			setHelioLayout(helio);
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}

}

void FermatLayout::adjustHelioLayout(vector<Heliostat*>& helios, const vector<vector<double>*>& field_args, const vector<Receiver*>& recvs)
{
	vector<double> recv_dis, angle_delta;
	vector<int> n_rows;
	calcCircleParams(recv_dis, n_rows, angle_delta, field_args);
	// 调整layout的参数
	layout_bound_pos = Vector3d(
		-recv_dis.back(),
		1,
		-recv_dis.back()
	);
	layout_size = Vector3d(
		2 * recv_dis.back(),
		helio_interval.y(),
		2 * recv_dis.back()
	);
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));


	// 调整helio位置参数
	setCircleHelios(0, recv_dis[0], helio_gap1, n_rows[0], angle_delta[0], helios, recvs);
	setCircleHelios(1, recv_dis[1], helio_gap2, n_rows[1], angle_delta[1], helios, recvs);
	setCircleHelios(2, recv_dis[2], helio_gap3, n_rows[2], angle_delta[2], helios, recvs);


#ifdef DEBUG
	fstream outFile("fermat_helio.scn", ios_base::out);
	for (auto&h : helios) {
		outFile << h->helio_pos.x() << ' ' << h->helio_pos.z() << endl;
	}
	outFile.close();
#endif // DEBUG

}

void FermatLayout::setCircleHelios(const int field_index, const double R, const double gap, const int rows,
	const double angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	int cnt = 0;
	int h_cnt = 2 * PI / angle_delta;	
	
	for (int i = 0; i < rows; i++) {
		double angle_d = 2 * PI / h_cnt;
		double start_angle = (i + 1) % 2 * (angle_d / 2);
		double start_r = R + i * gap;

		for (int h = 0; h < h_cnt; h++) {
			helio = new Heliostat((HelioType)helio_type, helio_gap, helio_matrix, helio_size, helios.size(),
				Vector3d(
					sin(start_angle + h*angle_d) * start_r,
					helio_pos.y(),
					cos(start_angle + h*angle_d) * start_r
				),
				Vector3d(
					start_angle + h*angle_d,
					helio_pos.y(),
					start_r
				)
			);
			helio->calcFluxParam(recvs);
			setHelioLayout(helio);
			helios.push_back(helio);
			helio = nullptr;
			cnt++;
		}
	}
	cout << "cnt: " << cnt << endl;
}

void FermatLayout::calcCircleParams(vector<double>& recv_dis, vector<int>& n_rows, vector<double>& angle_delta, 
	const vector<vector<double>*>& field_args)
{
	double dm = sqrt(pow(helio_size.x(), 2) + pow(helio_size.y(), 2)) + dsep;	// 定日镜对角线长度

	if (!field_args.empty()) {
		dsep = (*field_args[0])[0];					// 定日镜包围盒安全距离
		helio_recv_dis1 = (*field_args[1])[0];		// 第一个同心圆与接收器之间的距离
		helio_gap1 = (*field_args[2])[0] * dm;				// 第一个同心环中定日镜分布间隔
		helio_gap2 = (*field_args[3])[0] * dm;				// 第二个同心环中定日镜分布间隔
		helio_gap3 = (*field_args[4])[0] * dm;				// 第三个同心换中定日镜分布间隔
		helio_interval = Vector3d(dm, dm, dm);				// 定日镜包围盒间距
	}
	recv_dis.push_back(helio_recv_dis1);
	recv_dis.push_back(2 * helio_recv_dis1);
	recv_dis.push_back(4 * helio_recv_dis1);
	recv_dis.push_back(8 * helio_recv_dis1);
	angle_delta.push_back(dm / helio_recv_dis1);
	angle_delta.push_back(angle_delta.back() / 2);
	angle_delta.push_back(angle_delta.back() / 4);
	n_rows.push_back(int(recv_dis[1] - recv_dis[0]) / helio_gap1);
	n_rows.push_back(int(recv_dis[2] - recv_dis[1]) / helio_gap2);
	n_rows.push_back(int(recv_dis[3] - recv_dis[2]) / helio_gap3);
}

MatrixXd FermatLayout::getCircleHelioIndex(int& start_index, const double R, const double gap, const int rows, 
	const double angle_delta)
{
	int h_cnt = 2 * PI / angle_delta;
	MatrixXd helio_index_store(rows, h_cnt);

	for (int i = 0; i < rows; i++) 
		for (int h = 0; h < h_cnt; h++) 
			(helio_index_store)(i, h) = start_index + i*h_cnt;
	start_index += rows * h_cnt;
}

vector<MatrixXd> FermatLayout::getHelioIndex() {
	vector<MatrixXd> helio_index_store;
	vector<double> recv_dis, angle_delta;
	vector<int> n_rows;
	int start_index = 0;
	calcCircleParams(recv_dis, n_rows, angle_delta);
	helio_index_store.push_back(getCircleHelioIndex(start_index, recv_dis[0], helio_gap1, n_rows[0], angle_delta[0]));
	helio_index_store.push_back(getCircleHelioIndex(start_index, recv_dis[1], helio_gap2, n_rows[1], angle_delta[1]));
	helio_index_store.push_back(getCircleHelioIndex(start_index, recv_dis[2], helio_gap3, n_rows[2], angle_delta[2]));

	return helio_index_store;
}