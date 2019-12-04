//
// Created by Amber on 2018/4/3.
//

#include "Layout.h"


void Layout::initLayoutParams() {
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.y() / helio_interval.y();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_interval.y() = layout_size.y() / layout_row_col.x();
	helio_interval.x() = layout_size.x() / layout_row_col.y();
	helio_layout.clear();
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));
}

void Layout::loadFieldArgs(ArgumentParser& argumentParser, json& field_args, double& z_start, int& rows, int& cols) {
	rows = argumentParser.getConfig()["Layout"]["rows"].as<int>();
	cols = argumentParser.getConfig()["Layout"]["cols"].as<int>();
	z_start = field_args["z_start"].as_double();
	double interval_ratio_x = field_args["interval_ratio_x"].as_double();
	double interval_ratio_y = field_args["interval_ratio_y"].as_double();
	Vector3d helio_size = argumentParser.getHelio().helio_size;
	double dm = sqrt(pow(helio_size.x(), 2) + pow(helio_size.z(), 2));
	helio_interval = Vector2d(dm*interval_ratio_x, dm*interval_ratio_y);
	layout_bound_pos = Vector2d(-(cols / 2.0) *  helio_interval.x(), z_start - (rows+0.5)*helio_interval.y());
	layout_size = Vector2d(helio_interval.x()*cols, helio_interval.y()*rows);
}

void Layout::storeHelioToLayoutCore(Heliostat* helio) {
	vector<int> row_col(2, 0);
	set<vector<int>> pos;
	for (auto&v : helio->vertex) {
		row_col[0] = (v.z() - layout_bound_pos.y()) / helio_interval.y();	// smaller z is, smaller row is
		row_col[1] = (v.x() - layout_bound_pos.x()) / helio_interval.x();	// smaller x is, smaller col is
		if (pos.count(row_col) == 0) {
			pos.insert(row_col);
			helio_layout[row_col[0]][row_col[1]].push_back(helio);
		}
	}
}

void Layout::storeHelioToLayout(vector<Heliostat*>& helios) {
	helio_layout.clear();
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));

	for (int i = 0; i < helios.size(); i++)
		storeHelioToLayoutCore(helios[i]);
}

void Layout::initLayout(fstream & in, InputMode & input_mode){
	stringstream line_stream;
	string line, word, _;
	int helio_num, helio_type;
	while (getline(in, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> layout_bound_pos.x() >> _ >> layout_bound_pos.y();
		else if (word == "size")
			line_stream >> layout_size.x() >> _ >> layout_size.y();
		else if (word == "inter")
			line_stream >> helio_interval.x() >> _ >> helio_interval.y();
		else if (word == "n")
			line_stream >> real_helio_num;
		else
			line_stream >> helio_type;
	}
	initLayoutParams();
}

void Layout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	// 1. Load field arguments
	double z_start;
	int rows, cols;
	loadFieldArgs(argumentParser, field_args, z_start, rows, cols);

	// 2. Initialize layout parameters
	initLayoutParams();

	// 3. Put heliostat  into layout
	int cnt = 0;
	vector<int> row_col(2, 0);
	Heliostat *helio;
	Heliostat h_tmp = argumentParser.getHelio();
	for (int i = 0; i < rows; i++) {
		z_start -= helio_interval.y();
		double x_start = -(cols / 2.0 - 0.5)*helio_interval.x();
		for (int j = 0; j < cols; j++) {
			helio = new Heliostat(h_tmp);
			helio->helio_index = helios.size();
			helio->helio_pos = Vector3d(x_start, h_tmp.helio_pos.y(), z_start);
			helio->initFluxParam(argumentParser.getReceivers());
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}
	real_helio_num = helios.size();
}


void CrossRectLayout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	// 1. Load field arguments
	double z_start;
	int rows, cols;
	loadFieldArgs(argumentParser, field_args, z_start, rows, cols);

	// 2. Initialize layout parameters
	initLayoutParams();

	// 3. Put heliostat  into layout
	Heliostat* helio;
	Heliostat h_tmp = argumentParser.getHelio();
	for (int i = 0; i < rows; i++) {
		int tmp_col = cols - i % 2;
		z_start -= helio_interval.y();
		double x_start = -(tmp_col / 2.0)*helio_interval.x();
		for (int j = 0; j < tmp_col; j++) {
			helio = new Heliostat(h_tmp);
			helio->helio_index = helios.size();
			helio->helio_pos = Vector3d(x_start, h_tmp.helio_pos.y(), z_start);
			helio->initFluxParam(argumentParser.getReceivers());
			helios.push_back(helio);
			helio = nullptr;
			x_start += helio_interval.x();
		}
	}
	real_helio_num = helios.size();
}


void RadialStaggerLayout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	vector<double> recv_dis;
	vector<int> n_rows, n_cols;
	Heliostat h_tmp  = argumentParser.getHelio();
	dsep = field_args["dsep"].as_double();
	double dm = sqrt(pow(h_tmp.helio_size.x(), 2) + pow(h_tmp.helio_size.z(), 2)) * dsep;	// 定日镜对角线长度
	helio_interval = Vector2d(dm, dm);
	real_helio_num = argumentParser.getNumOfHelio();
	calcCircleParams(recv_dis, n_rows, n_cols, field_args, dm);

	for (int i = 0; i < n_rows.size(); ++i) {
		bool status = setCircleHelios(h_tmp, i, recv_dis, n_rows, helio_gap[i], n_cols[i], helios, argumentParser.getReceivers());
		if (!status)
			break;
	}

	layout_bound_pos = Vector2d(-recv_dis.back() - dm / 2., -recv_dis.back() - dm / 2.);
	layout_size = Vector2d(2 * recv_dis.back() + dm, 2 * recv_dis.back() + dm);
	initLayoutParams();
}


bool RadialStaggerLayout::setCircleHelios(Heliostat& h_tmp, const int idx, vector<double>& recv_dis, vector<int>& rows, const double gap,
	const int col, vector<Heliostat*>& helios, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	double angle_delta = 2.*PI / col;

	for (int i = 0; i < rows[idx]; i++) {
		double start_angle = (i + 1) % 2 * (angle_delta / 2.0);
		double start_r = recv_dis[idx] + i * gap;

		for (int h = 0; h < col; h++) {
			helio = new Heliostat(h_tmp);
			helio->helio_index = helios.size();
			helio->helio_pos = Vector3d(sin(start_angle + h*angle_delta) * start_r, h_tmp.helio_pos.y(), cos(start_angle + h*angle_delta) * start_r);
			helio->initFluxParam(recvs);
			helios.push_back(helio);
			helio = nullptr;
		}

		if (helios.size() >= real_helio_num) {
			recv_dis.resize(idx + 1);
			recv_dis.push_back(recv_dis[idx] + (i+1) * gap);
			rows.back() = i + 1;
			return false;;
		}
	}
	return true;
}

void RadialStaggerLayout::calcCircleParams(vector<double>& recv_dis, vector<int>& n_rows, vector<int>& n_cols, json& field_args, double dm) {
	recv_dis.push_back(field_args["helio_recv_dis"].as_double());
	helio_gap = field_args["gap"].as<vector<double>>();
	int h_cnt = 0;
	int sz = helio_gap.size();
	for (int i = 0; i < sz || h_cnt < real_helio_num; ++i) {
		if (i>=sz && h_cnt < real_helio_num) {
			field_args["gap"].push_back(helio_gap.back()/dm);
			helio_gap.push_back(helio_gap.back());
		}
		else
			helio_gap[i] *= dm;
		if (n_cols.empty()) {
			n_cols.push_back(int(2 * PI / (dm / recv_dis.front())));
			if (n_cols.front() == 0)
				throw runtime_error("[INFO] error field parameters!!!");
		}
		else {
			n_cols.push_back(n_cols.back() * 2);
		}
		recv_dis.push_back(2 * recv_dis.back());
		n_rows.push_back(int(recv_dis[i + 1] - recv_dis[i]) / helio_gap[i]);
		h_cnt +=  n_cols.back() * n_rows.back();
	}
}

void SpiralLayout::createHelioAndLayout(ArgumentParser & argumentParser, json & field_args, vector<Heliostat*>& helios)
{
	double a = field_args["a"].as_double();
	double b = field_args["b"].as_double();
	int test_helio_num = field_args["test_helio_num"].as<int>();
	real_helio_num = argumentParser.getNumOfHelio();

	int i = 1;
	double eta = pow(2. / (1 + sqrt(5)), 2);
	vector<Receiver*>& recvs = argumentParser.getReceivers();
	vector<Vector3d>& fc_centers = argumentParser.getReceivers()[0]->getFocusCenter();
	vector<Vector3d>& recv_norms = argumentParser.getReceivers()[0]->getNormalList();
	Heliostat* helio;
	Heliostat h_tmp = argumentParser.getHelio();
	Vector2d min_pos(INT_MAX, INT_MAX);
	Vector2d max_pos(INT_MIN, INT_MIN);

	Vector3d helio_size = argumentParser.getHelio().helio_size;
	double dm = sqrt(pow(helio_size.x(), 2) + pow(helio_size.z(), 2));

	while (helios.size() < test_helio_num) {
		double r = a* pow(i, b);
		double theta = 2 * PI * eta * i;
		if (i > 1) {
			double front_r = a*pow(i - 1, b);
			double t = 2 * PI*eta;
			double dis = sqrt(front_r*front_r + r*r - 2 * front_r*r*cos(t));
			if (dis < dm) {
				++i;
				continue;
			}
		}
		double pos_x = r*cos(theta);
		double pos_z = r*sin(theta);
		bool status = false;
		for (int j = 0; j < fc_centers.size(); ++j) {
			Vector2d dir(pos_x - fc_centers[j].x(), pos_z - fc_centers[j].z());
			Vector2d n(recv_norms[j].x(), recv_norms[j].z());
			dir = dir.normalized();
			n = n.normalized();
			if (dir.dot(n) > Epsilon) {
				status = true;
				break;
			}
		}
		if (status) {
			helio = new Heliostat(h_tmp);
			helio->helio_index = helios.size();
			helio->helio_pos = Vector3d(pos_x, h_tmp.helio_pos.y(), pos_z);
			helio->initFluxParam(recvs);
			helios.push_back(helio);
			helio = nullptr;
			min_pos[0] = (min_pos[0] > pos_x ? pos_x : min_pos[0]);
			min_pos[1] = (min_pos[1] > pos_z ? pos_z : min_pos[1]);
			max_pos[0] = (max_pos[0] < pos_x ? pos_x : max_pos[0]);
			max_pos[1] = (max_pos[1] < pos_z ? pos_z : max_pos[1]);
		}
		++i;
	}

	helio_interval = Vector2d(dm, dm);
	layout_bound_pos = min_pos - helio_interval/2.;
	layout_size = max_pos - min_pos + helio_interval;

	initLayoutParams();
}
