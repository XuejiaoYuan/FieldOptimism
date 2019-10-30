//
// Created by Amber on 2018/4/3.
//

#include "Layout.h"


void Layout::initLayoutParams() {
	layout_first_helio_center = layout_bound_pos + 0.5 * helio_interval;
	layout_row_col.x() = layout_size.z() / helio_interval.z();		//row
	layout_row_col.y() = layout_size.x() / helio_interval.x();		//col
	helio_layout.clear();
	helio_layout.resize(layout_row_col.x(), vector<vector<Heliostat*>>(layout_row_col.y()));
}

void Layout::loadFieldArgs(json& field_args, double& z_start, int& rows, int& cols) {
	z_start = field_args["z_start"].as_double();
	rows = field_args["rows"].as<int>();
	cols = field_args["cols"].as<int>();
	helio_interval = Vector3d(field_args["interval"]["x"].as_double(), field_args["interval"]["y"].as_double(), field_args["interval"]["z"].as_double());
	layout_bound_pos = Vector3d(-(cols / 2.0) *  helio_interval.x(), 1, z_start - (rows + 0.5)*helio_interval.z());
	layout_size = Vector3d(helio_interval.x()*cols, helio_interval.y(), helio_interval.z()*rows);
}

void Layout::storeHelioToLayoutCore(Heliostat* helio) {
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

void Layout::storeHelioToLayout(vector<Heliostat*>& helios) {
	for (int i = 0; i < helios.size(); i++)
		storeHelioToLayoutCore(helios[i]);
}

void Layout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	// 1. Load field arguments
	double z_start;
	int rows, cols;
	loadFieldArgs(field_args, z_start, rows, cols);

	// 2. Initialize layout parameters
	initLayoutParams();

	// 3. Put heliostat  into layout
	int cnt = 0;
	vector<int> row_col(2, 0);
	Heliostat *helio;
	Heliostat h_tmp = argumentParser.getHelio();
	for (int i = 0; i < rows; i++) {
		z_start -= helio_interval.z();
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
}


void CrossRectLayout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	// 1. Load field arguments
	double z_start;
	int rows, cols;
	loadFieldArgs(field_args, z_start, rows, cols);

	// 2. Initialize layout parameters
	initLayoutParams();

	// 3. Put heliostat  into layout
	Heliostat* helio;
	Heliostat h_tmp = argumentParser.getHelio();
	for (int i = 0; i < rows; i++) {
		int tmp_col = cols - i % 2;
		z_start -= helio_interval.z();
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
}


void FermatLayout::createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios) {
	vector<double> recv_dis, angle_delta;
	vector<int> n_rows;
	Heliostat h_tmp  = argumentParser.getHelio();
	double dm = sqrt(pow(h_tmp.helio_size.x(), 2) + pow(h_tmp.helio_size.y(), 2)) + dsep;	// 定日镜对角线长度
	helio_interval = Vector3d(dm, dm, dm);
	dsep = field_args["dsep"].as_double();
	real_helio_num = argumentParser.getNumOfHelio();

	calcCircleParams(recv_dis, n_rows, angle_delta, field_args, dm);

	layout_bound_pos = Vector3d(-recv_dis.back(), 1, -recv_dis.back());
	layout_size = Vector3d(2 * recv_dis.back(), helio_interval.y(), 2 * recv_dis.back());
	initLayoutParams();

	for (int i = 0; i < helio_gap.size(); ++i)
		setCircleHelios(h_tmp, recv_dis[i], helio_gap[i], n_rows[i], angle_delta[i], helios, argumentParser.getReceivers());
}


void FermatLayout::setCircleHelios(Heliostat& h_tmp, const double R, const double gap, const int rows,
	const double angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs)
{
	Heliostat* helio;
	int h_cnt = 2 * PI / angle_delta;	
	
	for (int i = 0; i < rows; i++) {
		double angle_d = 2 * PI / h_cnt;
		double start_angle = (i + 1) % 2 * (angle_d / 2);
		double start_r = R + i * gap;

		for (int h = 0; h < h_cnt; h++) {
			helio = new Heliostat(h_tmp);
			helio->helio_index = helios.size();
			helio->helio_pos = Vector3d(sin(start_angle + h*angle_d) * start_r, h_tmp.helio_pos.y(), cos(start_angle + h*angle_d) * start_r);
			helio->initFluxParam(recvs);
			helios.push_back(helio);
			helio = nullptr;
		}
	}
}

void FermatLayout::calcCircleParams(vector<double>& recv_dis, vector<int>& n_rows, vector<double>& angle_delta, json& field_args, double dm) {
	recv_dis.push_back(field_args["helio_recv_dis"].as_double());
	int i = 1;
	int h_cnt = 0;
	auto iter = field_args["gap"].object_range().begin(), end = field_args["gap"].object_range().end();
	while (iter != end || h_cnt < real_helio_num) {
		if (h_cnt < real_helio_num) {
			field_args["gap"].push_back(helio_gap.back());
			helio_gap.push_back(helio_gap.back());
		}
		else
			helio_gap.push_back(iter->value().as<double>());
		if (angle_delta.empty())
			angle_delta.push_back(dm / recv_dis.front());
		else
			angle_delta.push_back(angle_delta.back() / pow(2, i - 1));
		n_rows.push_back(int(recv_dis[i] - recv_dis[i - 1]) / helio_gap[i - 1]);
		h_cnt += 2 * PI / angle_delta.back();
		++i;
		++iter;
	}
}

