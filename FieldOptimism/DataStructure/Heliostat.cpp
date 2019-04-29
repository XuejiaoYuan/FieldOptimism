//
// Created by Amber on 2018/4/3.
//
#include "Heliostat.h"

bool Heliostat::initSurfaceNormal(const vector<Vector3d> &focus_center, const Vector3d &sunray_dir) {
	double dis_min = INT_MAX;
	max_rela_dis = INT_MIN;
	min_rela_dis = INT_MAX;
	approx_rela_dis = INT_MIN;
	for (int i = 0; i < focus_center.size(); i++) {
		double dis = (focus_center[i] - helio_pos).norm();
		if (dis < dis_min) {
			dis_min = dis;
			focus_center_index = i;
		}
	}
	Vector3d reflectray_dir = focus_center[focus_center_index] - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();

	setHelioVertex();
	return true;
}

void Heliostat::changeSurfaceNormal(const vector<Vector3d>& focus_center, const Vector3d & sunray_dir)
{
	max_rela_dis = INT_MIN;
	min_rela_dis = INT_MAX;
	approx_rela_dis = INT_MIN;
	Vector3d reflectray_dir = focus_center[focus_center_index] - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();
	cos_w = (-sunray_dir).dot(helio_normal);

	setHelioVertex();
	calcFluxParam(focus_center[focus_center_index]);
}

void Heliostat::changeSubHelio(const Vector3d & focus_center, const Vector3d & sunray_dir)
{
	vector<Vector3d> root_dir(2);
	root_dir[0] = (vertex[1] - vertex[0]).normalized();
	root_dir[1] = (vertex[3] - vertex[0]).normalized();
	int subHelio_num = helio_matrix.x()*helio_matrix.y();
	int i = 0;
	for (auto&subHelio : subhelios) {
		if (helio_type == RectangularHelioType)
			subHelio->helio_normal = helio_normal;
		subHelio->setVertex(this, root_dir, i, focus_center, sunray_dir, true);
		i++;
	}

}

double Heliostat::calcSunHelioAngle(const Vector3d & sunray_dir)
{
	Vector3d reverse_sunray_dir = -sunray_dir;
	return reverse_sunray_dir.dot(this->helio_normal);
}

void Heliostat::initHeliostat(stringstream& line_stream, fstream& inFile, LayoutType layout_type, const Vector2d& helio_gap,
	const Vector2i& helio_matrix, const vector<Vector3d>& focus_center, const Vector3d& sunray_dir)
{
	string line;
	if (layout_type == FermatLayoutType) {
		line_stream >> helio_poly_pos.x() >> helio_poly_pos.y() >> helio_poly_pos.z();
		helio_pos.x() = cos(helio_poly_pos.x())*helio_poly_pos.z();
		helio_pos.z() = sin(helio_poly_pos.x())*helio_poly_pos.z();
		helio_pos.y() = helio_poly_pos.y();
	}
	else
		line_stream >> helio_pos.x() >> helio_pos.y() >> helio_pos.z();
	getline(inFile, line);
	line_stream.clear();
	line_stream.str(line);
	line_stream >> helio_size.x() >> helio_size.y() >> helio_size.z();
	this->helio_gap = helio_gap;
	this->helio_matrix = helio_matrix;
	if(sunray_dir != Vector3d(0,0,0))
		bool flag = initSurfaceNormal(focus_center, sunray_dir);
}

void Heliostat::initFluxParam(const vector<Receiver*>& recvs)
{
	double dis = set_focus_center_index(recvs);
	if (dis <= 1000)
		mAA = (double)(0.99321 - 0.0001176 * dis + 1.97 * 1e-8 * dis * dis);      //d<1000
	else
		mAA = exp(-0.0001106 * dis);

	S = helio_size.x() * helio_size.z();

	sigma_list.clear();
	sigma_list.resize(6);
	sigma_list[0] = dis;								// dis
	sigma_list[1] = SIGMA_SUN;							// sigma_sun
	sigma_list[2] = 2*SIGMA_S;					// sigma_bq
	sigma_list[3] = 0;									// sigma_ast
	sigma_list[4] = 0;									// sigma_t
	sigma_list[5] = abs(cos_phi[focus_center_index]);	// cos_rev
}

void Heliostat::calcFluxParam(const Vector3d& focus_center)
{
	Vector3d reverse_sunray_dir = (helio_pos - focus_center).normalized();

	// l_w_ratio 对总能量无影响
	vector<Vector3d> inter_v(3);
	for (int i = 0; i < 3; i++) {
		// 计算image plane顶点
		GeometryFunc::calcIntersection(reverse_sunray_dir, focus_center, vertex[i], -reverse_sunray_dir, inter_v[i]);
	}
	double ip_w = (inter_v[1] - inter_v[0]).norm();
	double ip_l = (inter_v[2] - inter_v[1]).norm();
	l_w_ratio = ip_l / ip_w;
	//l_w_ratio = 1;

	//sigma_list[3] = sqrt(S) * (1 - cos_w) / (4 * sigma_list[0]);
	//sigma = calcSigma();
	//sigma = 1.46;
	flux_param = 0.5 * S * cos_w * rou * l_w_ratio * mAA / PI;
}



void Heliostat::getSubHelioVertex(vector<Vector3d>& subhelio_vertex)
{
	if (helio_matrix.x() == 1 && helio_matrix.y() == 1)
		subhelio_vertex = vertex;

	for (auto&sub : subhelios) {
		for (auto&v : sub->vertex)
			subhelio_vertex.push_back(v);
	}
}


void Heliostat::setHelioVertex()
{
	GeometryFunc::setLocalVertex(helio_size.x(), helio_size.z(), vertex);

	GeometryFunc::getHelioMatrix(helio_normal, helio_pos, local2worldM, world2localM);

	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);

}

double Heliostat::set_focus_center_index(const vector<Receiver*>& recvs)
{
	int recv_index = 0;
	int ret_index = 0;
	double min_d = INT_MAX;
	for (int i = 0; i < recvs.size(); i++) {
		for (int j = 0; j < recvs[i]->focus_center.size(); j++) {
			Vector3d fc = recvs[i]->focus_center[j];
			double dis = sqrt(pow(fc.x() - helio_pos.x(), 2)
				+ pow(fc.y()-helio_pos.y(), 2)
				+ pow(fc.z() - helio_pos.z(), 2));
			if (dis < min_d) {
				recv_index = i;
				min_d = dis;
				ret_index = j;
			}
		}
	}
	focus_center_index = ret_index;

	Vector3d image_plane_normal = (helio_pos - recvs[recv_index]->focus_center[ret_index]).normalized();
	for (int i = 0; i < recvs.size(); i++) {
		for (int j = 0; j < recvs[i]->focus_center.size(); j++) {
			cos_phi.push_back(recvs[i]->recv_normal_list[j].dot(image_plane_normal));
		}
	}

	return min_d;
}

void Heliostat::calcLsfParam()
{
	lsf_param_M = MatrixXd(6, 6);
	lsf_param_v = MatrixXd(1, 6);

	double pos_x = helio_pos.x();
	double pos_y = helio_pos.z();
	lsf_param_v(0, 0) = 1;
	lsf_param_v(0, 1) = pos_x;
	lsf_param_v(0, 2) = pos_y;
	lsf_param_v(0, 3) = pos_x*pos_x;
	lsf_param_v(0, 4) = pos_x*pos_y;
	lsf_param_v(0, 5) = pos_y*pos_y;

	lsf_param_M.row(0) = lsf_param_v;
	lsf_param_M.row(1) = pos_x*lsf_param_v;
	lsf_param_M.row(2) = pos_y*lsf_param_v;
	lsf_param_M.row(3) = pos_x*pos_x*lsf_param_v;
	lsf_param_M.row(4) = pos_x*pos_y*lsf_param_v;
	lsf_param_M.row(5) = pos_y*pos_y*lsf_param_v;
}

double Heliostat::calcSigma()
{
	return sqrt(pow(sigma_list[0], 2) * (pow(sigma_list[1], 2) + pow(sigma_list[2], 2) + pow(sigma_list[3], 2) + pow(sigma_list[4], 2))) / sqrt(sigma_list[5]);
}

void Heliostat::initializeSubHelio(const Vector3d&focus_center, const Vector3d&sunray_dir)
{
	Vector2d subhelio_gap(0, 0);
	Vector2i subhelio_matrix(1, 1);

	int row = helio_matrix.x();
	int col = helio_matrix.y();
	Vector3d subhelio_size;
	subhelio_size.x() = (helio_size.x() - (col - 1)*helio_gap.x()) / col;
	subhelio_size.y() = helio_size.y();
	subhelio_size.z() = (helio_size.z() - (row - 1)*helio_gap.y()) / row;

	vector<Vector3d> root_dir(2);
	root_dir[0] = (vertex[1] - vertex[0]).normalized();
	root_dir[1] = (vertex[3] - vertex[0]).normalized();
	SubHelio* subHelio_ptr = nullptr;
	int subHelio_num = helio_matrix.x()*helio_matrix.y();
	for (int i = 0; i < subHelio_num; i++) {
		subHelio_ptr = new SubHelio();
		subHelio_ptr->helio_type = SubHelioType;
		subHelio_ptr->helio_gap = subhelio_gap;
		subHelio_ptr->helio_matrix = subhelio_matrix;
		subHelio_ptr->helio_size = subhelio_size;
		if (helio_type == RectangularHelioType)
			subHelio_ptr->helio_normal = helio_normal;
		subHelio_ptr->setVertex(this, root_dir, i, focus_center, sunray_dir);
		subhelios.push_back(subHelio_ptr);
		subHelio_ptr = nullptr;
	}
}

void SubHelio::setVertex(const Heliostat* root_helio, const vector<Vector3d>&root_dir, const int sub_index,
	const Vector3d&focus_center, const Vector3d&sunray_dir, const bool init)
{
	helio_index = sub_index;

	int row = root_helio->helio_matrix.x();
	int col = root_helio->helio_matrix.y();

	int current_row = sub_index / col;
	int current_col = sub_index % col;

	helio_pos = root_helio->vertex[0]
		+ current_row*(helio_size.z() + root_helio->helio_gap.y())*root_dir[0] + 0.5*helio_size.z()*root_dir[0]
		+ current_col*(helio_size.x() + root_helio->helio_gap.x())*root_dir[1] + 0.5*helio_size.x()*root_dir[1];

	if (root_helio->helio_type != RectangularHelioType) {
		Vector3d reflectray_dir = focus_center - helio_pos;
		reflectray_dir = reflectray_dir.normalized();
		helio_normal = (reflectray_dir - sunray_dir).normalized();
	}

	// TODO: consider sub helios
	//GeometryFunc::setWorldVertex();

}