//
// Created by Amber on 2018/4/3.
//
#include "Heliostat.h"


void Heliostat::changeSurfaceNormal(const Vector3d & sunray_dir, const ModelType& type, const bool& calc_sigma)
{
	Vector3d reflectray_dir = focus_center - helio_pos;
	reflectray_dir = reflectray_dir.normalized();
	helio_normal = (reflectray_dir - sunray_dir).normalized();
	cos_w = (-sunray_dir).dot(helio_normal);

	setHelioVertex();
	calcFluxParam(type, calc_sigma);
}


void Heliostat::initHelio(json& config) {
	helio_size = Vector3d(config["size"]["x"].as_double(), config["size"]["y"].as_double(), config["size"]["z"].as_double());
	helio_gap = Vector2d(0, 0);
	helio_matrix = Vector2i(1, 1);
	helio_pos = Vector3d(0, config["height"].as_double(), 0);
}

void Heliostat::initFluxParam(const vector<Receiver*>& recvs)
{
	vector<double> res = setFocusCenterIndex(recvs);
	double dis = res[0];

	if (dis <= 1000)
		mAA = (double)(0.99321 - 0.0001176 * dis + 1.97 * 1e-8 * dis * dis);      //d<1000
	else
		mAA = exp(-0.0001106 * dis);

	S = helio_size.x() * helio_size.z();

	sigma_list.clear();
	sigma_list.resize(6);
	sigma_list[0] = dis;											// dis
	sigma_list[1] = SIGMA_SUN;										// sigma_sun
	//sigma_list[2] = 2*SIGMA_S;									// sigma_bq
	//sigma_list[2] = sqrt(2 * (1 + pow(cos_w,2)))*SIGMA_S;
	//sigma_list[3] = sqrt(S) * (1 - cos_w) / (4 * sigma_list[0]);	// sigma_ast
	sigma_list[4] = 0;												// sigma_t
	sigma_list[5] = abs(res[1]);									// cos_rev
}

void Heliostat::calcFluxParam(const ModelType& type, const bool& calc_sigma)
{
	if (type == bHFLCAL) {
		// 测试改进模型时
		// 1. l_w_ratio对模型分布有影响
		// 2. rotate_theta对模型转向有影响
		Vector3d reverse_sunray_dir = (helio_pos - focus_center).normalized();
		vector<Vector3d> inter_v(3);
		Vector3d row_dir = (vertex[1] - vertex[0]).normalized();
		Vector3d col_dir = (vertex[2] - vertex[1]).normalized();

		double t1 = GeometryFunc::calcIntersection(reverse_sunray_dir, focus_center, helio_pos, -reverse_sunray_dir, inter_v[0]);
		double t2 = GeometryFunc::calcIntersection(reverse_sunray_dir, focus_center, helio_pos + row_dir, -reverse_sunray_dir, inter_v[1]);
		double t3 = GeometryFunc::calcIntersection(reverse_sunray_dir, focus_center, helio_pos + col_dir, -reverse_sunray_dir, inter_v[2]);

		Vector3d img_x_axis = reverse_sunray_dir.cross(Vector3d(0, 1, 0)).normalized();
		Vector3d img_y_axis = img_x_axis.cross(reverse_sunray_dir).normalized();

		col_dir = (inter_v[2] - inter_v[0]).normalized();
		row_dir = (inter_v[1] - inter_v[0]).normalized();
		double cos_eta_col = col_dir.dot(img_x_axis);
		double cos_eta_row = row_dir.dot(img_y_axis);
		rotate_theta = acos(abs(cos_eta_col));
		//double distortion_theta = acos(abs(cos_eta_row));
		double distortion_theta = 0;

		if (col_dir.dot(img_y_axis) < 0) rotate_theta = -rotate_theta;
		if (row_dir.dot(img_x_axis) < 0) distortion_theta = -distortion_theta;

		rotate_theta += distortion_theta;

		double ip_w = (inter_v[1] - inter_v[0]).norm();
		double ip_l = (inter_v[2] - inter_v[0]).norm();
	
		double w_l = ip_w / ip_l;
		l_w_ratio = 1 +  abs(log10(2 - 0.5*w_l));
		//l_w_ratio = 1 + abs(log10(w_l));

		//l_w_ratio = 1 + log10(w_l);
		//l_w_ratio = 1;

	}
	else if (type == iHFLCAL) {
		rotate_theta = 0;

		vector<Vector3d> inter_v(3);
		Vector3d reverse_sunray_dir = (helio_pos - focus_center).normalized();
		for(int i=0; i<3; ++i)
			double t = GeometryFunc::calcIntersection(reverse_sunray_dir, focus_center, vertex[i], -reverse_sunray_dir, inter_v[i]);
		double ip_w = (inter_v[1] - inter_v[0]).norm();
		double ip_l = (inter_v[2] - inter_v[0]).norm();

		double w_l = ip_l / ip_w;
		l_w_ratio = 1 + log10(w_l);
	}
	else {
		// 计算全镜场能量时：
		// 1. l_w_ratio对总能量不产生影响
		// 2. rotate_theta对总能量不产生影响
		// 3. 使用公式计算sigma代替峰值拟合
		l_w_ratio = 1;
		rotate_theta = 0;
	}
	if (calc_sigma) {
		sigma_list[2] = sqrt(2 * (1 + pow(cos_w, 2)))*SIGMA_S;
		sigma_list[3] = sqrt(S) * (1 - cos_w) / (4 * sigma_list[0]);	// sigma_ast
		sigma = calcSigma();
	}
	flux_param = 0.5 * S * cos_w * HELIOSTAT_REFLECTIVITY * l_w_ratio * mAA / PI;
}


void Heliostat::setHelioVertex()
{
	GeometryFunc::setLocalVertex(helio_size, vertex);

	GeometryFunc::getHelioMatrix(helio_normal, helio_pos, local2worldM, world2localM);

	vertex[0] = GeometryFunc::mulMatrix(vertex[0], local2worldM);
	vertex[1] = GeometryFunc::mulMatrix(vertex[1], local2worldM);
	vertex[2] = GeometryFunc::mulMatrix(vertex[2], local2worldM);
	vertex[3] = GeometryFunc::mulMatrix(vertex[3], local2worldM);

}

vector<double> Heliostat::setFocusCenterIndex(const vector<Receiver*>& recvs)
{
	// 暂时只处理一个接收器

	if (recvs.size() < 1) {
		throw std::runtime_error("[Error] Field without receivers!!!\n");
	}
	else if (recvs.size() > 1)
		throw std::runtime_error("[Error Heliostat] Can't deal with more than one receiver!!!\n");

	int recv_index = 0;
	int ret_index = 0;
	double min_d = INT_MAX;
	double ret_cos = 1;
	switch (recvs[0]->recv_type)
	{
	case RectangularRecvType:
	case PolyhedronRecvType: {
		for (int i = 0; i < recvs.size(); i++) {
			vector<Vector3d> fc_centers = recvs[i]->getFocusCenter();
			for (int j = 0; j < fc_centers.size(); j++) {
				Vector3d fc = fc_centers[j];
				double dis = sqrt(pow(fc.x() - helio_pos.x(), 2)
					+ pow(fc.y() - helio_pos.y(), 2)
					+ pow(fc.z() - helio_pos.z(), 2));
				if (dis < min_d) {
					recv_index = i;
					min_d = dis;
					ret_index = j;
				}
			}
		}
		focus_center_index = ret_index;
		focus_center = recvs[recv_index]->getFocusCenter()[ret_index];
		Vector3d image_plane_normal = (helio_pos - focus_center).normalized();
		for (int i = 0; i < recvs.size(); i++) {
			vector<Vector3d> fc_centers = recvs[i]->getFocusCenter();
			for (int j = 0; j < fc_centers.size(); j++) {
				cos_phi.push_back(recvs[i]->getNormalList()[j].dot(image_plane_normal));
			}
		}
		ret_cos = cos_phi[ret_index];
		break;
	}
	case CylinderRecvType: {
		focus_center_index = 0;		// 每个定日镜对应一个唯一聚焦点
		CylinderRecv* recv = dynamic_cast<CylinderRecv*>(recvs[0]);
		focus_center = recv->getFocusCenter(helio_pos);
		Vector3d image_plane_normal = (helio_pos - focus_center).normalized();
		cos_phi.push_back((focus_center - recv->recv_pos).normalized().dot(image_plane_normal));
		min_d = (helio_pos - focus_center).norm();
		ret_cos = cos_phi.back();
		break;
	}
	default:
		throw std::runtime_error("[Error Heliostat] Wrong receiver type!!!\n");
		break;
	}

	return { min_d, ret_cos };
}


double Heliostat::calcSigma()
{
	return sqrt(pow(sigma_list[0], 2) * (pow(sigma_list[1], 2) + pow(sigma_list[2], 2) + pow(sigma_list[3], 2) + pow(sigma_list[4], 2))) / sqrt(sigma_list[5]);
}

