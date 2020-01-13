#include "Receiver.h"

void Receiver::initRecv(json& config) {
	readRecvFromJson(config);
	initRecvCore();
}

void Receiver::initRecv(fstream & in, InputMode& input_mode)
{
	readRecvFromScn(in, input_mode);
	initRecvCore();
}

void Receiver::initRecvCore() {
	focus_center.push_back(recv_pos + Vector3d(recv_normal.array() * recv_size.array() / 2));
	recv_normal_list.push_back(recv_normal);
	vector<Vector3d> vertex = getRecvVertexCore(focus_center[0], recv_size.y() / 2.0, recv_size.x() / 2.0, recv_normal);
	recv_vertex.push_back(vertex);
}


void Receiver::readRecvFromJson(json& config) {
	recv_face_num = config["faceNum"].as<int>();
	recv_normal = Vector3d(config["norm"]["x"].as_double(), config["norm"]["y"].as_double(), config["norm"]["z"].as_double());
	recv_size = Vector3d(config["size"]["x"].as_double(), config["size"]["y"].as_double(), config["size"]["z"].as_double());
	recv_pos = Vector3d(config["pos"]["x"].as_double(), config["pos"]["y"].as_double(), config["pos"]["z"].as_double());
	if (config.contains("receiver_pixel_length")) {
		double recv_pixel_length = config["receiver_pixel_length"].as<double>();
		rows_cols = Vector2i(recv_size.y() / recv_pixel_length + 0.5, recv_size.x() / recv_pixel_length + 0.5);
	}
}

void Receiver::readRecvFromScn(fstream& in, InputMode& input_mode) {
	string line, word;
	stringstream line_stream;
	while (getline(in, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;
		if (word == "end") {
			input_mode = Initial;
			break;
		}
		else if (word == "pos")
			line_stream >> recv_pos.x() >> recv_pos.y() >> recv_pos.z();
		else if (word == "size")				// 接收器边长，高度，厚度
			line_stream >> recv_size.x() >> recv_size.y() >> recv_size.z();
		else if (word == "num")
			line_stream >> recv_face_num;
		else if (word == "norm")
			line_stream >> recv_normal.x() >> recv_normal.y() >> recv_normal.z();
		else
			line_stream >> recv_face;
	}
}

vector<Vector3d> Receiver::getRecvVertexCore(Vector3d & center, double half_l, double half_w, Vector3d & recv_normal)
{
	Vector3d down_cor(0, -1, 0);
	Vector3d cor_dir = recv_normal.cross(down_cor).normalized();
	vector<Vector3d> vertex = {
		(center - down_cor* half_l - cor_dir*half_w),				// 0 - 3
		(center + down_cor* half_l - cor_dir*half_w),				// |   |	
		(center + down_cor* half_l + cor_dir*half_w),				// 1 - 2
		(center - down_cor* half_l + cor_dir*half_w),
	};
	return vertex;
}


void PolyhedronRecv::initRecv(json& config) {
	readRecvFromJson(config);
	initRecvCore();
}

//void PolyhedronRecv::initRecv(fstream& in, InputMode& input_mode) {
//
//}

void PolyhedronRecv::initRecvCore() {
	Matrix3d m;
	double delta_angle = 2 * PI / recv_face_num;
	double pos = recv_size.y() / 2 / tan(delta_angle / 2);
	double half_l = recv_size.y() / 2.0;
	double half_w = recv_size.x() / 2.0;

	for (int i = 0; i < recv_face_num; i++) {
		m << cos(i*delta_angle), 0, sin(i*delta_angle),
			0, 1, 0,
			-sin(i*delta_angle), 0, cos(i*delta_angle);
		recv_normal_list.push_back(m * recv_normal);
		focus_center.push_back(recv_pos + recv_normal_list[i] * pos);

		vector<Vector3d> vertex = getRecvVertexCore(focus_center[i], half_l, half_w, recv_normal_list[i]);
		recv_vertex.push_back(vertex);
	}
}


void CylinderRecv::initRecv(json& config) {
	readRecvFromJson(config);
}

void CylinderRecv::initRecv(fstream& in, InputMode& input_mode) {
	readRecvFromScn(in, input_mode);
}

Vector3d CylinderRecv::getFocusCenter(Vector3d& helio_pos)
{
	Vector3d dir = (helio_pos - recv_pos).normalized();
	double r = recv_size.x() / (Vector2d(dir.x(), dir.z()).norm());
	Vector3d fc_center = Vector3d(recv_pos.x() + dir.x()*r, recv_pos.y(), recv_pos.z() + dir.z()*r);
	return fc_center;
}

vector<vector<Vector3d>> CylinderRecv::getRecvVertex(Vector3d& focus_center)
{
	if ((focus_center - Vector3d(0, 0, 0)).norm() < Epsilon)
		throw runtime_error("[Error CylinderRecv] Get receiver vertex without focus center!!!\n");
	vector<vector<Vector3d>> recv_vertex;
	Vector3d recv_normal = (focus_center - recv_pos).normalized();
	vector<Vector3d> vertex = getRecvVertexCore(recv_pos, recv_size.y() / 2.0, recv_size.x(), recv_normal);
	recv_vertex.push_back(vertex);
	return recv_vertex;
}

void CylinderRecv::readRecvFromJson(json & config)
{
	recv_size = Vector3d(config["size"]["radius"].as_double(), config["size"]["height"].as_double(), 0);
	recv_pos = Vector3d(config["pos"]["x"].as_double(), config["pos"]["y"].as_double(), config["pos"]["z"].as_double());
	if (config.contains("receiver_pixel_length")) {
		double recv_pixel_length = config["receiver_pixel_length"].as<double>();
		rows_cols = Vector2i(recv_size.y() / recv_pixel_length + 0.5, PI*recv_size.x() / recv_pixel_length + 0.5);
	}
}
