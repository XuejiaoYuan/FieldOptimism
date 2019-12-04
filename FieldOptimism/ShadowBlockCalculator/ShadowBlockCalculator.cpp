#include"ShadowBlockCalculator.h"
#include <random>
#include <ctime>
#include "../Tool/Timer/Timer.h"
#include <direct.h>


//
// [多边形裁剪] 处理阴影和遮挡，考虑阴影遮挡重叠部分
//
double SdBkCalc::helioClipper(Heliostat * helio, const vector<Vector3d>& dir, const vector<unordered_set<int>>& estimate_grids)
{
	vector<Heliostat*>& helios = solar_scene->helios;
	vector<Vector3d> helio_v, tmp_v(4);
	vector<float2> local_v;
	double t;
	Paths subj(1), clips;
	helio_v = helio->vertex;

	for (int i = 0; i < helio_v.size(); i++) {
		Vector3d tmp_v = GeometryFunc::mulMatrix(helio_v[i], helio->world2localM);
		local_v.push_back(make_float2(tmp_v.x(), tmp_v.z()));
		subj[0] << IntPoint(VERTEXSCALE*local_v[i].x, VERTEXSCALE*local_v[i].y);
	}

	for (int index = 0; index < 2; index++) {
		Vector3d reverse_dir = -dir[index];
		vector<Vector3d> pro(4), tmp_pro(4);
		set<Heliostat*> relative_helio_set;
		for (auto& rela_index: estimate_grids[index]) {
			Heliostat* relative_helio = helios[rela_index];
			helio_v = relative_helio->vertex;
			int cnt = 0;
			for (int i = 0; i < helio_v.size(); i++) {
				t = GeometryFunc::calcIntersection(helio->helio_normal, helio->helio_pos, helio_v[i], reverse_dir, pro[i]);
				if (t > Epsilon) 
					cnt++;
			}
			if (cnt > 0) {
				Path clip;
				for (auto v : pro) {
					v = GeometryFunc::mulMatrix(v, helio->world2localM);
					clip << IntPoint(VERTEXSCALE *v.x(), VERTEXSCALE * v.z());
				}
				clips.push_back(clip);
			}
		}
	}

	Clipper c;
	Paths solution;													// solution represents the shadowing / blocking area
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clips, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);

	// 计算剩余面积
	c.Clear();
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(solution, ptClip, true);
	c.Execute(ctDifference, solution, pftNonZero, pftNonZero);

	Vector3d h_local_centerBias;
	double area = PolygonCenterCalculator::calcPolygonCenter(solution, h_local_centerBias, calcCenterMode);
	if (calcCenterMode) {
		h_local_centerBias = GeometryFunc::mulMatrix(h_local_centerBias, helio->local2worldM);
		helio->centerBias = make_float3(h_local_centerBias.x(), h_local_centerBias.y(), h_local_centerBias.z());
	}
	double res = area / helio->S;

#ifdef CLIPPER	
	unordered_set<int> test_helio_index = { 75, 1699, 2600, 4000, 5297, 6000, 7252 };
	Receiver* recv = solar_scene->recvs[0];
	Matrix4d rlocal2worldM, rworld2localM;
	GeometryFunc::getHelioMatrix(recv->get_normal_list()[0], recv->get_focus_center()[0], rlocal2worldM, rworld2localM);

	if (test_helio_index.count(helio->helio_index)) {	
		fstream outFile("Outputfiles/rectField/Clipper/" + to_string(helio->helio_index) + ".txt", ios_base::out | ios_base::app);
		for (int i = 0; i < solution.size(); i++) {
			int n = solution[i].size();
			for (int j = 0; j < n; j++) {
				Vector3d v(solution[i][j].X / (double)VERTEXSCALE, 0, solution[i][j].Y / (double)VERTEXSCALE);
				v = GeometryFunc::mulMatrix(v, helio->local2worldM);
				GeometryFunc::calcIntersection(recv->get_normal_list()[0], recv->get_focus_center()[0], v, dir[1], v);
				v = GeometryFunc::mulMatrix(v, rworld2localM);
				outFile << v.x() << ' ' << v.z() << endl;
			}
		}
		GeometryFunc::calcIntersection(recv->get_normal_list()[0], recv->get_focus_center()[0], h_local_centerBias, dir[1], h_local_centerBias);
		h_local_centerBias = GeometryFunc::mulMatrix(h_local_centerBias, rworld2localM);
		outFile << h_local_centerBias.x() << ' ' << h_local_centerBias.z() << endl;
		outFile << "#" << endl;
		outFile.close();

		outFile.open("Outputfiles/rectField/Clipper/v_" + to_string(helio->helio_index) + ".txt", ios_base::out | ios_base::app);
		c.Clear();
		c.AddPaths(subj, ptSubject, true);
		clips.clear();
		c.AddPaths(clips, ptClip, true);
		c.Execute(ctDifference, solution, pftNonZero, pftNonZero);
		vector<vector<Vector3d>> v_list(solution.size());
		for (int i = 0; i < solution.size(); i++) {
			int n = solution[i].size();
			for (int j = 0; j < n; j++) {
				Vector3d v(solution[i][j].X / (double)VERTEXSCALE, 0, solution[i][j].Y / (double)VERTEXSCALE);
				v = GeometryFunc::mulMatrix(v, helio->local2worldM);
				GeometryFunc::calcIntersection(recv->get_normal_list()[0], recv->get_focus_center()[0], v, dir[1], v);
				v = GeometryFunc::mulMatrix(v, rworld2localM);
				outFile << v.x() << ' ' << v.z() << endl;
			}
		}
		outFile.close();
	}
#endif // CLIPPER

	return 1.0 - fabs(res);
}



//
//	[采样计算通量密度] 以接收器中心点通量密度代表区域通量平均通量密度
//		将接收器网格中心点
//	TODO: 当接收器为圆柱体时，由于圆柱体的grid与image plane夹角均不同，因此不能使用，需要进一步修正
float SdBkCalc::calcHelio2RecvEnergy(vector<Vector3d>& recv_v, Vector3d& recv_n, Vector2i& rows_cols, 
	Heliostat* helio, Vector3d& fc_center, double DNI, double cos_phi) {
	int rows = rows_cols.x();
	int cols = rows_cols.y();
	//int rows = 240;
	//int cols = 240;
	Vector3d row_gap = (recv_v[1] - recv_v[0]) / rows;
	Vector3d col_gap = (recv_v[3] - recv_v[0]) / cols;

	Vector3d reverse_dir = (helio->helio_pos - fc_center).normalized();
	Matrix4d world2local, local2world;
	GeometryFunc::getImgPlaneMatrixs(reverse_dir, fc_center, local2world, world2local, 1);

	double sum = 0;
	double grid_area = 0.0;

	Vector3d i_center_bias(0, 0, 0);
	if (calcCenterMode) {
		Vector3d h_center_bias(helio->centerBias.x, helio->centerBias.y, helio->centerBias.z);
		GeometryFunc::calcIntersection(reverse_dir, fc_center, h_center_bias, -reverse_dir, i_center_bias);
		i_center_bias = GeometryFunc::mulMatrix(i_center_bias, world2local);
	}

	double sigma = helio->sigma;
	Matrix4d hWorld2local, hLocal2World;
	ModelType model_type = solar_scene->getModelType();
	vector<Vector3d> proj_v(4);
	if (model_type == Gracia) {
		GeometryFunc::getHelioMatrix(helio->helio_normal, helio->helio_pos, hLocal2World, hWorld2local);
		for (int i = 0; i < 4; ++i) {
			GeometryFunc::calcIntersection(helio->helio_normal, helio->helio_pos, recv_v[i], reverse_dir, proj_v[i]);
			proj_v[i] = GeometryFunc::mulMatrix(proj_v[i], hWorld2local);
		}
		for (int i = 0; i < 4; ++i)
			grid_area += proj_v[i].x()*proj_v[(i + 1) % 4].z() - proj_v[i].z()*proj_v[(i + 1) % 4].x();
		grid_area = fabs(grid_area / 2.0 / rows / cols);
	}
	else if (model_type == HFLCAL) {
		GeometryFunc::getHelioMatrix(recv_n, fc_center, hLocal2World, hWorld2local);
		grid_area = (recv_v[1] - recv_v[0]).norm() * (recv_v[3] - recv_v[0]).norm() / rows / cols;
	}
	else {
		for (int i = 0; i < 4; ++i) {
			GeometryFunc::calcIntersection(reverse_dir, fc_center, recv_v[i], reverse_dir, proj_v[i]);
			proj_v[i] = GeometryFunc::mulMatrix(proj_v[i], world2local);
		}

		for (int i = 0; i < 4; ++i)
			grid_area += proj_v[i].x()*proj_v[(i + 1) % 4].z() - proj_v[i].z()*proj_v[(i + 1) % 4].x();
		grid_area = fabs(grid_area / 2.0 / rows / cols);
	}

	vector<string> model_name = {"m1", "m2", "m3", "m4"};
	_mkdir(output_path.c_str());
	fstream outFile(output_path + model_name[model_type] + '_' + to_string(helio->helio_index) + ".txt", ios_base::out);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			Vector3d start_v = recv_v[0] + i*row_gap + j*col_gap;
			Vector3d trans_v, inter_v;

			double res = 0;

			if(model_type == Gracia) {
				GeometryFunc::calcIntersection(helio->helio_normal, helio->helio_pos, start_v, -reverse_dir, inter_v);
				trans_v = GeometryFunc::mulMatrix(inter_v, hWorld2local);
				res = DNI*(1 - helio->sd_bk)*helio->mAA*helio->S*HELIOSTAT_REFLECTIVITY / 2. / PI*gl->flux_func(trans_v.x(), trans_v.z(), sigma, 1)*cos_phi;
			}
			else if (model_type == HFLCAL) {
				trans_v = GeometryFunc::mulMatrix(start_v, hWorld2local);
				res = DNI*helio->flux_param * gl->flux_func(trans_v.x(), trans_v.z(), sigma, 1);
			}
			else {
				GeometryFunc::calcIntersection(reverse_dir, fc_center, start_v, reverse_dir, inter_v);
				Vector3d proj_v = GeometryFunc::mulMatrix(inter_v, world2local) - i_center_bias;
				//helio->rotate_theta = 0;
				trans_v.x() = proj_v.x()*cos(helio->rotate_theta) + proj_v.z()*sin(helio->rotate_theta);
				trans_v.z() = proj_v.z()*cos(helio->rotate_theta) - proj_v.x()*sin(helio->rotate_theta);
				res = DNI*(1 - helio->sd_bk)*cos_phi*helio->flux_param * gl->flux_func(trans_v.x(), trans_v.z(), sigma, helio->l_w_ratio);
			}
			sum += res *grid_area /cos_phi;
			outFile << start_v.x() << ' ' << start_v.y() << ' ' << res << endl;
		}
	}

	outFile << i_center_bias.x() << ' ' << i_center_bias.z() << endl;
	outFile.close();
	return sum;
}


inline double calc_mAA(double dis) {
	double mAA;
	if (dis <= 1000)
		mAA = (double)(0.99321 - 0.0001176 * dis + 1.97 * 1e-8 * dis * dis);      //d<1000
	else
		mAA = exp(-0.0001106 * dis);
	return mAA;
}

 
inline double random_double() {
	static std::default_random_engine e(time(NULL));
	static std::uniform_real_distribution<double> u(0, 1);
	return u(e);
}


//
// [计算定日镜的阴影遮挡] 计算目标定日镜的阴影与遮挡结果
//
double SdBkCalc::calcHelioShadowBlock(int helio_index)
{
	Heliostat* helio = solar_scene->helios[helio_index];
	vector<unordered_set<int>> estimate_index(2);
	Vector3d reflect_dir = (helio->focus_center - helio->helio_pos).normalized();
	GridDDA dda_handler;
	dda_handler.rayCastForShadow(solar_scene, helio, -solar_scene->sunray_dir, estimate_index[0]);

	if (!block_grid_init[helio_index]) {
		dda_handler.rayCastForBlock(solar_scene, helio, rela_block_grid_index[helio_index]);
		block_grid_init[helio_index] = true;
	}
	dda_handler.getBlockHelioFromGrid(solar_scene, rela_block_grid_index[helio_index], estimate_index[1], helio);

	vector<Vector3d> dir = { -solar_scene->sunray_dir, reflect_dir };
	helio->sd_bk = helioClipper(helio, dir, estimate_index);

	return helio->sd_bk;
}

//
// [计算定日镜的通量密度和] 计算目标定日镜投射到接收器面上的能量和
//
void SdBkCalc::calcSceneFluxDistribution(vector<int>& test_helio_index, const double DNI, json& config) {
	gl = new GaussLegendreCPU(config["M"].as<int>(), config["N"].as<int>(), config["m"].as<int>(), config["n"].as<int>());
	Vector2i rows_cols = solar_scene->recvs[0]->rows_cols;
	double sum = 0;
	for (int h_index : test_helio_index) {

		Heliostat* helio = solar_scene->helios[h_index];
		//if (helio->sd_bk > Epsilon) {
		//	system("pause");
		//	cout << helio->sd_bk << endl;
		//}
		int fc_index = helio->focus_center_index;
		vector<Receiver*> recvs = solar_scene->recvs;
		Vector3d reverse_dir = (helio->helio_pos - helio->focus_center).normalized();		// The normal of image plane
		
		Matrix4d world2localM, local2worldM;
		GeometryFunc::getImgPlaneMatrixs(reverse_dir, helio->focus_center, local2worldM, world2localM, 1);
		fstream outFile("Outputfiles/RectField/flux_param_100.txt", ios_base::out | ios_base::app);
		outFile << h_index << endl;
		outFile << "0 0" << endl;
		for (int i = 0; i < 4; ++i) {
			Vector3d inter_v;
			GeometryFunc::calcIntersection(reverse_dir, helio->focus_center, helio->vertex[i], -reverse_dir, inter_v);
			inter_v = GeometryFunc::mulMatrix(inter_v, world2localM);
			outFile << inter_v.x() << ' ' << inter_v.z() << endl;
		}
		outFile.close();
		
		vector<vector<Vector3d>> recv_vertex = recvs[0]->getRecvVertex(helio->focus_center);
		vector<Vector3d> recv_normal = recvs[0]->getNormalList(helio->focus_center);
		double res = 0;
		for (int i = 0; i < helio->cos_phi.size(); i++) {
			if (helio->cos_phi[i] > Epsilon) {
				res += calcHelio2RecvEnergy(recv_vertex[i], recv_normal[i], rows_cols, helio, helio->focus_center, DNI, helio->cos_phi[i]);
			}
		}
		sum += res;
		cout << res << endl;
	}

	cout << sum << endl;
	delete gl;
}


void SdBkCalc::calcSceneShadowBlock()
{
	vector<Heliostat*> helios = solar_scene->helios;

#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++) {
		calcHelioShadowBlock(i);
	}
}


//
// [测试镜场所有定日镜阴影遮挡预处理函数]
//
void SdBkCalcTest::totalHeliosTest(const string& _save_path) {
	save_path = _save_path;
	readRayTracingRes();
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		singleHelioTest(i);
	}
}

//
// [测试单个定日镜阴影遮挡预处理函数]
//
void SdBkCalcTest::singleHelioTest(const int _helio_index) {
	int fc_index = solar_scene->helios[_helio_index]->focus_center_index;
	Vector3d reflect_dir = solar_scene->helios[_helio_index]->focus_center - solar_scene->helios[_helio_index]->helio_pos;
	setDir(-solar_scene->sunray_dir, reflect_dir.normalized());
	setTestIndex(_helio_index);

	//rayTracingSdBk();
	readRayTracingRes(_helio_index);
	//normalSdBk();
	boundingSphereSdBk();
	//neighRowSdBk();
	//improvedNeighSdBk();
	use3dddaSdBk();
}

//
// [ray tracing] 辅助函数，读取ray tracing文件结果
//
void SdBkCalcTest::readRayTracingRes(){
	readRayTracingCore("RC_sdbk_index.txt", v_gt_sd_helio_index, v_gt_bk_helio_index);
}

//
// [ray tracing] 輔助函數，讀取文件內容
//
void SdBkCalcTest::readRayTracingCore(string file_name, vector<unordered_set<int>>& sd_index_set, vector<unordered_set<int>>& bk_index) {
	fstream inFile(save_path + file_name);
	string line;
	while (getline(inFile, line)) {
		stringstream ss(line);
		string word;
		//getline(ss, line, ' ');
		unordered_set<int> tmp;
		while (ss >> word) {
			tmp.insert(stoi(word));
		}
		sd_index_set.push_back(tmp);

		getline(inFile, line);
		ss.clear();
		ss.str(line);
		//getline(ss, line, ' ');
		tmp.clear();
		while (ss >> word) {
			tmp.insert(stoi(word));
		}
		bk_index.push_back(tmp);
	}
	inFile.close();
}

//
// [ray tracing] 辅助函数，读取ray tracing文件结果
//
void SdBkCalcTest::readRayTracingRes(int index) {
	gt_sd_helio_index.clear();
	gt_bk_helio_index.clear();
	gt_sd_helio_index = v_gt_sd_helio_index[index];
	gt_bk_helio_index = v_gt_bk_helio_index[index];
}



// [ray tracing] 处理阴影或遮挡
//	使用ray tracing方式计算blocking和shadowing
//	将定日镜均匀剖分，每个小区域中心点 发出一根光线，与周围定日镜求交
//	记录产生blocking和shadowing的定日镜的编号
//	version: CPU
void SdBkCalcTest::rayTracingSdBk() {
	cout << "Ray Tracing" << endl;
	Timer::resetStart();
	gt_sd_helio_index.clear();
	gt_bk_helio_index.clear();

	int helio_inner_rows = 20;
	int helio_inner_cols = 20;

	double total_sum = helio_inner_rows *helio_inner_cols;
	int cnt = 0;

	auto helio = solar_scene->helios[helio_index];
	auto helio_v = helio->vertex;
	Vector3d row_dir = (helio_v[1] - helio_v[0]) / helio_inner_rows;
	Vector3d col_dir = (helio_v[3] - helio_v[0]) / helio_inner_cols;

	for (int i = 0; i <= helio_inner_rows; i++) {
		for (int j = 0; j <= helio_inner_cols; j++) {
			Vector3d ori_v = helio_v[0] + i*row_dir + j*col_dir;
			if (calcIntersect(ori_v, sd_dir, gt_sd_helio_index)) {
				++cnt;
			}
			else if (calcIntersect(ori_v, bk_dir, gt_bk_helio_index)) {
				++cnt;
			}
		}
	}

	double res = cnt / total_sum;

	double time = Timer::getDuration();
	fstream outFile(save_path + "/rayTracing_time.txt", ios_base::out | ios_base::app);
	outFile << time << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_sdbk.txt", ios_base::out | ios_base::app);
	outFile << helio_index<< ' ' << res << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_sd_index.txt", ios_base::out | ios_base::app);
	outFile << helio_index << ' ';
	for (auto&n : gt_sd_helio_index)
		outFile << n << ' ';
	outFile << endl;
	outFile.close();
	outFile.open(save_path + "/rayTracing_bk_index.txt", ios_base::out | ios_base::app);
	outFile << helio_index << ' ';
	for (auto&n : gt_bk_helio_index)
		outFile << n << ' ';
	outFile << endl;
	outFile.close();
}

//
// [ray tracing] 计算阴影或遮挡
//		返回当前起始点发射的光线是否与定日镜相交
bool SdBkCalcTest::calcIntersect(Vector3d& ori_v, Vector3d& dir, unordered_set<int>& index_set) {
	double tMin = INT_MAX;
	Heliostat* hNear = nullptr;

	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
			vector<Vector3d> neigh_v = h->vertex;
			Vector3d inter_v;
			double t = GeometryFunc::calcIntersection(h->helio_normal, h->helio_pos, ori_v, dir, inter_v);
			if (t < Epsilon) continue;
			if (GeometryFunc::inProjArea(neigh_v, inter_v)) {
				if (t < tMin) {
					tMin = t;
					hNear = h;
				}
			}
		}
	}
	if (hNear != nullptr) {
		index_set.insert(hNear->helio_index);
		return true;
	}
	return false;
}


//
// [法向剔除] 计算两个定日镜之间的法向点积并进行剔除
//		设置相关定日镜的最远距离，以减少相关定日镜数量
void SdBkCalcTest::normalSdBk() {
	Timer::resetStart();

	Heliostat* cur_h = solar_scene->helios[helio_index];
	unordered_set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
//#pragma omp parallel for
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
			if (h->helio_normal.dot(sd_dir) > Epsilon) {
				if (checkHelioDis(sd_dir, cur_h, cur_h->helio_pos, h->helio_pos)) {
					//#pragma omp critical
					if (gt_sd_helio_index.count(h->helio_index)) {
					//#pragma omp critical
						++sd_ac;
						
					}
					sd_set.insert(h->helio_index);
				}
			}
			if (h->helio_normal.dot(bk_dir) > Epsilon) {
				if (checkHelioDis(bk_dir, cur_h, cur_h->helio_pos, h->helio_pos)) {
					//#pragma omp critical
					if (gt_bk_helio_index.count(h->helio_index)) {
					//#pragma omp critical
						++bk_ac;
						
					}
					bk_set.insert(h->helio_index);
				}
			}
		}
	}

	Timer::printDuration("normalSdBk");
	saveTestRes("normal", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [法向剔除] 法向剔除輔助函數，根據距離排除無關定日鏡
//
bool SdBkCalcTest::checkHelioDis(Vector3d& dir, Heliostat* helio, Vector3d& Hloc, Vector3d& HIloc) {
	double tanpi2zen = dir.y() / sqrt(dir.x()*dir.x() + dir.z()*dir.z());
	Vector3d HIt = helio->helio_normal;
	double HIzen = acos(HIt.y());
	double HIh = helio->helio_size.x();
	double l_max = (HIloc.y() - Hloc.y() + HIh*sin(HIzen)) / tanpi2zen + HIh*HIt.y();
	double interaction_limit = 100;
	l_max = fmin(l_max, interaction_limit*HIh);	//limit to a reasonable number
	double hdist = (HIloc - Hloc).norm();
	if (hdist < l_max) return true;
	else return false;
}

//
// [bounding sphere] 通過計算兩個定日鏡包圍盒之間的距离判断是否相关
//
void SdBkCalcTest::boundingSphereSdBk() {
	Timer::resetStart();
	Heliostat* cur_h = solar_scene->helios[helio_index];
	unordered_set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
	for (int i = 0; i < solar_scene->helios.size(); ++i) {
		auto h = solar_scene->helios[i];
		if (h->helio_index != helio_index) {
			if(checkBoundingBox(cur_h->helio_pos, h->helio_pos, sd_dir, diameter)) {
				if (gt_sd_helio_index.count(h->helio_index)) {
					++sd_ac;
				}
				sd_set.insert(h->helio_index);
			}
			if (checkBoundingBox(cur_h->helio_pos, h->helio_pos, bk_dir, diameter)) {
				if (gt_bk_helio_index.count(h->helio_index)) {
					++bk_ac;
				}
				bk_set.insert(h->helio_index);

			}
		}
	}

	Timer::printDuration("Bounding Sphere");
	saveTestRes("boundingBox", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [bounding box] bounding box辅助函数，判断是否有相交的可能性
//
bool SdBkCalcTest::checkBoundingBox(Vector3d& Hloc, Vector3d& HIloc, Vector3d& dir, double diameter) {
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}


//
// [neighbor Row] 对当前定日镜仅考虑当前行和相邻两行的定日镜，由东西、南北方向判断是否有可能存在阴影遮挡
//		论文表述不够清晰，所以在使用东西南北方向判断时根据bounding box
void SdBkCalcTest::neighRowSdBk() {
	Timer::resetStart();
	Heliostat* cur_h = solar_scene->helios[helio_index];
	set<int> sdbk_set;
	int sdbk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
	int start = 0;
	int end = 0;
	getStartEndIndex(cur_h, start, end);
//#pragma omp parallel for
	for (; start <= end; ++start) {
		auto h = solar_scene->helios[start];
		if (start != helio_index) {
			if (checkEffectRegion(Vector3d( 1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
				checkEffectRegion(Vector3d(-1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
				checkEffectRegion(Vector3d( 0, 0, 1), cur_h->helio_pos, h->helio_pos, diameter) ||
				checkEffectRegion(Vector3d(0, 0, -1), cur_h->helio_pos, h->helio_pos, diameter)) {
				//#pragma omp critical
				
				if (gt_sd_helio_index.count(h->helio_index) || gt_bk_helio_index.count(h->helio_index)) {
					//#pragma omp critical
					++sdbk_ac;
				}
				sdbk_set.insert(h->helio_index);
			}
		}
	}
	Timer::printDuration("neighRowSdBk");
	fstream outFile(save_path + "/neighborRow_pr.txt", ios_base::app | ios_base::out);
	outFile << helio_index << ' ' << sdbk_ac << ' ' << sdbk_set.size() << endl;
	outFile.close();

	outFile.open(save_path + "/neighborRow_sdbk_index.txt", ios_base::app | ios_base::out);
	outFile << helio_index << ' ';
	for (auto&n : sdbk_set)
		outFile << n << ' ';
	outFile << endl;
	outFile.close();

}


//
// [neighbor row] neighbor row的辅助函数，用于确定当前定日镜当前及相邻两行的定日镜起始和终止坐标
//
void SdBkCalcTest::getStartEndIndex(Heliostat* helio, int& start, int& end) {
	vector<Heliostat*>& helios = solar_scene->helios;
	double cur_row = helio->helio_pos.z();

	switch (solar_scene->layouts[0]->layout_type)
	{
	case RectLayoutType:
	case CrossRectLayoutType: {
		start = helio->helio_index;
		double gap = 2 * solar_scene->layouts[0]->helio_interval.y();
		while (start >= 0 && helios[start]->helio_pos.z() <= cur_row + gap) --start;
		if (start != helio->helio_index)
			++start;
		end = helio->helio_index;
		while (end < helios.size() && helios[end]->helio_pos.z() >= cur_row - gap) ++end;
		if (end != helio->helio_index)
			--end;
		break;
	}
		
	case RadialStaggerLayoutType: {
		RadialStaggerLayout* layout = dynamic_cast<RadialStaggerLayout*>(solar_scene->layouts[0]);
		vector<MatrixXd> helio_index_store = layout->getHelioIndex();
		int cur_region = 0;
		int sum = 0;
		while (helio->helio_index > sum - 1) {
			int size = helio_index_store[cur_region].rows() *helio_index_store[cur_region].cols();
			sum += size;
			++cur_region;
		}
		sum -= helio_index_store[cur_region].rows() *helio_index_store[cur_region].cols();
		--cur_region;
		start = helio->helio_index - sum + 1;
		int row = start / helio_index_store[cur_region].cols();
		if (row < 2 && cur_region == 0) start = 0;
		else {
			row -= 2;
			if (row < 0) {
				--cur_region;
				row %= helio_index_store[cur_region].rows();
			}
			start = (helio_index_store[cur_region])(row, 0);
		}
		if (row > helio_index_store[cur_region].rows() - 3 && cur_region == helio_index_store.size() - 1) {
			end = helios.size() - 1;

		}
		else {
			row += 2;
			if (row > helio_index_store[cur_region].rows() - 1) {
				row %= helio_index_store[cur_region].rows();
				++cur_region;
			}
			int col = helio_index_store[cur_region].cols() - 1;
			end = (helio_index_store[cur_region])(row, col);
		}
		break;
	}
		
	case SpiralLayoutType:

		break;
	default:
		break;
	}
}


//
// [neighbor row] neighbor row的辅助函数，通过东西和南北的方向确定相关定日镜
//
bool SdBkCalcTest::checkEffectRegion(Vector3d dir, Vector3d& Hloc, Vector3d& HIloc, double diameter) {
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}

//
// [improved neighbor row] 在neighbor row操作之后，再进行一次剔除操作
//
void SdBkCalcTest::improvedNeighSdBk() {
	Timer::resetStart();
	Heliostat* cur_h = solar_scene->helios[helio_index];
	unordered_set<int> sd_set, bk_set;
	int sd_ac = 0, bk_ac = 0;
	double diameter = sqrt(pow(cur_h->helio_size.x(), 2) + pow(cur_h->helio_size.z(), 2));
	int start = 0;
	int end = 0;
	getStartEndIndex(cur_h, start, end);
//#pragma omp parallel for
	for (; start <= end; ++start) {
		auto h = solar_scene->helios[start];
		if (checkEffectRegion(Vector3d( 1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(-1, 0, 0), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d( 0, 0, 1), cur_h->helio_pos, h->helio_pos, diameter) ||
			checkEffectRegion(Vector3d(0, 0, -1), cur_h->helio_pos, h->helio_pos, diameter)) {
			if (checkCenterDist(cur_h, h, sd_dir)) {
				if (gt_sd_helio_index.count(start)) {
//#pragma omp critical
					++sd_ac;
				}
				sd_set.insert(start);

			}
			if (checkCenterDist(cur_h, h, bk_dir)) {
				if (gt_bk_helio_index.count(start)) {
//#pragma omp critical
					++bk_ac;
				}	
				bk_set.insert(start);

			}
		}

	}

	Timer::printDuration("improvedNeighSdBk");
	saveTestRes("iNeighRow", sd_ac, bk_ac, sd_set, bk_set);
}


//
// [improved neighbor row] improved neighbor row辅助函数，用于计算投影坐标点，并剔除无关定日镜
//
bool SdBkCalcTest::checkCenterDist(Heliostat* H, Heliostat* HI, Vector3d& dir) {
	vector<Vector3d> helio_v, local_v, tmp_v(4);
	vector<Vector2d> project_v(4);
	double t;
	helio_v = H->vertex;
	Vector3d reverse_dir = Vector3d(-dir.x(), -dir.y(), -dir.z());

	for (int i = 0; i < helio_v.size(); i++)
		local_v.push_back(GeometryFunc::mulMatrix(helio_v[i], H->world2localM));


	vector<Vector3d> pro(4);
	helio_v = HI->vertex;
	int cnt = 0;
	for (int i = 0; i < helio_v.size(); i++) {
		t = GeometryFunc::calcIntersection(H->helio_normal, H->helio_pos, helio_v[i], reverse_dir, pro[i]);
		if (t < Epsilon)
			return false;
		//pro[i].x() = helio_v[i].x() + t*reverse_dir.x();
		//pro[i].y() = helio_v[i].y() + t*reverse_dir.y();
		//pro[i].z() = helio_v[i].z() + t*reverse_dir.z();
	}
	double phi_x = INT_MIN;
	double phi_y = INT_MIN;
	vector<Vector3d> local_proj_v;
	for (auto v : pro) {
		local_proj_v.push_back(GeometryFunc::mulMatrix(v, H->world2localM));
		phi_x = max(local_proj_v.back().x(), phi_x);
		phi_y = max(local_proj_v.back().z(), phi_y);
	}
	
	double h_l = H->helio_size.x();
	double h_w = H->helio_size.z();
	phi_x = 2 * phi_x / h_l;
	phi_y = 2 * phi_y / h_w;
	for (auto& v : local_proj_v) {
		if (abs(v.x()) > (phi_x + 1) / 2 * h_l || abs(v.z()) > (phi_y + 1) / 2 * h_w) return false;
	}
	return true;
		
}


//
// [3DDDA] 使用3DDDA计算光线与周围定日镜的关系
//		使用3DDDA确定相关定日镜所在网格，使用bounding sphere剔除无关定日镜
void SdBkCalcTest::use3dddaSdBk() {
	Timer::resetStart();
	unordered_set<int> sd_set, bk_set;
	int sd_ac = 0;
	int bk_ac = 0;
	
	checkEstimateHelio(sd_dir, sd_set, sd_ac, gt_sd_helio_index, true);
	checkEstimateHelio(bk_dir, bk_set, bk_ac, gt_bk_helio_index, false);

	Timer::printDuration("use3dddaSdBk");
	saveTestRes("3DDDA_bk", sd_ac, bk_ac, sd_set, bk_set);
}

//
// [3DDDA] 3DDDA辅助函数
//
void SdBkCalcTest::checkEstimateHelio(Vector3d& dir, unordered_set<int>& helio_set, int& ac, unordered_set<int>& gt_helio_set, bool shadowDir) {
	Heliostat* helio = solar_scene->helios[helio_index];
	GridDDA dda_handler;
	if(shadowDir)
		dda_handler.rayCastForShadow(solar_scene, helio, -solar_scene->sunray_dir, helio_set);
	else {
		set<vector<int>> rela_grid_label;
		dda_handler.rayCastForBlock(solar_scene, helio, rela_grid_label);
		dda_handler.getBlockHelioFromGrid(solar_scene, rela_grid_label, helio_set, helio);
	}
	for (auto& index : helio_set) {
		if (gt_helio_set.count(index)) ++ac;
	}
}

void SdBkCalcTest::saveTestRes(const string& file_name, const int sd_ac, const int bk_ac, const unordered_set<int>& sd_set, const unordered_set<int>& bk_set) {
	fstream outFile(save_path + "/" + file_name + "_pr.txt", ios_base::app | ios_base::out);
	outFile << helio_index << ' ' << sd_ac << ' ' << sd_set.size()  << ' ' << gt_sd_helio_index.size() << ' ' << bk_ac << ' ' << bk_set.size() << ' ' << gt_bk_helio_index.size() << endl;
	outFile.close();

	outFile.open(save_path + "/" + file_name + "_sd_index.txt", ios_base::app | ios_base::out);
	outFile << helio_index << ' ';
	for (auto&n : sd_set)
		outFile << n << ' ';
	outFile << endl;
	outFile.close();

	outFile.open(save_path + "/" +file_name + "_bk_index.txt", ios_base::app | ios_base::out);
	outFile << helio_index << ' ';
	for (auto&n : bk_set)
		outFile << n << ' ';
	outFile << endl;
	outFile.close();
}