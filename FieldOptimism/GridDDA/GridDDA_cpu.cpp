#include "GridDDA.h"


void GridDDA::predictRelatedHelio(SolarScene* solar_scene, RayCastHelioDeviceArgument& h_args, bool shadowDir)
{
	// 0. 参数准备
	Vector3d &sunray_dir = solar_scene->sunray_dir;
	vector<Heliostat*>& helios = solar_scene->helios;
	int* h_rela_index = new int[helios.size()*h_args.helio_list_size];

	// 1. 計算定日鏡的相關定日鏡
#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i) {
		unordered_set<int> helio_label;
		rayCastGridDDA(solar_scene, helios[i], -sunray_dir, helio_label, shadowDir);
		int j = 0;
		for (auto& iter : helio_label) {
			h_rela_index[i*h_args.helio_list_size + j] = iter;
			++j;
		}
		if (j < h_args.helio_list_size) h_rela_index[i*h_args.helio_list_size + j] = -1;		// 设置相关定日镜个数
	}

	// 2. 将数据拷贝至GPU用于后续ray cast计算
	int*& d_rela_index = shadowDir ? h_args.d_rela_shadow_helio_index : h_args.d_rela_block_helio_index;
	cudaMalloc((void**)&d_rela_index, sizeof(int)*helios.size()*h_args.helio_list_size);
	cudaMemcpy(d_rela_index, h_rela_index, sizeof(int)*helios.size()*h_args.helio_list_size, cudaMemcpyHostToDevice);
	delete[] h_rela_index;
}


bool GridDDA::checkBoundingBox(const Vector3d & Hloc, const Vector3d & Hnormal, const Vector3d & HIloc, const Vector3d & dir, double diameter)
{
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}

void GridDDA::testHandler(SolarScene* solar_scene)
{

	RayCastHelioDeviceArgument h_args;
	Vector3d helio_size = solar_scene->helios[0]->helio_size;
	h_args.setHelioDeviceOrigins(0.01, helio_size.x(), helio_size.z());		// 第一次调用cuda程序耗时较长
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);

	predictRelatedHelio(solar_scene, h_args, true);
	predictRelatedHelio(solar_scene, h_args, false);

}

void GridDDA::rayCastGridDDA(SolarScene* solar_scene, Heliostat * helio, Vector3d dir, unordered_set<int>& rela_helio_index, bool shadowDir)
{
	// 1. 获取光线方向
	Vector3d sunray_dir = -dir;
	Layout*& layout = solar_scene->layouts[0];
	Vector3d curNormal = helio->helio_normal;
	if (!shadowDir) {
		int fc_index = helio->focus_center_index;
		dir = (helio->focus_center - helio->helio_pos).normalized();
	}
	Vector2d ray_dir(dir.x(), dir.z());
	ray_dir.normalized();
	Vector2d tangent_dir(ray_dir.y(), -ray_dir.x());

	// 2. 获取光线起点
	Vector3d Hloc = helio->helio_pos;
	Vector3d helio_v[4];
	double radius = sqrt(pow(helio->helio_size.x(), 2) + pow(helio->helio_size.z(), 2)) / 2.0;
	Vector2d proj_origin[2] = {
		Vector2d(Hloc.x(), Hloc.z()) + radius*tangent_dir,
		Vector2d(Hloc.x(), Hloc.z()) - radius*tangent_dir
	};

	// 3. 确定光线在场地中y反向移动距离及最终离开网格位置
	double upper_y = layout->layout_size.y() + layout->layout_bound_pos.y();

	double dis = (upper_y - helio->helio_pos.y()) / dir.y();						// TODO: check if the vertex is right
	Vector2d upper_v[2] = {
		Vector2d((radius + dis)*dir.x(), (radius + dis)*dir.z()) + proj_origin[0],
		Vector2d((radius + dis)*dir.x(), (radius + dis)*dir.z()) + proj_origin[1]
	};

	// 4. 确定定日镜反射光线在镜场中的grid范围
	Vector2d boundBox[4] = { 
		GeometryFunc::fminf(proj_origin[0], upper_v[0]), GeometryFunc::fmaxf(proj_origin[0], upper_v[0]),
		GeometryFunc::fminf(proj_origin[1], upper_v[1]), GeometryFunc::fmaxf(proj_origin[1], upper_v[1])
	};		// min boundary, max boundary

																							// 5. 确定layout下各网格的矩阵间隔
	Vector2d cellDimension = Vector2d(layout->helio_interval.x(), layout->helio_interval.z());

	int helio_col = static_cast<int>((Hloc.x() - layout->layout_first_helio_center.x()) / cellDimension.x());			// smaller x is, smaller col is
	int helio_row = static_cast<int>((Hloc.z() - layout->layout_first_helio_center.z()) / cellDimension.y());			// smaller z is, smaller row is


	// 4. DDA求交
	set<vector<int>> relative_helio_label;
	for (int i = 0; i < 2; ++i) {
		// 4.0 设置光线范围
		int minCol = static_cast<int>((boundBox[2*i].x() - layout->layout_first_helio_center.x()) / cellDimension.x());
		int minRow = static_cast<int>((boundBox[2*i].y() - layout->layout_first_helio_center.z()) / cellDimension.y());
		int maxCol = static_cast<int>((boundBox[2*i+1].x() - layout->layout_first_helio_center.x()) / cellDimension.x() + 0.5);
		int maxRow = static_cast<int>((boundBox[2*i+1].y() - layout->layout_first_helio_center.z()) / cellDimension.y() + 0.5);

		minCol = std::max(0, minCol);
		minRow = std::max(0, minRow);
		maxCol = std::min(layout->layout_row_col.y() - 1, maxCol);
		maxRow = std::min(layout->layout_row_col.x() - 1, maxRow);

		// 4.1 计算起始参数
		Vector2d deltaT;
		Vector2d t;

		Vector2d o_grid(
			proj_origin[i].x() - layout->layout_bound_pos.x(),
			proj_origin[i].y() - layout->layout_bound_pos.z()
		);

		if (ray_dir.x() < 0) {
			deltaT.x() = -cellDimension.x() / ray_dir.x();
			t.x() = (floor(o_grid.x() / cellDimension.x())*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}
		else {
			deltaT.x() = cellDimension.x() / ray_dir.x();
			t.x() = (floor((o_grid.x() / cellDimension.x()) + 1)*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}

		if (ray_dir.y() < 0) {
			deltaT.y() = -cellDimension.y() / ray_dir.y();
			t.y() = (floor(o_grid.y() / cellDimension.y())*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}
		else {
			deltaT.y() = cellDimension.y() / ray_dir.y();
			t.y() = (floor((o_grid.y() / cellDimension.y()) + 1)*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}

		Vector2i grid_label(
			static_cast<int>((proj_origin[i].y() - layout->layout_bound_pos.z()) / cellDimension.y()),	// smaller z is, smaller row is
			static_cast<int>((proj_origin[i].x() - layout->layout_bound_pos.x()) / cellDimension.x())	// smaller x is, smaller col is
		);


		// 4.2 光线遍历格子
		while (1) {
			if (grid_label.x() < minRow || grid_label.x() > maxRow ||
				grid_label.y() < minCol || grid_label.y() > maxCol)
				break;
			else
				relative_helio_label.insert({ grid_label.x(), grid_label.y() });

			if (t.x() < t.y()) {
				t.x() += deltaT.x();
				if (ray_dir.x() < 0)
					grid_label.y()--;
				else
					grid_label.y()++;
			}
			else {
				t.y() += deltaT.y();
				if (ray_dir.y() < 0)		// smaller z is, smaller row is
					grid_label.x()--;
				else
					grid_label.x()++;
			}
		}
	}

	// 5. 计算相关定日镜
	for (auto& iter = relative_helio_label.begin(); iter != relative_helio_label.end(); ++iter) {
		for (auto& rela_helio : solar_scene->layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
			if (rela_helio == helio || rela_helio_index.count(rela_helio->helio_index)) continue;
			if (checkBoundingBox(helio->helio_pos, helio->helio_normal, rela_helio->helio_pos, dir, 2*radius)) {
				rela_helio_index.insert(rela_helio->helio_index);
			}
		}
	}

}

void GridDDA::saveTestRes(const string& file_name, int helioNum, int* d_related_index, int list_size)
{
	fstream outFile(file_name + "_gpu.txt", ios_base::app | ios_base::out);
	for (int i = 0; i < helioNum; ++i) {
		outFile << i << ' ';
		for (int j = 0; j < list_size && d_related_index[i*list_size+j]!=-1; ++j)
			outFile << d_related_index[i * list_size + j] << ' ';
		outFile << endl;
	}
	outFile.close();

}

