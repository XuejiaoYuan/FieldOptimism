#include "GridDDA.h"


void GridDDA::predictRelatedHelio(SolarScene * solar_scene, bool shadowDir)
{
	// 0. 每次重新計算將統計值清空
	if (shadowDir) shadowIndexSum = 0;
	else blockIndexSum = 0;

	vector<Heliostat*>& helios = solar_scene->helios;
	int heliosNum = solar_scene->helios.size();
	vector<set<set<pair<int, int>>> > predict_block_index(heliosNum), predict_shadow_index(heliosNum);
	Vector3d dir = -solar_scene->sunray_dir;
	
	// 1. 3DDDA + bounding sphere確定每個定日鏡的相關定日鏡
	vector<vector<int>> relaIndex(helios.size());
	for (int i = 0; i < heliosNum; ++i) {
		if (!shadowDir) {
			set<pair<int, int>> estimate_block, estimate_shadow;
			int fc_index = helios[i]->focus_center_index;
			dir =  (solar_scene->recvs[0]->focus_center[fc_index] - helios[i]->helio_pos).normalized();
		}
		calcIntersection3DDDA(solar_scene, i, dir, relaIndex[i]);
		if (shadowDir) 
			shadowIndexSum += relaIndex.size();
		else
			blockIndexSum += relaIndex.size();
	}
	
	// 2. 分配GPU上每個定日鏡的對應結果
	if (shadowDir) setRelatedHelioIndex(relaIndex, d_shadow_rela_index, d_start_shadow_index, shadowIndexSum);
	else setRelatedHelioIndex(relaIndex, d_block_rela_index, d_start_block_index, blockIndexSum);

}

 
void GridDDA::calcIntersection3DDDA(SolarScene * solar_scene, const int helio_index, const Vector3d & dir, vector<int>& relaIndex)
{
	Heliostat*& helio = solar_scene->helios[helio_index];
	vector<Layout*>& layouts = solar_scene->layouts;
	vector<Vector3d> helio_v;
	vector<Vector2d> project_helio_v;
	set<pair<int, int>> relative_helio_label;

	helio_v = helio->vertex;
	for (int i = 0; i < 4; i++) {
		project_helio_v.push_back(Vector2d(helio_v[i].x(), helio_v[i].z()));
	}

	// 1. 確定光線在場地中y方向移動距離及最終離開網格的位置
	double upper_y = layouts[0]->layout_size.y() + layouts[0]->layout_bound_pos.y();
	Vector2d upper_v[4];
	for (int i = 0; i < 4; i++) {
		double dis = (upper_y - helio_v[i].y()) / dir.y();
		upper_v[i] = Vector2d(dis*dir.x(), dis*dir.z()) + project_helio_v[i];
	}

	// 2. 確定定日鏡反射光柱在鏡場中的矩形包圍盒範圍
	vector<Vector2d> boundBox(2);
	boundBox[0] = Vector2d(INT_MAX, INT_MAX);		// min boundary
	boundBox[1] = Vector2d(INT_MIN, INT_MIN);		// max boundary
	for (int i = 0; i < 4; i++) {
		boundBox[0].x() = fmin(boundBox[0].x(), fmin(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[0].y() = fmin(boundBox[0].y(), fmin(project_helio_v[i].y(), upper_v[i].y()));
		boundBox[1].x() = fmax(boundBox[1].x(), fmax(project_helio_v[i].x(), upper_v[i].x()));
		boundBox[1].y() = fmax(boundBox[1].y(), fmax(project_helio_v[i].y(), upper_v[i].y()));
	}
	Vector2d ray_dir = Vector2d(dir.x(), dir.z()).normalized();
	Vector2d deltaT;
	Vector2d t;
	Vector2d origs[2] = {
		Vector2d(INT_MAX,INT_MAX),		// min vertex
		Vector2d(INT_MIN,INT_MIN)		// max vertex
	};


	// 3. 確定layout下各網格的矩陣間隔
	Vector2d cellDimension = Vector2d(layouts[0]->helio_interval.x(), layouts[0]->helio_interval.z());

	double colLength = boundBox[1].x() - boundBox[0].x();
	double rowLength = boundBox[1].y() - boundBox[0].y();

	int helio_col = (helio->helio_pos.x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x();			// smaller x is, smaller col is
	int helio_row = (helio->helio_pos.z() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y();			// smaller z is, smaller row is

	int minCol = (boundBox[0].x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x();
	int minRow = (boundBox[0].y() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y();
	int maxCol = (boundBox[1].x() - layouts[0]->layout_first_helio_center.x()) / cellDimension.x() + 0.5;
	int maxRow = (boundBox[1].y() - layouts[0]->layout_first_helio_center.z()) / cellDimension.y() + 0.5;

	minCol = max(0, minCol);
	minRow = max(0, minRow);
	maxCol = min(layouts[0]->layout_row_col.y() - 1, maxCol);
	maxRow = min(layouts[0]->layout_row_col.x() - 1, maxRow);

	// 4. 確定定日鏡包圍盒範圍
	Vector2d tangent_dir(ray_dir.y(), -ray_dir.x());
	Vector2d center_v(helio->helio_pos.x(), helio->helio_pos.z());
	double diameter = sqrt(pow(helio->helio_size.x(), 2) + pow(helio->helio_size.z(), 2)) / 2.0;
	vector<Vector2d> start_v = { center_v + diameter*tangent_dir, center_v - diameter*tangent_dir };

	for (auto&orig : start_v) {
		Vector2d o_grid(
			orig.x() - layouts[0]->layout_bound_pos.x(),
			orig.y() - layouts[0]->layout_bound_pos.z()
		);

		if (ray_dir.x() < 0) {
			deltaT.x() = -cellDimension.x() / ray_dir.x();
			t.x() = (floor(o_grid.x() / cellDimension.x())*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}
		else {
			deltaT.x() = cellDimension.x() / ray_dir.x();
			t.x() = ((floor(o_grid.x() / cellDimension.x()) + 1)*cellDimension.x() - o_grid.x()) / ray_dir.x();
		}

		if (ray_dir.y() < 0) {
			deltaT.y() = -cellDimension.y() / ray_dir.y();
			t.y() = (floor(o_grid.y() / cellDimension.y())*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}
		else {
			deltaT.y() = cellDimension.y() / ray_dir.y();
			t.y() = ((floor(o_grid.y() / cellDimension.y()) + 1)*cellDimension.y() - o_grid.y()) / ray_dir.y();
		}

		double tmp = 0;
		Vector2i grid_label(
			(orig.y() - layouts[0]->layout_bound_pos.z()) / cellDimension.y(),	// smaller z is, smaller row is
			(orig.x() - layouts[0]->layout_bound_pos.x()) / cellDimension.x()	// smaller x is, smaller col is
		);

		while (1) {
			if (grid_label.x() < minRow || grid_label.x() > maxRow ||
				grid_label.y() < minCol || grid_label.y() > maxCol)
				break;
			else if (grid_label.x() != helio_row || grid_label.y() != helio_col) {
				pair<int, int> res((int)grid_label.x(), (int)grid_label.y());
				relative_helio_label.insert(res);
			}

			if (t.x() < t.y()) {
				tmp = t.x();
				t.x() += deltaT.x();
				if (ray_dir.x() < 0)
					grid_label.y()--;
				else
					grid_label.y()++;
			}
			else {
				tmp = t.y();
				t.y() += deltaT.y();
				if (ray_dir.y() < 0)		// smaller z is, smaller row is
					grid_label.x()--;
				else
					grid_label.x()++;
			}
		}
	}

	// 5. 使用bounding sphere剔除無關定日鏡
	set<int> helio_set;
	for (auto& iter = relative_helio_label.begin(); iter != relative_helio_label.end(); ++iter) {
		for (auto& rela_helio : solar_scene->layouts[0]->helio_layout[(*iter).first][(*iter).second]) {
			if (rela_helio == helio) continue;
			if (checkBoundingBox(helio->helio_pos, rela_helio->helio_pos, dir, diameter)) {
				if (!helio_set.count(rela_helio->helio_index)) {
					helio_set.insert(rela_helio->helio_index);
					relaIndex.push_back(rela_helio->helio_index);
				}
			}
		}
	}
}


bool GridDDA::checkBoundingBox(const Vector3d& Hloc, const Vector3d& HIloc, const Vector3d& dir, double diameter) {
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}

__global__ void rayGridDDA(const Vector3d & dir)
{
	long long myId = GeometryFunc::getThreadId();
	
	return __global__ void();
}

void GridDDA::setRelatedHelioIndex(vector<vector<int>>& rela_index, int *& d_rela_index, int *& d_start_index, int indexSum)
{
	int heliosNum = rela_index.size();
	if (d_rela_index) cudaFree(d_rela_index);
	cudaMallocManaged(&d_rela_index, sizeof(int) * (indexSum+1));						// 最後一位用於標識結束
	if (!d_start_index) cudaMallocManaged(&d_start_index, sizeof(int)*(heliosNum+1));

	int startIndex = 0;
	for (int i = 0; i < heliosNum; ++i) {
		d_start_index[i] = startIndex;
		for (int j = 0; j < rela_index[i].size(); ++j)
			d_rela_index[startIndex + j] = rela_index[i][j];
		startIndex += rela_index[i].size();
	}

}
