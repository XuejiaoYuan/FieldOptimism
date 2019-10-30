#include "GridDDA.h"


bool GridDDA::checkBoundingBox(const Vector3d & Hloc, const Vector3d & Hnormal, const Vector3d & HIloc, const Vector3d & dir, double diameter)
{
	Vector3d dist = HIloc - Hloc;
	double proj = dist.dot(dir);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}

void GridDDA::rayCastForBlock(SolarScene * solar_scene, Heliostat * helio, set<vector<int>>& rela_grid_index)
{
	Vector3d bk_dir = (helio->focus_center - helio->helio_pos).normalized();
	double radius = GridDDACore(bk_dir, helio, solar_scene->layouts[0], rela_grid_index);
}

void GridDDA::getBlockHelioFromGrid(SolarScene * solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Heliostat* helio)
{
	Vector3d dir = (helio->focus_center - helio->helio_pos).normalized();
	getHelioFromGridCore(solar_scene, rela_grid_label, rela_helio_index, dir, helio);
}

void GridDDA::rayCastForShadow(SolarScene* solar_scene, Heliostat * helio, Vector3d dir, unordered_set<int>& rela_helio_index)
{
	set<vector<int>> relative_grid_label;
	double radius = GridDDACore(dir, helio, solar_scene->layouts[0], relative_grid_label);
	getHelioFromGridCore(solar_scene, relative_grid_label, rela_helio_index, dir, helio);
}


double GridDDA::GridDDACore(Vector3d& dir, Heliostat* helio, Layout* layout, set<vector<int>>& relative_grid_label)
{
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
	//set<vector<int>> relative_helio_label;
	for (int i = 0; i < 2; ++i) {
		// 4.0 设置光线范围
		int minCol = static_cast<int>((boundBox[2 * i].x() - layout->layout_first_helio_center.x()) / cellDimension.x());
		int minRow = static_cast<int>((boundBox[2 * i].y() - layout->layout_first_helio_center.z()) / cellDimension.y());
		int maxCol = static_cast<int>((boundBox[2 * i + 1].x() - layout->layout_first_helio_center.x()) / cellDimension.x() + 0.5);
		int maxRow = static_cast<int>((boundBox[2 * i + 1].y() - layout->layout_first_helio_center.z()) / cellDimension.y() + 0.5);

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
				relative_grid_label.insert({ grid_label.x(), grid_label.y() });

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

	return radius;
}

void GridDDA::getHelioFromGridCore(SolarScene * solar_scene, set<vector<int>>& rela_grid_label, unordered_set<int>& rela_helio_index, Vector3d & dir, Heliostat * helio)
{
	double radius = sqrt(pow(helio->helio_size.x(), 2) + pow(helio->helio_size.z(), 2)) / 2.0;
	for (auto& iter = rela_grid_label.begin(); iter != rela_grid_label.end(); ++iter) {
		for (auto& rela_helio : solar_scene->layouts[0]->helio_layout[(*iter)[0]][(*iter)[1]]) {
			if (rela_helio == helio || rela_helio_index.count(rela_helio->helio_index)) continue;
			if (checkBoundingBox(helio->helio_pos, helio->helio_normal, rela_helio->helio_pos, dir, 2 * radius)) {
				rela_helio_index.insert(rela_helio->helio_index);
			}
		}
	}
}

