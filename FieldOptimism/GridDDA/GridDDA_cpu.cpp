#include "GridDDA.h"
#include "RayCastGridDDA.cuh"


void GridDDA::predictRelatedHelio(SolarScene* solar_scene, HeliostatDeviceArgument& h_args, LayoutDeviceArgument l_args, bool shadowDir)
{
	// 0. ����׼��
	Vector3d &sunray_dir = solar_scene->sunray_dir;
	vector<Heliostat*>& helios = solar_scene->helios;
	auto layout = solar_scene->layouts[0];
	double diameter = sqrt(pow(h_args.helio_size.x, 2) + pow(h_args.helio_size.y, 2));
	DeviceList<int2>* d_related_grid_list = new DeviceList<int2>[h_args.numberOfHeliostats * 2];		// ���Ӧ����δ������Ը��Ƶ�host����
	DeviceList<int>*& d_related_index = shadowDir ? d_shadow_related_index : d_block_related_index;
	

	// 1. Ӌ�㶨���R�����P�����R
	Timer t;
	rayCastGridDDA(sunray_dir, h_args, l_args, d_related_grid_list, shadowDir);

	t.printDuration("grid dda");

	// 2. ����ÿ�����վ�����ض��վ��洢�ռ�
	if (!d_related_index)
		d_related_index = new DeviceList<int>[h_args.numberOfHeliostats];


	// 3. �_��ÿ�������R�����P�����R
	t.resetStart();
	int listSize = d_related_index[0].listSize;
	for (int i = 0; i < h_args.numberOfHeliostats; ++i) {
		unordered_set<int> related_index;
		set<pair<int, int>> visited_grid;
		Heliostat* helio = helios[i];
		int& index_cnt = d_related_index[i].listLength = 0;

		for (int k = 0; k < 2; ++k) {
			int numOfPredictGrid = d_related_grid_list[2 * i + k].listLength;

			for (int j = 0; j < numOfPredictGrid; ++j) {
				pair<int, int> grid_index(d_related_grid_list[2 * i + k].d_list[j].x, d_related_grid_list[2 * i + k].d_list[j].y);
				if (!visited_grid.count(grid_index)) {
					visited_grid.insert(grid_index);

					// 3.1 ȡÿ��grid�ж�Ӧ�Ķ��վ�
					for (auto& rela_helio : solar_scene->layouts[0]->helio_layout[grid_index.first][grid_index.second]) {
						int r_index = rela_helio->helio_index;
						if (r_index == i || related_index.count(r_index)) continue;

						// 3.2 ����Ƿ����bounding sphere��Ҫ��
						if (checkBoundingBox(helio->helio_pos, helio->helio_normal, rela_helio->helio_pos, -sunray_dir, diameter, shadowDir)) {
							related_index.insert(r_index);
							if (index_cnt < listSize)
								d_related_index[i].d_list[index_cnt++] = r_index;
						}
					}
				}
			}
		}
	}
	t.printDuration("bouding box test");

	// 4. ����^��׃��
	delete[] d_related_grid_list;
}


bool GridDDA::checkBoundingBox(const Vector3d & Hloc, const Vector3d & Hnormal, const Vector3d & HIloc, const Vector3d & dir, double diameter, bool shadowDir)
{
	Vector3d dist = HIloc - Hloc;
	double proj;
	if (shadowDir) proj = dist.dot(dir);
	else proj = GeometryFunc::reflect(-dir, Hnormal).dot(dist);
	if (proj < Epsilon) return false;
	if (sqrt(pow(dist.norm(), 2) - pow(proj, 2)) > diameter) return false;
	return true;
}

void GridDDA::testHandler(SolarScene* solar_scene)
{
	Timer t;
	HeliostatDeviceArgument h_args;
	Vector2d helio_size(solar_scene->helios[0]->helio_size.x(), solar_scene->helios[0]->helio_size.z());
	h_args.setHelioDeviceOrigins(100, 100,  GeometryFunc::convert2(helio_size));
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);
	h_args.setBoundingOrigins(solar_scene->helios, GeometryFunc::convert3(-solar_scene->sunray_dir), true, true);
	h_args.setBoundingOrigins(solar_scene->helios, GeometryFunc::convert3(-solar_scene->sunray_dir), false, true);
	t.printDuration("set helio device argument");

	t.resetStart();
	auto layout = solar_scene->layouts[0];
	LayoutDeviceArgument l_args(
		GeometryFunc::convert3(layout->layout_size), 
		GeometryFunc::convert3(layout->layout_bound_pos),
		GeometryFunc::convert3(layout->helio_interval), 
		GeometryFunc::convert3(layout->layout_first_helio_center), 
		GeometryFunc::convert2(layout->layout_row_col)
	);
	t.printDuration("set layout device argument");

	predictRelatedHelio(solar_scene, h_args, l_args, true);
	predictRelatedHelio(solar_scene, h_args, l_args, false);

	saveTestRes("SdBkRes/DDA_shadow", h_args.numberOfHeliostats, d_shadow_related_index);
	saveTestRes("SdBkRes/DDA_block", h_args.numberOfHeliostats, d_block_related_index);

}

void GridDDA::saveTestRes(const string& file_name, int helioNum, DeviceList<int>* d_related_index)
{
	fstream outFile(file_name + "_gpu.txt", ios_base::app | ios_base::out);
	for (int i = 0; i < helioNum; ++i) {
		outFile << i << ' ';
		for (int j = 0; j < d_related_index[i].listLength; ++j)
			outFile << d_related_index[i].d_list[j] << ' ';
		outFile << endl;
	}
	outFile.close();

}

