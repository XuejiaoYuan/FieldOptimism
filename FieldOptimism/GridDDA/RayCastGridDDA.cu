#include "RayCastGridDDA.cuh"
#include "../DataStructure/Timer.h"

__global__ void rayCastGridDDACore(float3 dir, HeliostatDeviceArgument h_args, LayoutDeviceArgument l_args, DeviceList<int2>* d_related_grid_list, bool shadowDir)
{
	//long long myId = GeometryFunc::getThreadId();

	//if (myId >= h_args.numberOfHeliostats * 2) return;

	//// 1. 获取光线方向
	//int curIndex = myId / 2;
	//float3 curNormal = h_args.d_helio_normals[curIndex];
	//if (!shadowDir) dir = reflect(-dir, curNormal);
	//float2 ray_dir = normalize(make_float2(dir.x, dir.z));
	//float2 tangent_dir = myId % 2 ? make_float2(ray_dir.y, -ray_dir.x) : make_float2(-ray_dir.y, ray_dir.x);
	//
	//// 2. 获取光线起点
	////float2 proj_origin = shadowDir ? h_args.d_helio_shadow_bounding_origins[myId] : h_args.d_helio_block_bounding_origins[myId];
	//float3 Hloc = h_args.d_helio_pos[curIndex];
	//float3 helio_v[4];
	//double radius = sqrt(pow(h_args.helio_size.x, 2) + pow(h_args.helio_size.y, 2)) / 2.0;
	////float2 project_helio_v[4];
	////for (int i = 0; i < 4; i++) {
	////	helio_v[i] = h_args.d_helio_vertexes[4 * curIndex + i];
	////	project_helio_v[i] = make_float2(helio_v[i].x, helio_v[i].z);
	////}
	////double proj_dis = fmax(norm(project_helio_v[0], project_helio_v[2]), norm(project_helio_v[1], project_helio_v[3]))/2.0;
	////float2 proj_origin = make_float2(Hloc.x, Hloc.z) + proj_dis*tangent_dir;
	//float2 proj_origin = make_float2(Hloc.x, Hloc.z) + radius*tangent_dir;

	//// 3. 确定光线在场地中y反向移动距离及最终离开网格位置
	//double upper_y = l_args.layout_size.y + l_args.layout_bound_pos.y;
	//
	//double dis = (upper_y - h_args.d_helio_vertexes[4 * curIndex + 1].y) / dir.y;						// TODO: check if the vertex is right
	//float2 upper_v = make_float2((radius+dis)*dir.x, (radius+dis)*dir.z) + proj_origin;
	////float2 upper_v = make_float2((dis + proj_dis)*dir.x, (dis + proj_dis)*dir.z) + proj_origin;

	//// 4. 确定定日镜反射光线在镜场中的grid范围
	//float2 boundBox[2] = { fminf(proj_origin, upper_v), fmaxf(proj_origin, upper_v) };		// min boundary, max boundary

	//// 5. 确定layout下各网格的矩阵间隔
	//float2 cellDimension = make_float2(l_args.helio_interval.x, l_args.helio_interval.z);
	//
	//int helio_col = static_cast<int>((Hloc.x - l_args.layout_first_helio_center.x) / cellDimension.x);			// smaller x is, smaller col is
	//int helio_row = static_cast<int>((Hloc.z - l_args.layout_first_helio_center.z) / cellDimension.y);			// smaller z is, smaller row is

	//int minCol = static_cast<int>((boundBox[0].x - l_args.layout_first_helio_center.x) / cellDimension.x);
	//int minRow = static_cast<int>((boundBox[0].y - l_args.layout_first_helio_center.z) / cellDimension.y);
	//int maxCol = static_cast<int>((boundBox[1].x - l_args.layout_first_helio_center.x) / cellDimension.x + 0.5);
	//int maxRow = static_cast<int>((boundBox[1].y - l_args.layout_first_helio_center.z) / cellDimension.y + 0.5);

	//minCol = std::max(0, minCol);	// 为了保证搜索范围，向外扩大一圈
	//minRow = std::max(0, minRow);
	//maxCol = std::min(l_args.layout_row_col.y - 1, maxCol);
	//maxRow = std::min(l_args.layout_row_col.x - 1, maxRow);


	//// 4. DDA求交
	//// 4.1 计算起始参数
	//float2 deltaT;
	//float2 t;

	//float2 o_grid = make_float2(
	//	proj_origin.x - l_args.layout_bound_pos.x,
	//	proj_origin.y - l_args.layout_bound_pos.z
	//);

	//if (ray_dir.x < 0) {
	//	deltaT.x = -cellDimension.x / ray_dir.x;
	//	t.x = (floor(o_grid.x / cellDimension.x)*cellDimension.x - o_grid.x) / ray_dir.x;
	//}
	//else {
	//	deltaT.x = cellDimension.x / ray_dir.x;
	//	t.x = (floor((o_grid.x / cellDimension.x) + 1)*cellDimension.x - o_grid.x) / ray_dir.x;
	//}

	//if (ray_dir.y < 0) {
	//	deltaT.y = -cellDimension.y / ray_dir.y;
	//	t.y = (floor(o_grid.y / cellDimension.y)*cellDimension.y - o_grid.y) / ray_dir.y;
	//}
	//else {
	//	deltaT.y = cellDimension.y / ray_dir.y;
	//	t.y = (floor((o_grid.y / cellDimension.y) + 1)*cellDimension.y - o_grid.y) / ray_dir.y;
	//}

	//int2 grid_label = make_int2(
	//	static_cast<int>((proj_origin.y - l_args.layout_bound_pos.z) / cellDimension.y),	// smaller z is, smaller row is
	//	static_cast<int>((proj_origin.x - l_args.layout_bound_pos.x) / cellDimension.x)	// smaller x is, smaller col is
	//);
	//

	//// 4.2 光线遍历格子
	//while (1) {
	//	if (grid_label.x < minRow || grid_label.x > maxRow ||
	//		grid_label.y < minCol || grid_label.y > maxCol)
	//		break;
	//	else /*if(grid_label.x != helio_row || grid_label.y != helio_col) */{
	//		if (d_related_grid_list[myId].listLength < d_related_grid_list[myId].listSize)
	//			d_related_grid_list[myId].d_list[d_related_grid_list[myId].listLength++] = grid_label;
	//	}

	//	if (t.x < t.y) {
	//		t.x += deltaT.x;
	//		if (ray_dir.x < 0)
	//			grid_label.y--;
	//		else
	//			grid_label.y++;
	//	}
	//	else {
	//		t.y += deltaT.y;
	//		if (ray_dir.y < 0)		// smaller z is, smaller row is
	//			grid_label.x--;
	//		else
	//			grid_label.x++;
	//	}
	//}
return;
}



void rayCastGridDDA(Vector3d sunray_dir, HeliostatDeviceArgument h_args, LayoutDeviceArgument l_args, DeviceList<int2>* d_related_grid_list, bool shadowDir) {
	int nThreads = 1024;
	dim3 nBlocks;
	float3 _sunray_dir = GeometryFunc::convert3(-sunray_dir);
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, h_args.numberOfHeliostats * 2);
	rayCastGridDDACore << <nBlocks, nThreads >> > (_sunray_dir, h_args, l_args, d_related_grid_list, shadowDir);
	cudaDeviceSynchronize();

	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}
