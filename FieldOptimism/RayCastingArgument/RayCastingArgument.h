#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../Common/vector_arithmetic.cuh"

class HeliostatDeviceArgument {
public:
	float3* d_helio_origins;		// 离散网格采样点局部坐标 pixel_width x piexel_height
	float3* d_helio_normals;		// 各定日镜法向量 heliosNum
	float3* d_helio_vertexes;		// 各定日镜顶点坐标 heliosNum x 4
	float2* d_helio_shadow_bounding_origins;		// 各定日镜包围盒上下边界起点坐标（边界沿入射光线反方向）
	float2* d_helio_block_bounding_origins;		// 各定日镜包围盒上下边界起点坐标（边界沿反射光线方向）
	float3* d_helio_pos;
	int numberOfOrigions;
	int numberOfHeliostats;
	int helio_slice_length;
	int helio_slice_width;
	float2 helio_size;
	HeliostatDeviceArgument() : d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		d_helio_shadow_bounding_origins(nullptr), d_helio_block_bounding_origins(nullptr),
		numberOfHeliostats(0), numberOfOrigions(HELIOSTAT_SLICE_LENGTH* HELIOSTAT_SLICE_WIDTH), 
		helio_slice_length(HELIOSTAT_SLICE_LENGTH), helio_slice_width(HELIOSTAT_SLICE_WIDTH){}
	HeliostatDeviceArgument(int nOfHelios, int slice_length, int slice_width, float2& helio_size):
		d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		d_helio_shadow_bounding_origins(nullptr), d_helio_block_bounding_origins(nullptr),
		numberOfHeliostats(nOfHelios), numberOfOrigions(slice_length*slice_width), 
		helio_slice_length(slice_length), helio_slice_width(slice_width), helio_size(helio_size) {}
	void setHelioDeviceOrigins(const int slice_length, const int slice_width, float2 helio_size);
	void setHelioDevicePos(vector<Heliostat*>& helios);
	void setHelioDeviceArguments(vector<Heliostat*>& helios);
	void setBoundingOrigins(vector<Heliostat*>& helios, float3 dir, bool shadowDir, bool init = false);

	~HeliostatDeviceArgument() {
		cudaFree(d_helio_origins);
		cudaFree(d_helio_normals);
		cudaFree(d_helio_vertexes);
		cudaFree(d_helio_pos);
		cudaFree(d_helio_shadow_bounding_origins);
		cudaFree(d_helio_block_bounding_origins);

		d_helio_origins = nullptr;
		d_helio_normals = nullptr;
		d_helio_vertexes = nullptr;
		d_helio_pos = nullptr;
		d_helio_shadow_bounding_origins = nullptr;
		d_helio_block_bounding_origins = nullptr;
	}
};

class LayoutDeviceArgument {
public:
	float3 layout_size;
	float3 layout_bound_pos;
	float3 helio_interval;
	float3 layout_first_helio_center;
	int2 layout_row_col;

	LayoutDeviceArgument(float3 _layout_sz, float3 layout_bd_ps, float3 h_inter, float3 layout_h_c, int2 layout_row_col) :
		layout_size(_layout_sz), layout_bound_pos(layout_bd_ps), helio_interval(h_inter), layout_first_helio_center(layout_h_c), layout_row_col(layout_row_col) {};
};