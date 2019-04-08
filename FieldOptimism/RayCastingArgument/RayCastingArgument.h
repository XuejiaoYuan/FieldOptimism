#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../Common/vector_arithmetic.cuh"
#include "../DataStructure/Timer.h"

class HeliostatDeviceArgument {
public:
	int2* d_helio_origins;		// 离散网格采样点局部坐标 pixel_width x piexel_height
	float3* d_helio_normals;		// 各定日镜法向量 heliosNum
	float3* d_helio_vertexes;		// 各定日镜顶点坐标 heliosNum x 4
	float3* d_helio_pos;
	int* d_rela_shadow_helio_index;		// 各定日镜对应的相关定日镜 helioNum x helioListSize
	int* d_rela_block_helio_index;		// 各定日镜对应的相关定日镜 helioNum x helioListSize
	unsigned int* d_hit_cnt;						// 各定日镜上光线被阴影或遮挡的个数
	int numberOfOrigions;
	int numberOfHeliostats;
	int helio_slice_length;
	int helio_slice_width;
	int helio_list_size;	
	float2 helio_size;
	HeliostatDeviceArgument() : d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		d_rela_shadow_helio_index(nullptr), d_rela_block_helio_index(nullptr),
		numberOfHeliostats(0), numberOfOrigions(HELIOSTAT_SLICE_LENGTH* HELIOSTAT_SLICE_WIDTH), 
		helio_slice_length(HELIOSTAT_SLICE_LENGTH), helio_slice_width(HELIOSTAT_SLICE_WIDTH), helio_list_size(DEVICE_LIST_SIZE){}
	HeliostatDeviceArgument(int nOfHelios, int slice_length, int slice_width, float2& helio_size, int helio_list_size):
		d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		d_rela_shadow_helio_index(nullptr), d_rela_block_helio_index(nullptr),
		numberOfHeliostats(nOfHelios), numberOfOrigions(slice_length*slice_width), 
		helio_slice_length(slice_length), helio_slice_width(slice_width), helio_size(helio_size), helio_list_size(helio_list_size) {}
	void setHelioDeviceOrigins(const int slice_length, const int slice_width);
	void setHelioDevicePos(vector<Heliostat*>& helios);
	void setHelioDeviceArguments(vector<Heliostat*>& helios);

	~HeliostatDeviceArgument() {
		cudaFree(d_helio_origins);
		cudaFree(d_helio_normals);
		cudaFree(d_helio_vertexes);
		cudaFree(d_helio_pos);
		cudaFree(d_rela_shadow_helio_index);
		cudaFree(d_rela_block_helio_index);

		d_helio_origins = nullptr;
		d_helio_normals = nullptr;
		d_helio_vertexes = nullptr;
		d_helio_pos = nullptr;
		d_rela_shadow_helio_index = nullptr;
		d_rela_block_helio_index = nullptr;
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