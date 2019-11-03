#pragma once
#include "../../Common/CommonFunc.h"
#include "../../DataStructure/SolarScene.h"
#include "../../Common/vector_arithmetic.cuh"
#include "../../Tool/Timer/Timer.h"

class HeliostatDeviceArgument {
public:
	int2* d_helio_origins;			// 离散网格采样点局部坐标 pixel_width x piexel_height
	float3* d_helio_normals;		// 各定日镜法向量 heliosNum
	float3* d_helio_vertexes;		// 各定日镜顶点坐标 heliosNum x 4
	float3* d_helio_pos;
	float2 d_helio_size;

	int numberOfHeliostats;

	HeliostatDeviceArgument() : d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		numberOfHeliostats(0) {}
	HeliostatDeviceArgument(int nOfHelios, int slice_length, int slice_width, int helio_list_size) :
		d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		numberOfHeliostats(nOfHelios) {}
	void setHelioDevicePos(vector<Heliostat*>& helios, bool update = false);
	void setHelioDeviceArguments(vector<Heliostat*>& helios, bool update = false);

	~HeliostatDeviceArgument() {
		d_helio_origins = nullptr;
		d_helio_normals = nullptr;
		d_helio_vertexes = nullptr;
		d_helio_pos = nullptr;

	}

	void clear() {
		if (d_helio_origins) {
			cudaFree(d_helio_origins);
			d_helio_origins = nullptr;
		}
		if (d_helio_normals) {
			cudaFree(d_helio_normals);
			d_helio_normals = nullptr;
		}
		if (d_helio_vertexes) {
			cudaFree(d_helio_vertexes);
			d_helio_vertexes = nullptr;
		}
		if (d_helio_pos) {
			cudaFree(d_helio_pos);
			d_helio_pos = nullptr;
		}
	}
};


class RayCastHelioDeviceArgument :public HeliostatDeviceArgument {
public:
	int* d_rela_shadow_helio_index;		// 各定日镜对应的相关定日镜 helioNum x helioListSize
	int* d_rela_block_helio_index;		// 各定日镜对应的相关定日镜 helioNum x helioListSize
	unsigned int* d_hit_cnt;						// 各定日镜上光线被阴影或遮挡的个数
	int numberOfOrigions;
	int helio_slice_length;
	int helio_slice_width;
	int helio_list_size;

	RayCastHelioDeviceArgument() :d_rela_shadow_helio_index(nullptr), d_rela_block_helio_index(nullptr), d_hit_cnt(nullptr),
		numberOfOrigions(0), helio_slice_length(0), helio_slice_width(0), helio_list_size(DEVICE_LIST_SIZE) {}
	RayCastHelioDeviceArgument(double helio_slice, int helio_length, int helio_width, int helio_list_size) :
		d_rela_shadow_helio_index(nullptr), d_rela_block_helio_index(nullptr), d_hit_cnt(nullptr), helio_list_size(helio_list_size) {
		helio_slice_length = helio_length / helio_slice;
		helio_slice_width = helio_width / helio_slice;
		numberOfOrigions = helio_slice_length * helio_slice_width;
	}

	~RayCastHelioDeviceArgument() {
		d_hit_cnt = nullptr;
		d_rela_shadow_helio_index = nullptr;
		d_rela_block_helio_index = nullptr;
	}

	void setHelioDeviceOrigins(const double helio_slice, int helio_length, int helio_width, bool update = false);
	void clear() {
		HeliostatDeviceArgument::clear();
		if (d_rela_shadow_helio_index) {
			cudaFree(d_rela_shadow_helio_index); d_rela_shadow_helio_index = nullptr;
		}
		if (d_rela_block_helio_index) {
			cudaFree(d_rela_block_helio_index); d_rela_block_helio_index = nullptr;
		}
		if (d_hit_cnt) {
			cudaFree(d_hit_cnt); d_hit_cnt = nullptr;
		}
	}
};

class IntegralHelioDeviceArgumet : public HeliostatDeviceArgument {
public:
	int* d_focus_index;					// 各定日镜聚焦接收器平面序号 helioNum
	float4* d_imgplane_world2local;		// 各定日镜对应image plane的坐标变换结果
	float2* d_gauss_param;				// 各定日镜在image plane上的边长比, sigma
	float* d_factor;					// 各定日镜各项因子参数
	float3* d_center_bias;				// 各定日镜阴影遮挡后区域中心位置
	float sigma;						// iHFCAL积分参数
	float DNI;							// 当前时刻DNI
	IntegralHelioDeviceArgumet() : d_focus_index(nullptr), d_imgplane_world2local(nullptr), d_gauss_param(nullptr),
		d_factor(nullptr), d_center_bias(nullptr) {}
	~IntegralHelioDeviceArgumet() {
		clearArguments();
	}

	void setHelioRecvArguments(vector<Heliostat*>& helios, Receiver& recv, bool update = false);
	void setHelioCenterBias(vector<Heliostat*>& helios, bool update = false);
	void clearArguments() {
		d_focus_index = nullptr;
		d_imgplane_world2local = nullptr;
		d_gauss_param = nullptr;
		d_factor = nullptr;
		d_center_bias = nullptr;
	}
	void clear() {
		HeliostatDeviceArgument::clear();
		if (d_focus_index) {
			cudaFree(d_focus_index);
			d_focus_index = nullptr;
		}
		if (d_imgplane_world2local) {
			cudaFree(d_imgplane_world2local);
			d_imgplane_world2local = nullptr;
		}
		if (d_gauss_param) {
			cudaFree(d_gauss_param);
			d_gauss_param = nullptr;
		}
		if (d_factor) {
			cudaFree(d_factor);
			d_factor = nullptr;
		}
		if (d_center_bias) {
			cudaFree(d_center_bias);
			d_center_bias = nullptr;
		}
	}
};
