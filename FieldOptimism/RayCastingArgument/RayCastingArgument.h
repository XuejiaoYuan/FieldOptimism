#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/SolarScene.h"
#include "../Common/vector_arithmetic.cuh"
#include "../DataStructure/Timer.h"

class HeliostatDeviceArgument {
public:
	int2* d_helio_origins;			// ��ɢ���������ֲ����� pixel_width x piexel_height
	float3* d_helio_normals;		// �����վ������� heliosNum
	float3* d_helio_vertexes;		// �����վ��������� heliosNum x 4
	float3* d_helio_pos;
	
	int numberOfHeliostats;
	
	HeliostatDeviceArgument() : d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		numberOfHeliostats(0){}
	HeliostatDeviceArgument(int nOfHelios, int slice_length, int slice_width, int helio_list_size):
		d_helio_origins(nullptr), d_helio_normals(nullptr), d_helio_vertexes(nullptr), d_helio_pos(nullptr),
		numberOfHeliostats(nOfHelios) {}
	void setHelioDevicePos(vector<Heliostat*>& helios,  bool update = false);
	void setHelioDeviceArguments(vector<Heliostat*>& helios, bool update = false);

	~HeliostatDeviceArgument() {
		cudaFree(d_helio_origins);
		cudaFree(d_helio_normals);
		cudaFree(d_helio_vertexes);
		cudaFree(d_helio_pos);

		d_helio_origins = nullptr;
		d_helio_normals = nullptr;
		d_helio_vertexes = nullptr;
		d_helio_pos = nullptr;
	}

};


class RayCastHelioDeviceArgument :public HeliostatDeviceArgument {
public:
	int* d_rela_shadow_helio_index;		// �����վ���Ӧ����ض��վ� helioNum x helioListSize
	int* d_rela_block_helio_index;		// �����վ���Ӧ����ض��վ� helioNum x helioListSize
	unsigned int* d_hit_cnt;						// �����վ��Ϲ��߱���Ӱ���ڵ��ĸ���
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
		cudaFree(d_rela_shadow_helio_index);
		cudaFree(d_rela_block_helio_index);
		if (d_hit_cnt) {
			cudaFree(d_hit_cnt);
			d_hit_cnt = nullptr;
		}
		d_rela_shadow_helio_index = nullptr;
		d_rela_block_helio_index = nullptr;
	}

	void setHelioDeviceOrigins(const double helio_slice, int helio_length, int helio_width, bool update = false);

};

class IntegralHelioDeviceArgumet : public HeliostatDeviceArgument {
public:
	int* d_focus_index;					// �����վ��۽�������ƽ����� helioNum
	float4* d_imgplane_world2local;		// �����վ���Ӧimage plane������任���
	float* d_lw_ratio;					// �����վ���image plane�ϵı߳���
	float* d_factor;					// �����վ��������Ӳ���
	float sigma;						// iHFCAL���ֲ���
	float DNI;							// ��ǰʱ��DNI
	IntegralHelioDeviceArgumet() : d_focus_index(nullptr), d_imgplane_world2local(nullptr), 
		d_lw_ratio(nullptr), d_factor(nullptr){}
	~IntegralHelioDeviceArgumet() {
		clearArguments();
	}

	void setHelioRecvArguments(vector<Heliostat*>& helios, Receiver& recv, bool update = false);
	void clearArguments() {
		cudaFree(d_focus_index);
		cudaFree(d_imgplane_world2local);
		cudaFree(d_lw_ratio);
		cudaFree(d_factor);

		d_focus_index = nullptr;
		d_imgplane_world2local = nullptr;
		d_lw_ratio = nullptr;
		d_factor = nullptr;
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


class ReceiverDeviceArgument {
public:
	float3* d_recv_vertexes;		// ��������������
	float3* d_recv_focus_pos;		// ���������ĵ�����
	float3* d_recv_normal;			// ����������������
	int numberOfReceivers;				
	ReceiverDeviceArgument() :numberOfReceivers(0), d_recv_vertexes(nullptr), d_recv_focus_pos(nullptr), d_recv_normal(nullptr) {}
	~ReceiverDeviceArgument() {
		cudaFree(d_recv_vertexes);
		cudaFree(d_recv_focus_pos);
		cudaFree(d_recv_normal);

		d_recv_vertexes = nullptr;
		d_recv_focus_pos = nullptr;
		d_recv_normal = nullptr;
	}
	void setRecvDeviceArguments(Receiver& recv);
};