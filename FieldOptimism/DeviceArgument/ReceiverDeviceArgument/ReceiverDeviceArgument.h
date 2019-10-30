#pragma once
#include "../../Common/CommonFunc.h"
#include "../../DataStructure/SolarScene.h"
#include "../../Common/vector_arithmetic.cuh"
#include "../../Tool/Timer/Timer.h"

class ReceiverDeviceArgument {
public:
	float3* d_recv_vertexes;		// 接收器顶点坐标
	float3* d_recv_focus_pos;		// 接收器中心点坐标
	float3* d_recv_normal;			// 接收器法向量坐标
	int numberOfReceivers;
	ReceiverDeviceArgument() :numberOfReceivers(0), d_recv_vertexes(nullptr), d_recv_focus_pos(nullptr), d_recv_normal(nullptr) {}
	~ReceiverDeviceArgument() {
		d_recv_vertexes = nullptr;
		d_recv_focus_pos = nullptr;
		d_recv_normal = nullptr;
	}
	virtual void setRecvDeviceArguments(Receiver& recv, vector<Heliostat*>& helios);
	void clear() {
		if (d_recv_vertexes){
			cudaFree(d_recv_vertexes);
			d_recv_vertexes = nullptr;
		}
		if (d_recv_focus_pos) {
			cudaFree(d_recv_focus_pos);
			d_recv_focus_pos = nullptr;
		}
		if (d_recv_normal) {
			cudaFree(d_recv_normal);
			d_recv_normal = nullptr;
		}
	}
};