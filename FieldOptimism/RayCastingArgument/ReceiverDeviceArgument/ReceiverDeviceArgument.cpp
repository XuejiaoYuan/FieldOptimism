#include "ReceiverDeviceArgument.h"

void ReceiverDeviceArgument::setRecvDeviceArguments(Receiver& recv, vector<Heliostat*>& helios)
{
	vector<vector<Vector3d>> recv_vertex;
	vector<Vector3d> recv_normal_list, focus_center;

	recv_vertex = recv.get_recv_vertex();
	recv_normal_list = recv.get_normal_list();
	focus_center = recv.get_focus_center();
	numberOfReceivers = recv_vertex.size();

	cudaMalloc((void**)&d_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers);
	cudaMalloc((void**)&d_recv_focus_pos, sizeof(float3)*numberOfReceivers);
	cudaMalloc((void**)&d_recv_normal, sizeof(float3)*numberOfReceivers);

	float3* h_recv_vertexes = new float3[4 * numberOfReceivers];
	float3* h_recv_focus_pos = new float3[numberOfReceivers];
	float3* h_recv_normal = new float3[numberOfReceivers];

#pragma omp parallel for
	for (int i = 0; i < recv_vertex.size(); ++i) {
		for (int j = 0; j < 4; j++)
			h_recv_vertexes[4 * i + j] = GeometryFunc::convert3(recv_vertex[i][j]);
		h_recv_focus_pos[i] = GeometryFunc::convert3(focus_center[i]);
		h_recv_normal[i] = GeometryFunc::convert3(recv_normal_list[i]);
	}

	cudaMemcpy(d_recv_vertexes, h_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_focus_pos, h_recv_focus_pos, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_normal, h_recv_normal, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);

	delete[] h_recv_vertexes;
	delete[] h_recv_focus_pos;
	delete[] h_recv_normal;
}

