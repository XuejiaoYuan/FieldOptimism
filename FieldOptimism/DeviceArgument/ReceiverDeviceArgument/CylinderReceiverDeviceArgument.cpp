#include "CylinderReceiverDeviceArgument.h"

void CylinderReceiverDeviceArgument::setRecvDeviceArguments(Receiver & recv, vector<Heliostat*>& helios)
{
	numberOfReceivers = helios.size();

	float3* h_recv_vertexes = new float3[4 * numberOfReceivers];
	float3* h_recv_focus_pos = new float3[numberOfReceivers];
	float3* h_recv_normal = new float3[numberOfReceivers];

#pragma omp parallel for
	for (int i = 0; i < numberOfReceivers; ++i) {
		h_recv_focus_pos[i] = GeometryFunc::convert3(helios[i]->focus_center);
		h_recv_normal[i] = GeometryFunc::convert3((helios[i]->focus_center - recv.recv_pos).normalized());
		vector<Vector3d> recv_vertex = recv.getRecvVertex(helios[i]->focus_center)[0];
		for (int j = 0; j < 4; ++j)
			h_recv_vertexes[4 * i + j] = GeometryFunc::convert3(recv_vertex[j]);
	}

	cudaMalloc((void**)&d_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers);
	cudaMalloc((void**)&d_recv_focus_pos, sizeof(float3)*numberOfReceivers);
	cudaMalloc((void**)&d_recv_normal, sizeof(float3)*numberOfReceivers);

	cudaMemcpy(d_recv_vertexes, h_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_focus_pos, h_recv_focus_pos, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_normal, h_recv_normal, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);

	numberOfReceivers = 1;
	delete[] h_recv_vertexes;
	delete[] h_recv_focus_pos;
	delete[] h_recv_normal;

}
