#include "RayCastingArgument.h"

void HeliostatDeviceArgument::setHelioDeviceOrigins(const int width_slice, const int length_slice) {
	
	helio_slice_length = length_slice <= 0 ? HELIOSTAT_SLICE_LENGTH : length_slice;
	helio_slice_width = width_slice <= 0 ? HELIOSTAT_SLICE_WIDTH : width_slice;

	if (!d_helio_origins || 
		(d_helio_origins && numberOfOrigions != helio_slice_length * helio_slice_width)) {
		if(d_helio_origins) cudaFree(d_helio_origins);
		numberOfOrigions = helio_slice_length * helio_slice_width;
		cudaMalloc((void**)&d_helio_origins, sizeof(int2)*numberOfOrigions);
	}

	int2* host_origins = new int2[numberOfOrigions];
	for (int i = 0; i < helio_slice_length; ++i)
		for (int j = 0; j < helio_slice_width; ++j)
			host_origins[i*helio_slice_width + j] = make_int2(i, j);

	cudaMemcpy(d_helio_origins, host_origins, sizeof(int2)*numberOfOrigions, cudaMemcpyHostToDevice);

	delete[] host_origins;
}

void HeliostatDeviceArgument::setHelioDevicePos(vector<Heliostat*>& helios)
{
	if (d_helio_pos) cudaFree(d_helio_pos);
	
	cudaMalloc((void**)&d_helio_pos, sizeof(float3)*helios.size());
	
	float3* h_helio_pos = new float3[helios.size()];
	for (int i=0; i<helios.size(); ++i)
		h_helio_pos[i] = GeometryFunc::convert3(helios[i]->helio_pos);
	cudaMemcpy(d_helio_pos, h_helio_pos, sizeof(float3)*helios.size(), cudaMemcpyHostToDevice);

	delete[]h_helio_pos;
}

void HeliostatDeviceArgument::setHelioDeviceArguments(vector<Heliostat*>& helios, Receiver& recv)
{
	if (!d_helio_normals ||
		(d_helio_normals && numberOfHeliostats != helios.size())) {
		if (d_helio_normals) cudaFree(d_helio_normals);
		numberOfHeliostats = helios.size();
		cudaMalloc((void**)&d_helio_normals, sizeof(float3)*numberOfHeliostats);
		cudaMalloc((void**)&d_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats);
		cudaMalloc((void**)&d_focus_index, sizeof(int)*numberOfHeliostats);
		cudaMalloc((void**)&d_imgplane_proj, sizeof(float2) * 4 * recv.recv_vertex.size()* numberOfHeliostats);
	}

	float3* h_helio_normals = new float3[numberOfHeliostats];
	float3* h_helio_vertexes = new float3[4 * numberOfHeliostats];
	int* h_focus_index = new int[numberOfHeliostats];
	float4* h_imgplane_world2local = new float4[4 * numberOfHeliostats];
	for (int i = 0; i < helios.size(); ++i) {
		h_helio_normals[i] = GeometryFunc::convert3(helios[i]->helio_normal);
		h_focus_index[i] = helios[i]->focus_center_index;
		for (int j = 0; j < 4; ++j)
			h_helio_vertexes[4 * i + j] = GeometryFunc::convert3(helios[i]->vertex[j]);



		Vector3d focus_center = recv.focus_center[h_focus_index[i]];
		Vector3d reverse_dir = (helios[i]->helio_pos - focus_center).normalized();		// The normal of image plane
		Matrix4d world2localM, local2worldM;
		GeometryFunc::getImgPlaneMatrixs(reverse_dir, focus_center, local2worldM, world2localM, 1);
		for (int j = 0; j < 4; ++j)		// TODO: check whether the matrix is right
			h_imgplane_world2local[4 * i + j] = make_float4(world2localM(j, 0), world2localM(j, 1), world2localM(j, 2), world2localM(j, 3));
	}
	cudaMemcpy(d_helio_normals, h_helio_normals, sizeof(float3)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_helio_vertexes, h_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_focus_index, h_focus_index, sizeof(int)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_imgplane_proj, h_imgplane_world2local, sizeof(float2) * 4 * recv.recv_vertex.size()* numberOfHeliostats, cudaMemcpyHostToDevice);

	delete[] h_helio_normals;
	delete[] h_helio_vertexes;
	delete[] h_focus_index;
	delete[] h_imgplane_world2local;
}

void ReceiverDeviceArgument::setRecvDeviceArguments(Receiver& recv)
{
	numberOfReceivers = recv.recv_vertex.size();
	cudaMalloc((void**)&d_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers);
	cudaMalloc((void**)&d_recv_focus_pos, sizeof(float3)*numberOfReceivers);
	cudaMalloc((void**)&d_recv_normal, sizeof(float3)*numberOfReceivers);

	float3* h_recv_vertexes = new float3[4 * numberOfReceivers];
	float3* h_recv_focus_pos = new float3[numberOfReceivers];
	float3* h_recv_normal = new float3[numberOfReceivers];

	for (int i = 0; i < recv.recv_vertex.size(); ++i) {
		for (int j = 0; j < 4; j++)
			h_recv_vertexes[4 * i + j] = GeometryFunc::convert3(recv.recv_vertex[i][j]);
		h_recv_focus_pos[i] = GeometryFunc::convert3(recv.focus_center[i]);
		h_recv_normal[i] = GeometryFunc::convert3(recv.recv_normal_list[i]);
	}

	cudaMemcpy(d_recv_vertexes, h_recv_vertexes, sizeof(float3) * 4 * numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_focus_pos, h_recv_focus_pos, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);
	cudaMemcpy(d_recv_normal, h_recv_normal, sizeof(float3)*numberOfReceivers, cudaMemcpyHostToDevice);

	delete[] h_recv_vertexes;
	delete[] h_recv_focus_pos;
	delete[] h_recv_normal;
}
