#include "HeliostatDeviceArgument.h"

void RayCastHelioDeviceArgument::setHelioDeviceOrigins(const double helio_slice, int helio_length, int helio_width, bool update) {
	helio_slice_length = helio_length / helio_slice;
	helio_slice_width = helio_width / helio_slice;

	if (update && d_helio_origins) {
		cudaFree(d_helio_origins);
		d_helio_origins = nullptr;
	}
	numberOfOrigions = helio_slice_length * helio_slice_width;
	cudaMalloc((void**)&d_helio_origins, sizeof(int2)*numberOfOrigions);


	int2* host_origins = new int2[numberOfOrigions];
	for (int i = 0; i < helio_slice_length; ++i)
		for (int j = 0; j < helio_slice_width; ++j)
			host_origins[i*helio_slice_width + j] = make_int2(i, j);

	cudaMemcpy(d_helio_origins, host_origins, sizeof(int2)*numberOfOrigions, cudaMemcpyHostToDevice);

	delete[] host_origins;
}

void HeliostatDeviceArgument::setHelioDevicePos(vector<Heliostat*>& helios, bool update)
{
	if (update && d_helio_pos) {
		cudaFree(d_helio_pos);
		d_helio_pos = nullptr;
	}
	if (!d_helio_pos)
		cudaMalloc((void**)&d_helio_pos, sizeof(float3)*helios.size());

	float3* h_helio_pos = new float3[helios.size()];
#pragma omp parallel for
	for (int i = 0; i<helios.size(); ++i)
		h_helio_pos[i] = GeometryFunc::convert3(helios[i]->helio_pos);
	cudaMemcpy(d_helio_pos, h_helio_pos, sizeof(float3)*helios.size(), cudaMemcpyHostToDevice);

	delete[]h_helio_pos;
}

void HeliostatDeviceArgument::setHelioDeviceArguments(vector<Heliostat*>& helios, bool update)
{
	if (update && d_helio_normals) {
		cudaFree(d_helio_normals);
		d_helio_normals = nullptr;
	}
	if (!d_helio_normals) {
		numberOfHeliostats = helios.size();
		cudaMalloc((void**)&d_helio_normals, sizeof(float3)*numberOfHeliostats);
		cudaMalloc((void**)&d_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats);
	}

	float3* h_helio_normals = new float3[numberOfHeliostats];
	float3* h_helio_vertexes = new float3[4 * numberOfHeliostats];
	int* h_focus_index = new int[numberOfHeliostats];
	float4* h_imgplane_world2local = new float4[4 * numberOfHeliostats];
	Vector3d focus_center, reverse_dir;
	Matrix4d world2localM, local2worldM;
#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i) {
		h_helio_normals[i] = GeometryFunc::convert3(helios[i]->helio_normal);
		h_focus_index[i] = helios[i]->focus_center_index;
		for (int j = 0; j < 4; ++j)
			h_helio_vertexes[4 * i + j] = GeometryFunc::convert3(helios[i]->vertex[j]);
	}
	cudaMemcpy(d_helio_normals, h_helio_normals, sizeof(float3)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_helio_vertexes, h_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats, cudaMemcpyHostToDevice);
	d_helio_size = make_float2(helios[0]->helio_size.x(), helios[0]->helio_size.z());

	delete[] h_helio_normals;
	delete[] h_helio_vertexes;
}

void IntegralHelioDeviceArgumet::setHelioRecvArguments(vector<Heliostat*>& helios, Receiver& recv, bool update)
{
	if (update && d_focus_index) {
		clearArguments();
	}
	if (!d_focus_index) {
		cudaMalloc((void**)&d_focus_index, sizeof(int)*numberOfHeliostats);
		cudaMalloc((void**)&d_imgplane_world2local, sizeof(float4) * 4 * numberOfHeliostats);
		cudaMalloc((void**)&d_gauss_param, sizeof(float2)*numberOfHeliostats);
		cudaMalloc((void**)&d_factor, sizeof(float)*numberOfHeliostats);
	}

	int* h_focus_index = new int[numberOfHeliostats];
	float4* h_imgplane_world2local = new float4[4 * numberOfHeliostats];
	float2* h_gauss_param = new float2[numberOfHeliostats];
	float* h_factor = new float[numberOfHeliostats];
	Vector3d reverse_dir;
	Matrix4d world2localM, local2worldM;
#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i) {
		h_focus_index[i] = helios[i]->focus_center_index;
		h_gauss_param[i] = make_float2(helios[i]->l_w_ratio, helios[i]->sigma);
		h_factor[i] = helios[i]->flux_param * (1 - helios[i]->sd_bk);
	}

	for (int i = 0; i < helios.size(); ++i) {
		reverse_dir = (helios[i]->helio_pos - helios[i]->focus_center).normalized();		// The normal of image plane
		GeometryFunc::getImgPlaneMatrixs(reverse_dir, helios[i]->focus_center, local2worldM, world2localM, 1);
		for (int j = 0; j < 4; ++j) {
			h_imgplane_world2local[4 * i + j] = make_float4(world2localM(j, 0), world2localM(j, 1), world2localM(j, 2), world2localM(j, 3));
		}
	}
	cudaMemcpy(d_focus_index, h_focus_index, sizeof(int)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_imgplane_world2local, h_imgplane_world2local, sizeof(float4) * 4 * numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_gauss_param, h_gauss_param, sizeof(float2)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_factor, h_factor, sizeof(float)*numberOfHeliostats, cudaMemcpyHostToDevice);

	delete[] h_focus_index;
	delete[] h_imgplane_world2local;
	delete[] h_gauss_param;
	delete[] h_factor;
}

void IntegralHelioDeviceArgumet::setHelioCenterBias(vector<Heliostat*>& helios, bool update)
{
	if (update && d_center_bias) {
		cudaFree(d_center_bias);
		d_center_bias = nullptr;
	}
	if (!d_center_bias)
		cudaMalloc((void**)&d_center_bias, sizeof(float3)*numberOfHeliostats);
	float3* h_center_bias = new float3[numberOfHeliostats];
#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i)
		h_center_bias[i] = helios[i]->centerBias;

	cudaMemcpy(d_center_bias, h_center_bias, sizeof(float3)*numberOfHeliostats, cudaMemcpyHostToDevice);
	delete[] h_center_bias;
}
