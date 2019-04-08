#include "RayCastingArgument.h"

void HeliostatDeviceArgument::setHelioDeviceOrigins(const int width_slice, const int length_slice) {
	
	helio_slice_length = length_slice <= 0 ? HELIOSTAT_SLICE_LENGTH : length_slice;
	helio_slice_width = width_slice <= 0 ? HELIOSTAT_SLICE_WIDTH : width_slice;

	if (!d_helio_origins || 
		(d_helio_origins && numberOfOrigions != helio_slice_length * helio_slice_width)) {
		if(d_helio_origins) cudaFree(d_helio_origins);
		numberOfOrigions = helio_slice_length * helio_slice_width;
		//cudaMallocManaged(&d_helio_origins, sizeof(float3)*numberOfOrigions);
		//cudaMalloc((void**)&d_helio_origins, sizeof(float3)*numberOfOrigions);
		cudaMalloc((void**)&d_helio_origins, sizeof(int2)*numberOfOrigions);
	}

	//double length_interval = helio_size.x / helio_slice_length;
	//double width_interval = helio_size.y / helio_slice_width;

	//double start_x = -helio_size.x / 2.0;
	//double start_y = -helio_size.y / 2.0;
	//float3* host_origins = new float3[numberOfOrigions];
	int2* host_origins = new int2[numberOfOrigions];
	for (int i = 0; i < helio_slice_length; ++i)
		for (int j = 0; j < helio_slice_width; ++j)
			//d_helio_origins[i*helio_slice_width + j] = make_float3(start_x + i*length_interval, start_y + j*width_interval, 0);
			//host_origins[i*helio_slice_width + j] = make_float3(start_x + i*length_interval, start_y + j*width_interval, 0);
			host_origins[i*helio_slice_width + j] = make_int2(i, j);

	//cudaMemcpy(d_helio_origins, host_origins, sizeof(float3)*numberOfOrigions, cudaMemcpyHostToDevice);
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

void HeliostatDeviceArgument::setHelioDeviceArguments(vector<Heliostat*>& helios)
{
	if (!d_helio_normals ||
		(d_helio_normals && numberOfHeliostats != helios.size())) {
		if (d_helio_normals) cudaFree(d_helio_normals);
		numberOfHeliostats = helios.size();
		cudaMalloc((void**)&d_helio_normals, sizeof(float3)*numberOfHeliostats);
		cudaMalloc((void**)&d_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats);
	}

	float3* h_helio_normals = new float3[numberOfHeliostats];
	float3* h_helio_vertexes = new float3[4 * numberOfHeliostats];
	for (int i = 0; i < helios.size(); ++i) {
		h_helio_normals[i] = GeometryFunc::convert3(helios[i]->helio_normal);
		for (int j = 0; j < 4; ++j)
			h_helio_vertexes[4 * i + j] = GeometryFunc::convert3(helios[i]->vertex[j]);
	}
	cudaMemcpy(d_helio_normals, h_helio_normals, sizeof(float3)*numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaMemcpy(d_helio_vertexes, h_helio_vertexes, sizeof(float3) * 4 * numberOfHeliostats, cudaMemcpyHostToDevice);

	delete[] h_helio_normals;
	delete[] h_helio_vertexes;
}

