#include "RayCasterCore.cuh"

vector<double> rayCastCore(Vector3d sunray_dir, RayCastHelioDeviceArgument& h_args) {
	cudaError_t cudaStatus;
	cudaSetDevice(0);
	unsigned int* h_hit_cnt = new unsigned int[h_args.numberOfHeliostats];
	for (int i = 0; i < h_args.numberOfHeliostats; ++i) h_hit_cnt[i] = 0;

	int nThreads = 1024;
	dim3 nBlocks;
	if (!h_args.d_hit_cnt)
		cudaStatus = cudaMalloc((void**)&h_args.d_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "1 cudamemcpy launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}

	cudaMemcpy(h_args.d_hit_cnt, h_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats, cudaMemcpyHostToDevice);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "1 cudamemcpy launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, h_args.numberOfHeliostats * h_args.numberOfOrigions);
	rayCollisionCalc << <nBlocks, nThreads >> > (GeometryFunc::convert3(sunray_dir), h_args);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}

	cudaMemcpy(h_hit_cnt, h_args.d_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats, cudaMemcpyDeviceToHost);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudamemcpy launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}

	vector<double> sdbk_res;
	for (int i = 0; i < h_args.numberOfHeliostats; ++i)
		sdbk_res.push_back(h_hit_cnt[i] / (double)h_args.numberOfOrigions);

	delete[] h_hit_cnt;
	cudaFree(h_args.d_hit_cnt);
	h_args.d_hit_cnt = nullptr;
	return sdbk_res;
}

__global__ void rayCollisionCalc(float3 ray_dir, RayCastHelioDeviceArgument h_args) {
	long long myId = GeometryFunc::getThreadId();

	if (myId >= h_args.numberOfHeliostats*h_args.numberOfOrigions) return;

	// 1. 获取光线起始点坐标
	int helioIndex = myId / h_args.numberOfOrigions;
	int originIndex = myId % h_args.numberOfOrigions;
	int2 originRowCol = h_args.d_helio_origins[originIndex];
	float3 vertexes[4] = {
		h_args.d_helio_vertexes[4 * helioIndex],
		h_args.d_helio_vertexes[4 * helioIndex + 1],
		h_args.d_helio_vertexes[4 * helioIndex + 2],
		h_args.d_helio_vertexes[4 * helioIndex + 3]
	};

	float3 rowGap = (vertexes[1] - vertexes[0]) / h_args.helio_slice_length;
	float3 colGap = (vertexes[3] - vertexes[0]) / h_args.helio_slice_width;
	float3 originPos = vertexes[0] + (0.5+originRowCol.x)*rowGap + (0.5+originRowCol.y)*colGap;

	if (collision(helioIndex, originPos, ray_dir, h_args, true))
		return;
	else
		collision(helioIndex, originPos, ray_dir, h_args, false);
}


__device__ bool collision(int helioIndex, float3 originPos, float3 sunray_dir, RayCastHelioDeviceArgument& h_args, bool shadowDir) {
	int i = 0;
	float3 curNormal = h_args.d_helio_normals[helioIndex];
	float3 ray_dir = shadowDir ? -sunray_dir : reflect(sunray_dir, curNormal);

	int* d_rela_helio_index = shadowDir ? h_args.d_rela_shadow_helio_index : h_args.d_rela_block_helio_index;
	while (i < h_args.helio_list_size) {
		int relaIndex = d_rela_helio_index[helioIndex*h_args.helio_list_size + i];
		++i;
		if (relaIndex != -1) {
			float3 relaNormal = h_args.d_helio_normals[relaIndex];
			float3 relaPos = h_args.d_helio_pos[relaIndex];
			float3 inter_v;
			double t = GeometryFunc::calcIntersection(relaNormal, relaPos, originPos, ray_dir, inter_v);
			if (t < Epsilon) continue;
			if (GeometryFunc::inProjArea(&h_args.d_helio_vertexes[4 * relaIndex], inter_v)) {
				atomicAdd(&(h_args.d_hit_cnt[helioIndex]), 1);
				return true;
			}
		}
		else return false;
	}
	return false;
}