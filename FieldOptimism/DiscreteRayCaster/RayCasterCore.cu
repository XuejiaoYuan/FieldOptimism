#include "RayCasterCore.cuh"

vector<double> rayCastCore(Vector3d sunray_dir, RayCastHelioDeviceArgument& h_args, string save_path, CalcMode calc_mode) {
	unsigned int* h_hit_cnt = new unsigned int[h_args.numberOfHeliostats];
	int* h_hit_index = nullptr;
	int* d_hit_index = nullptr;
	for (int i = 0; i < h_args.numberOfHeliostats; ++i) h_hit_cnt[i] = 0;

	if (!h_args.d_hit_cnt)
		cudaMalloc((void**)&h_args.d_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats);
	cudaMemcpy(h_args.d_hit_cnt, h_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats, cudaMemcpyHostToDevice);

	if (calc_mode == TestMode) {
		h_hit_index = new int[2 * h_args.numberOfHeliostats * h_args.numberOfOrigions];
		for (int i = 0; i < 2 * h_args.numberOfHeliostats*h_args.numberOfOrigions; ++i) h_hit_index[i] = -1;
		cudaMalloc((void**)&d_hit_index, sizeof(int)*2*h_args.numberOfHeliostats*h_args.numberOfOrigions);
		cudaMemcpy(d_hit_index, h_hit_index, sizeof(int)*2*h_args.numberOfHeliostats*h_args.numberOfOrigions, cudaMemcpyHostToDevice);
	}

	int nThreads = 1024;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, h_args.numberOfHeliostats * h_args.numberOfOrigions);
	rayCollisionCalc << <nBlocks, nThreads >> > (GeometryFunc::convert3(sunray_dir), h_args, d_hit_index);
	cudaDeviceSynchronize();

	cudaMemcpy(h_hit_cnt, h_args.d_hit_cnt, sizeof(unsigned int)*h_args.numberOfHeliostats, cudaMemcpyDeviceToHost);

	vector<double> sdbk_res;
	for (int i = 40; i < h_args.numberOfHeliostats; ++i)
		sdbk_res.push_back(h_hit_cnt[i] / (double)h_args.numberOfOrigions);

	delete[] h_hit_cnt;
	cudaFree(h_args.d_hit_cnt);
	h_args.d_hit_cnt = nullptr;

	if (calc_mode == TestMode) {
		cudaMemcpy(h_hit_index, d_hit_index, sizeof(unsigned int)*2*h_args.numberOfHeliostats, cudaMemcpyDeviceToHost);
		fstream outFile(save_path + "RC_sdbk_gt.txt", ios_base::out);
		for (int i = 40; i < h_args.numberOfHeliostats; ++i) {
			unordered_set<int> shadow_index, block_index;
			for(int j=0; j<h_args.helio_slice_length; ++j)
				for (int k = 0; k < h_args.helio_slice_width; ++k) {
					int shadow_i = h_hit_index[i*h_args.numberOfOrigions*2 + j*h_args.helio_slice_width*2 + k];
					int block_i = h_hit_index[i*h_args.numberOfOrigions * 2 + j*h_args.helio_slice_width * 2 + k + 1];
					if (shadow_i!= -1) shadow_index.insert(shadow_i);
					if (block_i != -1) block_index.insert(block_i);
				}
			for (auto& cur_i : shadow_index)
				outFile << cur_i << ' ';
			outFile << endl;
			for (auto& cur_i : block_index)
				outFile << cur_i << ' ';
			outFile << endl;
		}
		outFile.close();

		delete[] h_hit_index;
		cudaFree(d_hit_index);
	}

	return sdbk_res;
}

__global__ void rayCollisionCalc(float3 ray_dir, RayCastHelioDeviceArgument h_args, int* d_hit_index) {
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

	if (collision(helioIndex, originPos, ray_dir, h_args, true, &(d_hit_index[2*myId+1])))
		return;
	else
		collision(helioIndex, originPos, ray_dir, h_args, false, &(d_hit_index[2*myId]));
}


///
// [碰撞检测] 测试下检测最近的碰撞定日镜；实际计算检测到碰撞立即返回
///
__device__ bool collision(int helioIndex, float3 originPos, float3 sunray_dir, RayCastHelioDeviceArgument& h_args, bool shadowDir, int* hit_index) {
	int i = 0;
	float3 curNormal = h_args.d_helio_normals[helioIndex];
	float3 ray_dir = shadowDir ? -sunray_dir : reflect(sunray_dir, curNormal);

	int* d_rela_helio_index = shadowDir ? h_args.d_rela_shadow_helio_index : h_args.d_rela_block_helio_index;
	double tmin = INT_MAX;
	int hNearIndex = -1;
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
				if (!hit_index) {
					atomicAdd(&(h_args.d_hit_cnt[helioIndex]), 1);
					return true;
				}
				else if (hit_index && tmin > t) {
					tmin = t;
					hNearIndex = relaIndex;
				}
			}
		}
		else break;
	}

	if (hit_index && hNearIndex != -1) {
		atomicAdd(&(h_args.d_hit_cnt[helioIndex]), 1);
		unsigned int i = h_args.d_hit_cnt[helioIndex];
		*hit_index = hNearIndex;
		return true;
	}
	return false;
}

