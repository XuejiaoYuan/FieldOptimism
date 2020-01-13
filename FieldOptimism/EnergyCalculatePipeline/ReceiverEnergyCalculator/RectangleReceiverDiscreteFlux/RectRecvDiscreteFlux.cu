#include "RectRecvDiscreteFlux.cuh"

void calcRectRecvFluxDistribution(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, float* d_recv_flux, int row, int col) {	
	int nThreads = 512;
	dim3 nBlocks;
	GeometryFunc::setThreadsBlocks(nBlocks, nThreads, r_args.numberOfReceivers*h_args.numberOfHeliostats*row*col);

	calcRectRecvDiscreteFlux << <nBlocks, nThreads >> > (h_args, r_args, gl, d_recv_flux, row, col);
	cudaDeviceSynchronize();

}


__global__ void calcRectRecvDiscreteFlux(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_recv_flux, int row, int col) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= r_args.numberOfReceivers*h_args.numberOfHeliostats*row*col) return;

	float res = calcRectRecvDiscreteFluxCore(h_args, r_args, gl, row, col);

	int recvIndex = (myId % (row*col*r_args.numberOfReceivers)) / (row*col);
	int row_col = (myId % (row*col*r_args.numberOfReceivers)) % (row*col);
	atomicAdd(d_recv_flux + recvIndex*row*col + row_col, res);
}

__device__ float calcRectRecvDiscreteFluxCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, int row, int col) {
	int myId = GeometryFunc::getThreadId();
	int helioIndex = myId / (row*col*r_args.numberOfReceivers);
	int recvIndex = (myId % (row*col*r_args.numberOfReceivers)) / (row*col);
	int row_col = (myId % (row*col*r_args.numberOfReceivers)) % (row*col);
	int i = row_col / col;
	int j = row_col % col;

	int focus_idx = h_args.d_focus_index[helioIndex];
	float3 focus_pos = r_args.d_recv_focus_pos[focus_idx];
	float3 recv_normal = r_args.d_recv_normal[recvIndex];
	float3 imgplane_normal = normalize(h_args.d_helio_pos[helioIndex] - focus_pos);
	float cos_phi = dot(recv_normal, imgplane_normal);
	if (cos_phi < Epsilon) return;

	float3 reverse_dir = imgplane_normal;		// The normal of image plane
	float3* recv_v = r_args.d_recv_vertexes + 4 * recvIndex;
	float4* imgplane_m = h_args.d_imgplane_world2local + 4 * helioIndex;
	float2 proj_v[4];
	float3 inter_v[4];

	float3 h_center_bias = make_float3(0, 0, 0);
	float3 i_center_bias = make_float3(0, 0, 0);
	float rotate_theta = 0;
	float shear_theta = 0;
	if (h_args.d_center_bias) {
		h_center_bias = h_args.d_center_bias[helioIndex];
		GeometryFunc::calcIntersection(reverse_dir, focus_pos, h_center_bias, -reverse_dir, i_center_bias);
		i_center_bias = GeometryFunc::multMatrix(i_center_bias, imgplane_m);
		rotate_theta = h_args.d_rotate_theta[helioIndex];
		shear_theta = h_args.d_shear_theta[helioIndex];
	}

	for (int i = 0; i < 4; ++i) {
		GeometryFunc::calcIntersection(reverse_dir, focus_pos, recv_v[i], reverse_dir, inter_v[i]);
		inter_v[i] = GeometryFunc::multMatrix(inter_v[i], imgplane_m);
		proj_v[i] = make_float2(inter_v[i].x, inter_v[i].z);
	}

	float2 row_gap = (proj_v[3] - proj_v[0]) / row;
	float2 col_gap = (proj_v[1] - proj_v[0]) / col;
	float2 gauss_param = h_args.d_gauss_param[helioIndex];
	float l_w_ratio = gauss_param.x;
	float sigma = gauss_param.y;
	float2 start_v = proj_v[0] + i*row_gap + j*col_gap;
	start_v -= make_float2(i_center_bias.x, i_center_bias.z);
	float2 trans_v;
	trans_v.x = start_v.x*cos(rotate_theta) + start_v.y*sin(rotate_theta);
	trans_v.y = start_v.y*cos(rotate_theta) - start_v.x*sin(rotate_theta);

	float res = h_args.d_factor[helioIndex] * gl.flux_func(trans_v.x, trans_v.y, sigma, l_w_ratio) *cos_phi;
	return res;
}