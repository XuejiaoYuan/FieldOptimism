#include "FluxIntegral.cuh"

__global__ void fluxIntegral(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats*r_args.numberOfReceivers) return;

	int helioIndex = myId / (m*n*r_args.numberOfReceivers);
	int recvIndex = (myId % (m*n*r_args.numberOfReceivers)) / (m*n);
	int row_col = (myId % (m*n*r_args.numberOfReceivers)) % (m*n);
	int i = row_col / n;
	int j = row_col % n;

	float3 recv_pos = r_args.d_recv_focus_pos[recvIndex];
	float3 recv_normal = r_args.d_recv_normal[recvIndex];
	float3 imgplane_normal = normalize(h_args.d_helio_pos[helioIndex] - recv_pos);
	float cos_phi = dot(recv_normal, imgplane_normal);
	if (cos_phi < Epsilon) return;

	float3 reverse_dir = imgplane_normal;		// The normal of image plane
	float3* recv_v = r_args.d_recv_vertexes + 4 * recvIndex;
	float4* imgplane_m = h_args.d_imgplane_world2local + 4 * helioIndex;
	float2 proj_v[4];
	float3 inter_v;
	for (int i = 0; i < 4; ++i) {
		GeometryFunc::calcIntersection(reverse_dir, recv_pos, recv_v[i], reverse_dir, inter_v);
		inter_v = GeometryFunc::multMatrix(inter_v, imgplane_m);
		proj_v[i] = make_float2(inter_v.x, inter_v.z);
	}

	float2 row_gap = (proj_v[3] - proj_v[0]) / m;
	float2 col_gap = (proj_v[1] - proj_v[0]) / n;

	float l_w_ratio = h_args.d_lw_ratio[helioIndex];
	float4 tmp_x = make_float4(
		(proj_v[0] + i*row_gap + j*col_gap).x,
		(proj_v[0] + (i + 1)*row_gap + j*col_gap).x,
		(proj_v[0] + (i + 1)*row_gap + (j + 1)*col_gap).x,
		(proj_v[0] + i*row_gap + (j + 1)*col_gap).x
	);

	float4 tmp_y = make_float4(
		(proj_v[0] + i*row_gap + j*col_gap).y,
		(proj_v[0] + (i + 1)*row_gap + j*col_gap).y,
		(proj_v[0] + (i + 1)*row_gap + (j + 1)*col_gap).y,
		(proj_v[0] + i*row_gap + (j + 1)*col_gap).y
	);

	float sum = gl.calcInte(tmp_x, tmp_y, h_args.sigma, l_w_ratio) * h_args.d_factor[helioIndex] * cos_phi;

	atomicAdd(d_total_energy, sum);
}