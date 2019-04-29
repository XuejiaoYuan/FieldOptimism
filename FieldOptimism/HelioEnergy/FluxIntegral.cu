#include "FluxIntegral.cuh"

__global__ void calcFieldFluxIntegral(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_total_energy, const int m, const int n) {
	float res = calcFluxIntegralCore(h_args, r_args, gl, m, n);
	if (res < Epsilon) return;
	atomicAdd(d_total_energy, res);
}

__global__ void calcHelioFluxIntegral(IntegralHelioDeviceArgumet h_args, ReceiverDeviceArgument r_args, GaussLegendre gl, float* d_helio_energy, const int m, const int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats*r_args.numberOfReceivers) return;

	float res = calcFluxIntegralCore(h_args, r_args, gl, m, n);

	int helioIndex = myId / (m*n*r_args.numberOfReceivers);
	atomicAdd(d_helio_energy + helioIndex, res);
}

//__device__ float calcSigma(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, int helioIndex, int recvIndex)
//{
//	float3 helio_pos = h_args.d_helio_pos[helioIndex];
//	float3 recv_pos = r_args.d_recv_focus_pos[recvIndex];
//	float dis = norm(helio_pos, recv_pos);
//
//	float3 helio_normal = h_args.d_helio_normals[helioIndex];
//	float3 recv_normal = r_args.d_recv_normal[recvIndex];
//	float3 reflect_dir = normalize(recv_pos - helio_pos);
//	float3 sunray_dir = -reflect(-reflect_dir, helio_normal);
//	float cos_w = abs(dot(sunray_dir, helio_normal));
//	float cos_rev = abs(dot(recv_normal, reflect_dir));
//	
//	float d = sqrt(h_args.d_helio_size.x * h_args.d_helio_size.y);
//	float Ht = d*(1 - cos_w);
//	float Ws = Ht;
//	float sigma_ast = Ht / (4 * dis);			// already simplified
//
//	float sigma_sun = SIGMA_SUN;
//	float sigma_s = SIGMA_S;
//	float sigma_bq = pow(2 * sigma_s, 2);
//	float sigma_t = 0;
//
//	float sigma_hf = sqrt(pow(dis, 2) * (pow(sigma_sun, 2) + pow(sigma_bq, 2) + pow(sigma_ast, 2) + pow(sigma_t, 2))) / sqrt(cos_rev);
//	return sigma_hf;
//}


__device__ float calcFluxIntegralCore(IntegralHelioDeviceArgumet& h_args, ReceiverDeviceArgument& r_args, GaussLegendre& gl, const int m, const int n) {
	int myId = GeometryFunc::getThreadId();
	if (myId >= m*n*h_args.numberOfHeliostats*r_args.numberOfReceivers) return -1;

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

	float2 gauss_param = h_args.d_gauss_param[helioIndex];
	float l_w_ratio = gauss_param.x;
	float sigma = gauss_param.y;

	//float l_w_ratio = norm(proj_v[1], proj_v[0]) / norm(proj_v[3], proj_v[0]);

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

	//float sigma = calcSigma(h_args, r_args, helioIndex, recvIndex);

	float sum = gl.calcInte(tmp_x, tmp_y, sigma, l_w_ratio) * h_args.d_factor[helioIndex];// *cos_phi;

	return sum;
}