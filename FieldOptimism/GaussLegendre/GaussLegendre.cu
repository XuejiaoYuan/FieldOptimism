#include "GaussLegendre.cuh"

void GaussLegendre::initNodeWeight(const int _M, const int _N)
{
	M = _M;
	N = _N;

	calcWeight(M, d_node_row, d_weight_row);
	calcWeight(N, d_node_col, d_weight_col);
}

__device__
float GaussLegendre::calcInte(const float4& x, const float4& y, const float sigma, const float ratio)
{
	
	float sum = 0.0;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			float2 map_v = map(x, y, d_node_row[i], d_node_col[j]);
			sum += d_weight_row[i] * d_weight_col[j] * jacobi(x, y, d_node_row[i], d_node_col[j])*flux_func(map_v.x, map_v.y, sigma, ratio);
		}
	}
	return sum;
}


void GaussLegendre::legendre(const float t, const float m, float&p, float& dp)
{
	float p0 = 1.0;
	float p1 = t;
	for (int k = 1; k < m; k++) {
		p = ((2.0*k + 1)*t*p1 - k*p0) / (1.0 + k);
		p0 = p1;
		p1 = p;
	}
	dp = m*(p0 - t*p1) / (1.0 - t*t);
}


///
//给定积分上下限x1和x2及阶数n，返回长度为n的x_list及w_list，
//其中分别存放n点gauss - legendre求积分公式的坐标点及权重
//	:param x1 : 积分下限
//	:param x2 : 积分上限
//	:param x : 求积分公式的坐标点
//	:param w : 求积分公式的权重
//	:param n : 高斯积分阶数
void GaussLegendre::calcWeight(const int n, float*& x, float*& w, const float a, const float b) {
	float * h_node = new float[n];
	float * h_weight = new float[n];
	cudaMalloc((void**)&x, sizeof(float)*n);
	cudaMalloc((void**)&w, sizeof(float)*n);

	int nRoots = int((n + 1) / 2);
	float p, dp, dt, t;
	for (int i = 0; i < nRoots; i++) {
		t = cos(PI*(i + 0.75) / (n + 0.5));
		while (true) {
			legendre(t, n, p, dp);
			dt = -p / dp;
			t += dt;
			if (abs(dt) < Epsilon) {
				h_node[i] = -t;
				h_node[n - 1 - i] = t;
				h_weight[i] = 2.0 / (1.0 - t*t) / (dp*dp);
				h_weight[n - i - 1] = h_weight[i];
				break;
			}
		}
	}

	cudaMemcpy(x, h_node, sizeof(float)*n, cudaMemcpyHostToDevice);
	cudaMemcpy(w, h_weight, sizeof(float)*n, cudaMemcpyHostToDevice);

	delete[] h_node;
	delete[] h_weight;
}

__device__ 
float GaussLegendre::jacobi(const float4& x, const float4& y, const float s, const float t) {
	float J00 = -(1.0 - t) * x.x + (1.0 - t) * x.y + (1.0 + t) * x.z - (1.0 + t) * x.w;
	float J01 = -(1.0 - t) * y.x + (1.0 - t) * y.y + (1.0 + t) * y.z - (1.0 + t) * y.w;
	float J10 = -(1.0 - s) * x.x - (1.0 + s) * x.y + (1.0 + s) * x.z + (1.0 - s) * x.w;
	float J11 = -(1.0 - s) * y.x - (1.0 + s) * y.y + (1.0 + s) * y.z + (1.0 - s) * y.w;
	return (J00*J11 - J01*J10) / 16.0;
}


__device__
float2 GaussLegendre::map(const float4&x, const float4&y, const float s, const float t) {
	float4 N;
	N.x = (1.0 - s)*(1.0 - t) / 4.0;
	N.y = (1.0 + s)*(1.0 - t) / 4.0;
	N.z = (1.0 + s)*(1.0 + t) / 4.0;
	N.w = (1.0 - s)*(1.0 + t) / 4.0;
	float2 map_v;
	map_v.x = dot(N, x);
	map_v.y = dot(N, y);
	return map_v;
}

