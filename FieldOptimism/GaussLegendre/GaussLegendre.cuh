#pragma once
#include "../Common/CommonFunc.h"

class GaussLegendre {
public:
	GaussLegendre() : d_node_row(nullptr), d_node_col(nullptr), d_weight_row(nullptr), d_weight_col(nullptr) {}
	void initNodeWeight(const int _M, const int _N);
	void calcWeight(const int n, float*& x, float*& w, const double a = -1, const double b = 1);
	__device__ __host__
	double calcInte(const float4& x, const float4& y, const double sigma, const double ratio);
	__device__ __host__
	void legendre(const double t, const double m, double&p, double& dp);
	__device__ __host__
	double jacobi(const float4& x, const float4& y, double s, double t);
	__device__ __host__
	float2 map(const float4& x, const float4& y, double s, double t);
	__device__ __host__
	double flux_func(double x, double y, const double sigma, const double ratio) {
		return exp(-0.5 / pow(sigma, 2)*(pow(x, 2) + pow(y * ratio, 2)));
	}

	float* d_node_row, *d_node_col;
	float* d_weight_row, *d_weight_col;
	int M, N;
	~GaussLegendre() {
		cudaFree(d_node_row);
		cudaFree(d_node_col);
		cudaFree(d_weight_row);
		cudaFree(d_weight_col);

		d_node_row = nullptr;
		d_node_col = nullptr;
		d_weight_row = nullptr;
		d_weight_col = nullptr;
	}
};