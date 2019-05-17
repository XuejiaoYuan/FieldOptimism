#pragma once
#include "../Common/CommonFunc.h"


class GaussLegendre {
public:
	GaussLegendre() : d_node_row(nullptr), d_node_col(nullptr), d_weight_row(nullptr), d_weight_col(nullptr) {}
	__device__
		float calcInte(const float4& x, const float4& y, const float sigma, const float ratio);
	void legendre(const float t, const float m, float&p, float& dp);
	__device__
		float jacobi(const float4& x, const float4& y, const float s, const float t);
	__device__
		float2 map(const float4& x, const float4& y, const float s, const float t);
	__device__
		float flux_func(float x, float y, const float sigma, const float ratio) {
		return exp(-0.5 / pow(sigma, 2)*(pow(x, 2) + pow(y * ratio, 2))) / pow(sigma, 2);
	}

	float* d_node_row, *d_node_col;
	float* d_weight_row, *d_weight_col;
	int M, N;
	~GaussLegendre() {
		d_node_row = nullptr;
		d_node_col = nullptr;
		d_weight_row = nullptr;
		d_weight_col = nullptr;
	}

	void clear();
	static GaussLegendre* m_instance;
	static GaussLegendre* getInstance(int _M, int _N);

private:
	void initNodeWeight(const int _M, const int _N);
	void calcWeight(const int n, float*& x, float*& w, const float a = -1, const float b = 1);

};