#pragma once
#ifndef GEOMETRY_FUNC_H
#define GEOMETRY_FUNC_H	
#include "CommonFunc.h"


namespace GeometryFunc
{
	inline Vector3d mulMatrix(const Vector3d&vertex, const Matrix4d& matrix, bool point = true) {
		RowVector4d v(vertex.x(), vertex.y(), vertex.z(), 1);		// vertex: 1; vector: 0 
		if (!point)
			v.w() = 0;

		RowVector4d res = v*matrix;
		if (res(3) > Epsilon)
			res /= res(3);
		return Vector3d(res(0), res(1), res(2));
	}

	inline Vector2d multMatrix(const Vector2d& vertex, const Matrix3d& matrix, bool point = true) {
		RowVector3d v(vertex.x(), vertex.y(), 1);
		if (!point)
			v.z() = 0;

		RowVector3d res = v*matrix;
		if (res(2) > Epsilon)
			res /= res(2);
		return Vector2d(res(0), res(1));
	}

	__device__ 
		inline float3 multMatrix(const float3& vertex, const float4* matrix, bool point = true) {
		double v[4] = { vertex.x, vertex.y, vertex.z, point ? 1 : 0 };
		float4 res = make_float4(0,0,0,0);

		for (int i = 0; i < 4; ++i)		// TODO check wether right
			res += v[i] * matrix[i];

		if (res.w > Epsilon) res /= res.w;
		return make_float3(res.x, res.y, res.z);
	}

	__device__
		inline float2 multMatrix(const float2& vertex, const float3* matrix, bool point = true) {
		double v[3] = { vertex.x, vertex.y, point ? 1 : 0 };
		float3 res = make_float3(0, 0, 0);

		for (int i = 0; i < 3; ++i)		// TODO check wether right
			res += v[i] * matrix[i];

		if (res.z > Epsilon) res /= res.z;
		return make_float2(res.x, res.y);

	}

	inline float3 convert3(const Vector3d& a) {
		return make_float3(a.x(), a.y(), a.z());
	}

	inline float2 convert2(const Vector2d& a) {
		return make_float2(a.x(), a.y());
	}

	inline int2 convert2(const Vector2i& a) {
		return make_int2(a.x(), a.y());
	}

	inline Vector2d fminf(const Vector2d& a, const Vector2d& b) {
		return Vector2d(min(a.x(), b.x()), min(a.y(), b.y()));
	}

	inline Vector2d fmaxf(const Vector2d& a, const Vector2d& b) {
		return Vector2d(max(a.x(), b.x()), max(a.y(), b.y()));
	}

	__device__ __host__
		inline double calcIntersection(float3& normal, float3& origin, float3& v, float3& dir, float3& inter_v) {
		double div = dot(dir, normal);
		double t = dot(origin - v, normal) / div;
		inter_v = v + dir*t;
		return t;
	}

	inline double calcIntersection(const Vector3d& normal, const Vector3d& origin_p, const Vector3d& v, const Vector3d& dir, Vector3d& inter_v) {
		double div = dir.dot(normal);
		double t = (origin_p - v).dot(normal) / div;
		inter_v = v + dir*t;
		return t;
	}


	//inline bool inProjArea(const vector<Vector3d>& v, const Vector3d& p) {
	//	Vector3d pre_n, cur_n, edg, line;
	//	int n = v.size();
	//	for (int l = 0; l < n; l++) {
	//		edg = v[(l + 1) % n] - v[l];
	//		line = p - v[l];
	//		if (l == 0) {
	//			pre_n = edg.cross(line);
	//			continue;
	//		}
	//		cur_n = edg.cross(line);
	//		if (pre_n.dot(cur_n) < Epsilon)
	//			return false;
	//		pre_n = cur_n;
	//	}
	//	return true;
	//}
	
	inline int dcmp(double x) {
		if (fabs(x) < 1e-6)
			return 0;
		else
			return x < 0 ? -1 : 1;
	}

	inline bool onSegment(const Vector3d& p1, const Vector3d& p2, const Vector3d& q) {
		Vector3d k1 = p1 - q;
		Vector3d k2 = p2 - q;
		return dcmp(k1.x()*k2.y() - k2.x()*k1.y()) == 0 && dcmp((p1 - q).dot(p2 - q)) <= 0;
	}

	inline bool inProjArea(const vector<Vector3d>& v, const Vector3d& p) {
		bool flag = false; 
		Vector3d P1, P2; 
		int n = v.size();
		for (int i = 0, j = n-1; i < n; j = i++)
		{
			P1 = v[i];
			P2 = v[j];
			if (onSegment(P1, P2, p)) return true; 
			if ((dcmp(P1.y() - p.y())>0 != dcmp(P2.y() - p.y())>0) && dcmp(p.x() - (p.y() - P1.y())*(P1.x() - P2.x()) / (P1.y() - P2.y()) - P1.x())<0)
				flag = !flag;
		}
		return flag;
	}

	__device__ __host__ 
	inline bool inProjArea(const float3* v, const float3& p) {
		float3 pre_n, cur_n, edg, line;
		for (int l = 0; l < 4; l++) {
			edg = v[(l + 1) % 4] - v[l];
			line = p - v[l];
			if (l == 0) {
				pre_n = cross(edg, line);
				continue;
			}
			cur_n = cross(edg, line);
			if (dot(pre_n, cur_n) < Epsilon)
				return false;
			pre_n = cur_n;
		}
		return true;
	}

	inline void setLocalVertex(const Vector3d& size, vector<Vector3d>& vertex) {
		double half_l = size.x() / 2.0;
		double half_w = size.z() / 2.0;
		double half_h = size.y() / 2;
		vertex.clear();
		vertex.push_back(Vector3d(-half_l, half_h, -half_w));
		vertex.push_back(Vector3d(-half_l, half_h, +half_w));
		vertex.push_back(Vector3d(+half_l, half_h, +half_w));
		vertex.push_back(Vector3d(+half_l, half_h, -half_w));
	}

	inline void getHelioMatrix(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM) {
		Vector3d u[3];	// could be shared

		u[1] = normal;

		//if (abs(u[1].x()) > abs(u[1].z()))
		//{
		//	u[2] = u[1].cross(Vector3d(0.0f, 1.0f, 0.0f)).normalized();
		//	u[0] = u[1].cross(u[2]).normalized();
		//}
		//else
		//{
		Vector3d tmp_u(0.0f, 1.0f, 0.0f);
		u[0] = tmp_u.cross(u[1]).normalized();
		u[2] = u[0].cross(u[1]).normalized();
		//}

		for (int i = 0; i < 3; i++) {
			local2worldM(i, 0) = u[i].x();
			local2worldM(i, 1) = u[i].y();
			local2worldM(i, 2) = u[i].z();
			local2worldM(i, 3) = 0;
		}
		local2worldM(3, 0) = origin_p.x();
		local2worldM(3, 1) = origin_p.y();
		local2worldM(3, 2) = origin_p.z();
		local2worldM(3, 3) = 1;

		world2localM = local2worldM.inverse();
	}

	inline void getImgPlaneMatrixs(const Vector3d& normal, const Vector3d& origin_p, Matrix4d& local2worldM, Matrix4d& world2localM, unsigned int mode = 0) {
		Vector3d u[3];	// could be shared
		Vector3d tmp_u(0.0f, 1.0f, 0.0f);
		u[1] = normal;
		// ����image plane����ϵ�任����
		u[0] = tmp_u.cross(u[1]).normalized();
		u[2] = u[0].cross(u[1]).normalized();

		for (int i = 0; i < 3; i++) {
			local2worldM(i, 0) = u[i].x();
			local2worldM(i, 1) = u[i].y();
			local2worldM(i, 2) = u[i].z();
			local2worldM(i, 3) = 0;
		}
		local2worldM(3, 0) = origin_p.x();
		local2worldM(3, 1) = origin_p.y();
		local2worldM(3, 2) = origin_p.z();
		local2worldM(3, 3) = 1;

		world2localM = local2worldM.inverse();
	}

	inline void getTransformMatrix(const vector<Vector2d>& origin_p, const vector<Vector2d>& trans_p, Matrix3d& transformM) {
		double delta = origin_p[0].x()*origin_p[1].y() - origin_p[1].x()*origin_p[0].y();
		transformM(1, 0) = (origin_p[0].x()*trans_p[1].x() - trans_p[0].x()*origin_p[1].x()) / delta;
		cout << transformM(1, 0) << endl;
		transformM(0, 0) = (trans_p[1].x() - transformM(1, 0)* origin_p[0].y()) / origin_p[0].x();
		transformM(1, 1) = (origin_p[0].x()*trans_p[1].y() - trans_p[0].y()*origin_p[1].x()) / delta;
		transformM(0, 1) = (trans_p[1].y() - transformM(1, 1)* origin_p[0].y()) / origin_p[0].x();
		transformM(2, 2) = 1;
	}

	__host__ __device__
		inline bool setThreadsBlocks(dim3 &nBlocks, int nThreads, size_t size, bool threadFixed) {
		if (size > MAX_ALL_THREADS) {
			printf("There are too many threads to cope with, please use less threads.\n");
			return false;
		}

		int block_lastDim = (size + nThreads - 1) / nThreads;
		if (block_lastDim < MAX_BLOCK_SINGLE_DIM) {
			nBlocks.x = block_lastDim;
			nBlocks.y = nBlocks.z = 1;
			return true;
		}

		block_lastDim = (block_lastDim + MAX_BLOCK_SINGLE_DIM - 1) / MAX_BLOCK_SINGLE_DIM;
		if (block_lastDim < MAX_BLOCK_SINGLE_DIM) {
			nBlocks.x = MAX_BLOCK_SINGLE_DIM;
			nBlocks.y = block_lastDim;
			nBlocks.z = 1;
			return true;
		}
		else {
			nBlocks.x = nBlocks.y = MAX_BLOCK_SINGLE_DIM;
			nBlocks.z = (block_lastDim + MAX_BLOCK_SINGLE_DIM - 1) / MAX_BLOCK_SINGLE_DIM;
			return true;
		}
	}

	__host__ __device__
		inline bool setThreadsBlocks(dim3 &nBlocks, int &nThreads, const size_t &size) {
		//nThreads = (MAX_THREADS <= size) ? MAX_THREADS : size;
		return setThreadsBlocks(nBlocks, nThreads, size, true);
	}

	__host__ __device__
		inline unsigned long long getThreadId() {
		// unique block index inside a 3D block grid
		const unsigned long long int blockId = blockIdx.x
			+ blockIdx.y * gridDim.x
			+ gridDim.x * gridDim.y *blockIdx.z;

		// global unique thread index, block dimesion uses only x-coordinate
		const unsigned long long int threadId = blockId * blockDim.x + threadIdx.x;
		return threadId;
	}

	__host__ __device__
		inline Vector3d reflect(const Vector3d& i, const Vector3d& n) {
		return (i - 2.0*n*n.dot(i)).normalized();
	}

	__host__ __device__
	inline double min(const double& a, const double& b) {
		return a < b ? a : b;
	}

	__host__ __device__
	inline double max(const double& a, const double& b) {
		return a < b ? b : a;
	}

	
};

#endif // GEOMETRY_FUNC_H

