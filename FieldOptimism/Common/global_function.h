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
		// double res[4] = { 0 };
		//for (int i = 0; i<4; i++)
		//	for (int j = 0; j<4; j++)
		//		res(i) += v(j) * matrix[j][i];

		if (res(3) > Epsilon)
			res /= res(3);
		return Vector3d(res(0), res(1), res(2));

	}

	inline double calcIntersection(const Vector3d& normal, const Vector3d& origin_p, const Vector3d& v, const Vector3d& dir, Vector3d& inter_v) {
		double div = dir.dot(normal);
		double t = (origin_p - v).dot(normal)/div;
		inter_v = v + dir*t;
		return t;
	}

	inline bool inProjArea(const vector<Vector3d>& v, const Vector3d& p) {
		Vector3d pre_n, cur_n, edg, line;
		for (int l = 0; l < 4; l++) {
			edg = v[(l + 1) % 4] - v[l];
			line = p-v[l];
			if (l == 0) {
				pre_n = edg.cross(line);
				continue;
			}
			cur_n = edg.cross(line);
			if (pre_n.dot(cur_n) < Epsilon)
				return false;
			pre_n = cur_n;
		}
		return true;
	}

	inline void setLocalVertex(const double l, const double w, vector<Vector3d>& vertex) {
		double half_l = l / 2.0;
		double half_w = w / 2.0;
		vertex.clear();
		vertex.push_back(Vector3d(-half_l, 0, -half_w));
		vertex.push_back(Vector3d(-half_l, 0, +half_w));
		vertex.push_back(Vector3d(+half_l, 0, +half_w));
		vertex.push_back(Vector3d(+half_l, 0, -half_w));
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

};

#endif // GEOMETRY_FUNC_H

