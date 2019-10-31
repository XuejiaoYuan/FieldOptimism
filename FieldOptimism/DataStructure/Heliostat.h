//
// Created by Amber on 2018/4/3.
//

#ifndef HELIOSHADOW_HELIOSTAT_H
#define HELIOSHADOW_HELIOSTAT_H
#pragma once

#include "../DataStructure/Receiver.h"


typedef enum {
	RectangularHelioType, ParaboloidHelioType
}HelioType;

typedef enum {
	HFLCAL, Gracia, iHFLCAL, bHFLCAL
}ModelType;

class Receiver;

class Heliostat {
public:
	Heliostat(const HelioType& _type) :helio_type(_type), helio_normal(Vector3d(0, 0, 0)), sd_bk(0), sigma(0) {};
	void initHelio(json& config);
	void initFluxParam(const vector<Receiver*>& recvs);
	void changeSurfaceNormal(const Vector3d& sunray_dir, const ModelType& type, const bool& calc_sigma);

	vector<Vector3d> vertex;		//Heliostat's vertex
	HelioType helio_type;			//Heliostat's type
	Vector3d helio_pos;           //The position of the heliostat's center
	Vector3d helio_poly_pos;		//The poly postion of the heliostat's center
	Vector3d helio_size;          //Heliostat's size:length, thickness, width 
	Vector2d helio_gap;           //Heliostat's slice gap: x, z
	Vector2i helio_matrix;          //Heliostat's slice matrix: row, col
	Vector3d helio_normal;          //Heliostat's surface normal
	unsigned int helio_index;			//Heliostat's index in the field
	unsigned int focus_center_index;	//Focus center index of reciver
	Matrix4d local2worldM;		//Heliostat's transform matrixs
	Matrix4d world2localM;
	double sd_bk;				// Heliostat's shadow and block ratio
	double mAA;					// Heliostat's atomospheric attenuation factor
	double cos_w;				// Heliostat's incidence cosine efficiency
	vector<double> cos_phi;		// Receiver cosine efficiencies
	double S;					// Heliostat's surface area
	double l_w_ratio;			// Heliostat's projection length and width on image plane
	double sigma;				// Heliostat's sigma
	double flux_param;			// flux_param = 0.5 * S * cos_w * rou * l_w_ration * mAA / pi
	float3 centerBias;			// 由于阴影遮挡导致的重心偏移位置
	double rotate_theta;			// 定日镜投影到image plane产生的轴旋转角度
	Vector3d focus_center;			// 定日镜在接收器平面聚焦中心

protected:
	vector<double> sigma_list;	// dis, sigma_sun, sigma_bq, sigma_ast, sigma_t, cos_rev

	void setHelioVertex();
	vector<double> setFocusCenterIndex(const vector<Receiver*>& recvs);
	void calcFluxParam(const ModelType& type, const bool& calc_sigma);
	double calcSigma();

};

class HeliostatCreator {
public:
	static Heliostat* getHeliostat(const HelioType& type) {
		return new Heliostat(type);
	}
};

#endif //HELIOSHADOW_HELIOSTAT_H
