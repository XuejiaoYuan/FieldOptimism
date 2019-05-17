#pragma once
//
// Created by Amber on 2018/4/3.
//
//  Receiver
//  Define the receiver in the solar system.
//

#ifndef HELIOSHADOW_RECEIVER_H
#define HELIOSHADOW_RECEIVER_H

#include "../Common/global_function.cuh"
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"

typedef enum {
	RectangularRecvType, CylinderRecvType, CircularTruncatedConeRecvType, PolyhedronRecvType
}ReceiverType;

class Receiver {
public:
	Receiver(const ReceiverType&_recv_type) {
		recv_type = _recv_type;
		recv_pos = Vector3d(0, 0, 0);
		recv_size = Vector3d(0, 0, 0);
		recv_normal = Vector3d(0, 0, 0);
		recv_face = 0;
	}

	ReceiverType recv_type;					//Receiver's type
	Vector3d recv_pos;                      //The position of receiver's center
	Vector3d recv_size;                     //Receiver's size
	Vector3d recv_normal;                   //The normal of the first receiver's face
	int recv_face;							//The index of the receiver's face
	int recv_face_num;
	virtual void init_recv(fstream& inFile, InputMode& input_mode);
	void readRecvFromFile(fstream& inFile, InputMode& input_mode);
	virtual vector<Vector3d> get_focus_center() {
		return focus_center;
	}
	virtual vector<vector<Vector3d>> get_recv_vertex(Vector3d& focus_center = Vector3d(0, 0, 0)) {
		return recv_vertex;
	}
	virtual vector<Vector3d> get_normal_list() {
		return recv_normal_list;
	}

protected:
	vector<Vector3d> focus_center;				//The focus center of the receiver
	vector<Vector3d> recv_normal_list;
	vector<vector<Vector3d>> recv_vertex;

	vector<Vector3d> getRecvVertex(Vector3d& center, double half_l, double half_w, Vector3d& recv_normal);
};

class RectangularRecv :public Receiver {
public:
	RectangularRecv() :Receiver(RectangularRecvType) {};
};

class CylinderRecv :public Receiver {
public:
	CylinderRecv() :Receiver(CylinderRecvType) {};
	void init_recv(fstream& inFile, InputMode& input_mode);
	Vector3d get_focus_center(Vector3d& helio_pos);
	vector<vector<Vector3d>> get_recv_vertex(Vector3d& focus_center = Vector3d(0, 0, 0));
};

class CircularTruncatedConeRecv :public Receiver {
public:
	CircularTruncatedConeRecv() :Receiver(CircularTruncatedConeRecvType) {};
};

class PolyhedronRecv :public Receiver {
public:
	PolyhedronRecv() : Receiver(PolyhedronRecvType) {};
	void init_recv(fstream& inFile, InputMode& input_mode);
};

class ReceiverCreator {
public:
	Receiver* getReceiver(const ReceiverType& recv_type) {
		switch (recv_type) {
		case RectangularRecvType:
			return new RectangularRecv();
		case CylinderRecvType:
			return new CylinderRecv();
		case CircularTruncatedConeRecvType:
			return new CircularTruncatedConeRecv();
		case PolyhedronRecvType:
			return new PolyhedronRecv();
		default:
			return new RectangularRecv();
		}
	}
};
#endif //HELIOSHADOW_RECEIVER_H
