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
	RectangularRecvType, CylinderRecvType, PolyhedronRecvType
}ReceiverType;

class Receiver {
public:
	Receiver() {}
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
	Vector2i rows_cols;
	virtual void initRecv(json& config);
	virtual void initRecv(fstream& in, InputMode& input_mode);
	
	virtual vector<Vector3d> getFocusCenter() {
		return focus_center;
	}
	virtual vector<vector<Vector3d>> getRecvVertex(Vector3d& focus_center = Vector3d(0, 0, 0)) {
		return recv_vertex;
	}
	virtual vector<Vector3d> getNormalList(Vector3d& focus_center = Vector3d(0,0,0)) {
		return recv_normal_list;
	}

protected:
	vector<Vector3d> focus_center;				//The focus center of the receiver
	vector<Vector3d> recv_normal_list;
	vector<vector<Vector3d>> recv_vertex;
	virtual void readRecvFromJson(json& config);
	virtual void readRecvFromScn(fstream& in, InputMode& input_mode);
	vector<Vector3d> getRecvVertexCore(Vector3d& center, double half_l, double half_w, Vector3d& recv_normal);
	virtual void initRecvCore();
};

class RectangularRecv :public Receiver {
public:
	RectangularRecv() :Receiver(RectangularRecvType) {};
};

class CylinderRecv :public Receiver {
public:
	CylinderRecv() :Receiver(CylinderRecvType) {};
	void initRecv(json& config);
	void initRecv(fstream& in, InputMode& input_mode);
	Vector3d getFocusCenter(Vector3d& helio_pos);
	vector<vector<Vector3d>> getRecvVertex(Vector3d& focus_center = Vector3d(0, 0, 0));
	Vector2i getRowsCols(const double recv_pixel_length) { return Vector2i(recv_size.x() / recv_pixel_length, PI*recv_size.z() / recv_pixel_length); }
	vector<Vector3d> getNormalList(Vector3d& focus_center = Vector3d(0, 0, 0)) { return{(focus_center - recv_pos).normalized() };}

private:
	virtual void readRecvFromJson(json& config);
};


class PolyhedronRecv :public Receiver {
public:
	PolyhedronRecv() : Receiver(PolyhedronRecvType) {};
	void initRecv(json& config);
	//void initRecv(fstream& in, InputMode& input_mode);

private:
	void initRecvCore();
};

class ReceiverCreator {
public:
	static Receiver* getReceiver(const ReceiverType& recv_type) {
		switch (recv_type) {
		case RectangularRecvType:
			return new RectangularRecv();
		case CylinderRecvType:
			return new CylinderRecv();
		case PolyhedronRecvType:
			return new PolyhedronRecv();
		default:
			return new RectangularRecv();
		}
	}
};
#endif //HELIOSHADOW_RECEIVER_H
