//
// Created by Amber on 2018/4/3.
//
// Layout
// Define the layout of the heliostats field
//

#ifndef HELIOSHADOW_LAYOUT_H
#define HELIOSHADOW_LAYOUT_H
#pragma once

#include "../DataStructure/Heliostat.h"
#include "../DataStructure/Receiver.h"

#include "../Tool/ArgumentParser/ArgumentParser.h"



class Layout {
public:
    Layout(const LayoutType&_layout_type){
        layout_type = _layout_type;
		helio_interval = Vector3d(0, 0, 0);
		layout_bound_pos = Vector3d(0, 0, 0);
        layout_size = Vector3d(0, 0, 0);
		layout_row_col = Vector2i(0, 0);
    }
	virtual void createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios);
	void storeHelioToLayout(vector<Heliostat*>& helios);

	Vector3d helio_interval;					//Interval between heliostat's center
    int real_helio_num;								//The real number of heliostat in the field.(Optimization result)
    LayoutType layout_type;						//Heliostat field's layout type
	Vector3d layout_bound_pos;					// The bounding box of layout
	Vector3d layout_first_helio_center;			// The first heliostat center's position in the field
	Vector3d layout_size;						//Size of the layout, length/thickness/width
	Vector2i layout_row_col;					//The rows and cols of the layout
	vector<vector<vector<Heliostat*>>> helio_layout;				//List the index of heliostats in the field

protected:
	void initLayoutParams();
	void storeHelioToLayoutCore(Heliostat* helio);
	virtual void loadFieldArgs(json& field_args, double& z_start, int& rows, int& cols);
};

//
// 180417 Only consider rectangular field in one side
//
class RectLayout:public Layout{
public:
    RectLayout():Layout(RectLayoutType){};
};

class CrossRectLayout :public Layout {
public:
	CrossRectLayout() :Layout(CrossRectLayoutType){};
	void createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios);
};

class FermatLayout:public Layout{
public:
    FermatLayout():Layout(FermatLayoutType){}
	void createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios);

	vector<MatrixXd> getHelioIndex() { cerr << "Not implement!" << endl; return{}; };
	double dsep;						// 定日镜包围盒安全距离
	vector<double> helio_gap;

private:
	void setCircleHelios(Heliostat& helio_temp, const double R, const double gap, const int rows, const double angle_delta, vector<Heliostat*>& helios, const vector<Receiver*>& recvs);
	void calcCircleParams(vector<double>& recv_dis, vector<int>& n_rows, vector<double>& angle_delta, json& field_args, double dm);
};

class RadialLayout:public Layout{
public:
    RadialLayout():Layout(RadialLayoutType){}
};

class LayoutCreator{
public:
    static Layout* getLayout(const LayoutType& layout_type){
        switch(layout_type){
            case RectLayoutType:
                return new RectLayout();
			case CrossRectLayoutType:
				return new CrossRectLayout();
            case FermatLayoutType:
                return new FermatLayout();
            case RadialLayoutType:
                return new RadialLayout();
            default:
                return nullptr;
        }
    }
};

#endif //HELIOSHADOW_LAYOUT_H
