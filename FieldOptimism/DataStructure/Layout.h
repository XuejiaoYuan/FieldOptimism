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
		helio_interval = Vector2d(0, 0);
		layout_bound_pos = Vector2d(0, 0);
        layout_size = Vector3d(0, 0, 0);
		layout_row_col = Vector2i(0, 0);
    }
	virtual void createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios);
	void storeHelioToLayout(vector<Heliostat*>& helios);
	void initLayout(fstream& in, InputMode& input_mode);

	Vector2d helio_interval;					//Interval between heliostat's center
    int real_helio_num;								//The real number of heliostat in the field.(Optimization result)
    LayoutType layout_type;						//Heliostat field's layout type
	Vector2d layout_bound_pos;					// The bounding box of layout
	Vector2d layout_first_helio_center;			// The first heliostat center's position in the field
	Vector3d layout_size;						//Size of the layout, length/width/thickness
	Vector2i layout_row_col;					//The rows and cols of the layout
	vector<vector<vector<Heliostat*>>> helio_layout;				//List the index of heliostats in the field

protected:
	void initLayoutParams();
	void storeHelioToLayoutCore(Heliostat* helio);
	virtual void loadFieldArgs(ArgumentParser& argumentParser, json& field_args, double& z_start, int& rows, int& cols);
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

class RadialStaggerLayout:public Layout{
public:
	RadialStaggerLayout():Layout(RadialStaggerLayoutType){}
	void createHelioAndLayout(ArgumentParser& argumentParser, json& field_args, vector<Heliostat*>& helios);

	vector<MatrixXd> getHelioIndex() { cerr << "Not implement!" << endl; return{}; };
	//double dsep;						// 定日镜包围盒安全距离
	vector<double> dsep;
	vector<double> helio_gap;

private:
	bool setCircleHelios(Heliostat& h_tmp, const int idx, vector<double>& recv_dis, vector<int>& rows, const double gap,
		const int col, vector<Heliostat*>& helios, const vector<Receiver*>& recvs, const int test_helio_num);
	void calcCircleParams(int real_helio_num, vector<double>& recv_dis, vector<int>& n_rows, vector<int>& n_cols, json& field_args, double dm);
};

class SpiralLayout:public Layout{
public:
    SpiralLayout():Layout(SpiralLayoutType){}
	void createHelioAndLayout(ArgumentParser& argumentPaser, json& field_args, vector<Heliostat*>& helios);
private:
	bool checkHelioDistance(double dm, double pre_theta, double cur_theta, double pre_r, double cur_r);
	double calcHeliostatHeight(double x, double y, double h);
};

class LayoutCreator{
public:
    static Layout* getLayout(const LayoutType& layout_type){
        switch(layout_type){
            case RectLayoutType:
                return new RectLayout();
			case CrossRectLayoutType:
				return new CrossRectLayout();
            case RadialStaggerLayoutType:
                return new RadialStaggerLayout();
            case SpiralLayoutType:
                return new SpiralLayout();
            default:
                return nullptr;
        }
    }
};

#endif //HELIOSHADOW_LAYOUT_H
