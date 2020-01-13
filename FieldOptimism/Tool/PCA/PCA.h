#pragma once

#include "../../Common/CommonFunc.h"

class PCAHandler
{
public:
	static double calcAxis(const vector<Vector2d>& pts, vector<Vector2d>& eigen_vecs, vector<double>& eigen_val);

};

