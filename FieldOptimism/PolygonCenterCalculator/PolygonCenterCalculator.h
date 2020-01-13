#pragma once
#include "../Common/CommonFunc.h"
#include "../Tool/Clipper/clipper.hpp"
using namespace ClipperLib;
#include "../Tool/Delaunay/triangle.h"
#include "../DataStructure/Heliostat.h"


class PolygonCenterCalculator {
public:
	static double calcPolygonCenter(int h_idx, const Paths& solution, Vector3d& center, bool calcCenter = true);
	static double calcPolygonCenter_Delaunay(Heliostat* helio, const Paths & solution, Vector3d& center, bool calcCenter);
	static double calcArea(const vector<float2>& vertexes);

private:
	static float2 calcTriangleCenter(const vector<float2>& vertexes);
};