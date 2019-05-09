#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Clipper/clipper.hpp"
using namespace ClipperLib;

class PolygonCenterCalculator {
public:
	static double calcPolygonCenter(const Paths& solution, Vector3d& center, bool calcCenter = true);
	static double calcArea(const vector<float2>& vertexes);

private:
	static float2 calcTriangleCenter(const vector<float2>& vertexes);
};