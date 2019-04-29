#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Clipper/clipper.hpp"
using namespace ClipperLib;

class PolygonCenterCalculator {
public:
	static double calcPolygonCenter(const Paths& solution, float2& center);
private:
	static double calcArea(const vector<float2>& vertexes);
	static float2 calcTriangleCenter(const vector<float2>& vertexes);
};