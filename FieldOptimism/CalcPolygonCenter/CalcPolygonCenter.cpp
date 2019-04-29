#include "CalcPolygonCenter.h"

double PolygonCenterCalculator::calcPolygonCenter(const Paths & solution, float2& center)
{
	vector<float2> v_list;	
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			v_list.push_back(make_float2(solution[i][j].X / (double)VERTEXSCALE, solution[i][j].Y / (double)VERTEXSCALE));
		}
	}
	double area = calcArea(v_list);
	
	vector<float2> tris(3);
	tris[0] = v_list[0];
	center.x = 0;
	center.y = 0;
	for (int i = 2; i < v_list.size(); ++i) {
		tris[1] = v_list[i - 1];
		tris[2] = v_list[i];
		double curS = calcArea(tris);
		if (curS) {
			float2 curC = calcTriangleCenter(tris);
			center += curS *curC;
		}
	}
	if (area) center /= area;

	double S = 0.;          //S面积,xy横纵坐标和
	float2 s = make_float2(0, 0);
	int k = v_list.size();
	for (int i = 0; i < k; i++) {
		double tS = (v_list[i].x*v_list[(i+1)%k].y - v_list[i].y*v_list[(i + 1) % k].x) / 2.;
		S += tS;
		s += tS*(v_list[i] + v_list[(i+1)%k]) / 3;
	}
	s /= S;

	return fabs(area);
}

double PolygonCenterCalculator::calcArea(const vector<float2>& vertexes)
{
	double sum = 0;
	int n = vertexes.size();
	for (int i = 0; i < n; ++i)
		sum += vertexes[i].x*vertexes[(i + 1) % n].y - vertexes[(i + 1) % n].x*vertexes[i].y;
	return sum / 2.0;
}

float2 PolygonCenterCalculator::calcTriangleCenter(const vector<float2>& vertexes)
{
	float2 v1 = vertexes[0];
	float2 v2 = (vertexes[1] + vertexes[2]) / 2.0;
	float2 v3 = vertexes[1];
	float2 v4 = (vertexes[2] + vertexes[0]) / 2.0;

	if (v1.x == v2.x && v1.y == v2.y) {
		v1 = vertexes[2];
		v2 = (vertexes[0] + vertexes[1]) / 2.0;
	}
	else if (v3.x == v4.x && v3.y == v4.y) {
		v3 = vertexes[2];
		v4 = (vertexes[0] + vertexes[1]) / 2.0;
	}

	double a1 = (v2.y - v1.y) / (v2.x - v1.x);
	double b1 = v1.y - a1*v1.x;

	double a2 = (v4.y - v3.y) / (v4.x - v3.x);
	double b2 = v3.y - a2*v3.x;
	
	double x = (b2 - b1) / (a1 - a2);
	double y = a1*x + b1;

	return make_float2(x, y);
}
