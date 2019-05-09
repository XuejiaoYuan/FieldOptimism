#include "CalcPolygonCenter.h"

double PolygonCenterCalculator::calcPolygonCenter(const Paths & solution, Vector3d& center, bool calcCenter)
{
	vector<vector<Vector3d>> v_list(solution.size());
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			v_list[i].push_back(Vector3d(solution[i][j].X / (double)VERTEXSCALE, 0, solution[i][j].Y / (double)VERTEXSCALE));
		}
	}
	double S = 0;
	center = Vector3d(0, 0, 0);
	for (int i = 0; i < solution.size(); ++i) {
		int k = v_list[i].size();
		for (int j = 0; j < k; ++j) {
			double tS = (v_list[i][j].x()*v_list[i][(j + 1) % k].z() - v_list[i][j].z()*v_list[i][(j + 1) % k].x()) / 2.;
			S += tS;
			if(calcCenter)
				center += tS*(v_list[i][j] + v_list[i][(j + 1) % k]) / 3;
		}
	}
	if(calcCenter) center /= S;

	return fabs(S);
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
