#include "PolygonCenterCalculator.h"
#include "../Tool/PCA/PCA.h"

//
// [多边形重心计算接口] 计算多边形重心
//
double PolygonCenterCalculator::calcPolygonCenter(int h_idx, const Paths & solution, Vector3d& center, bool calcCenter)
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

double PolygonCenterCalculator::calcPolygonCenter_Delaunay(Heliostat* helio, const Paths & solution, Vector3d& center, bool calcCenter)
{
	double S = 0;
	center = Vector3d(0, 0, 0);
	fstream outFile("centers_a0-3.txt", ios_base::out | ios_base::app);
	outFile << "#" << helio->helio_index << endl;
	char argv[] = "qpa0.3zD";
	vector<Vector2d> centers;
	outFile << "#v" << endl;
	for (int i = 0; i < solution.size(); i++) {
		int n = solution[i].size();
		triangulateio in, out;
		in.numberofpoints = n;
		in.numberofpointattributes = 0;
		in.pointmarkerlist = (int*)NULL;
		in.pointlist = (REAL*)malloc(in.numberofpoints * 2 * sizeof(REAL));
		in.numberofsegments = n;
		in.segmentlist = (int*)malloc(n * 2 * sizeof(int));
		in.segmentmarkerlist = (int*)malloc(n * sizeof(int));
		for (int j = 0; j < n; j++) {
			in.pointlist[2 * j] = solution[i][j].X / (double)VERTEXSCALE;
			in.pointlist[2 * j + 1] = solution[i][j].Y / (double)VERTEXSCALE;
			in.segmentlist[2 * j] = j;
			in.segmentlist[2 * j + 1] = (j + 1) % n;
			in.segmentmarkerlist[j] = 1;
		}
		
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = NULL;
		out.pointlist = (REAL*)NULL;
		out.segmentlist = (int*)NULL;
		out.segmentmarkerlist = (int*)NULL;
		out.pointmarkerlist = (int*)NULL;
		out.trianglelist = (int*)NULL;
		triangulate(argv, &in, &out, NULL);

		for (int k = 0; k < out.numberoftriangles; ++k) {
			vector<float2> vertexes;
			Vector3d c(0, 0, 0);
			for (int j = 0; j < out.numberofcorners; ++j) {
				int idx = out.trianglelist[k*out.numberofcorners + j];
				vertexes.push_back(make_float2(out.pointlist[2 * idx], out.pointlist[2 * idx + 1]));
				c += Vector3d(vertexes.back().x, 0, vertexes.back().y);
				outFile << vertexes.back().x << ' ' << vertexes.back().y << ' ';
			}
			double tS = calcArea(vertexes);
			center += tS*c / 3.;
			centers.push_back(Vector2d(tS* c.x() / 3., tS*c.z() / 3.));
			S += tS;
			outFile << '\n';
		}

		free(in.pointlist);
		free(in.segmentlist);
		free(in.segmentmarkerlist);
		free(out.pointlist);
		free(out.segmentlist);
		free(out.segmentmarkerlist);
		free(out.pointmarkerlist);
		free(out.trianglelist);
	}

	vector<vector<Vector3d>> v_list(solution.size());
	for (int i = 0; i < solution.size(); ++i) {
		int n = solution[i].size();
		for (int j = 0; j < n; j++) {
			v_list[i].push_back(Vector3d(solution[i][j].X / (double)VERTEXSCALE, 0, solution[i][j].Y / (double)VERTEXSCALE));
		}
	}

	//Vector3d start_v(-helio->helio_size.x() / 2., 0, -helio->helio_size.z() / 2.);
	//double discrete = 0.01;
	//int row = helio->helio_size.x() / discrete;
	//int col = helio->helio_size.z() / discrete;
	//for (int k = 0; k < solution.size(); ++k) {
	//	for (int i = 0; i < row; ++i) {
	//		for (int j = 0; j < col; ++j) {
	//			Vector3d v = start_v + Vector3d(i*discrete, 0, j*discrete);
	//			if (GeometryFunc::inProjArea(v_list[k], v)) {
	//				centers.push_back(Vector2d(v.x(), v.z()));
	//			}
	//		}
	//	}
	//}
	
	
 	if (calcCenter)
		center /= S;
	outFile << "#c\n" << center.x() << ' ' << center.z() << endl;

	for (int i = 0; i < centers.size(); ++i) {
		centers[i] /= S;
	}
	vector<Vector2d> eigen_vecs;
	vector<double> eigen_val;
	double angle = PCAHandler::calcAxis(centers, eigen_vecs, eigen_val);
	outFile << "#e" << endl;
	outFile << eigen_vecs[0].x() << ' ' << eigen_vecs[0].y() << ' ' << eigen_vecs[1].x() << ' ' << eigen_vecs[1].y() << endl;
	outFile << eigen_val[0] << ' ' << eigen_val[1] << endl;
	outFile.close();

	return S;
}
