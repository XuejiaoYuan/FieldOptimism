#include "PCA.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace cv;


double PCAHandler::calcAxis(const vector<Vector2d>& pts, vector<Vector2d>& eigen_vecs, vector<double>& eigen_val)
{
	int sz = static_cast<int>(pts.size());
	Mat data_pts = Mat(sz, 2, CV_64F);
	for (int i = 0; i < data_pts.rows; i++)
	{
		data_pts.at<double>(i, 0) = pts[i].x();
		data_pts.at<double>(i, 1) = pts[i].y();
	}
	//Perform PCA analysis
	PCA pca_analysis(data_pts, Mat(), PCA::DATA_AS_ROW);

	//Store the center of the object
	Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)),
		static_cast<int>(pca_analysis.mean.at<double>(0, 1)));

	//Store the eigenvalues and eigenvectors
	eigen_vecs.resize(2);
	eigen_val.resize(2);
	for (int i = 0; i < 2; i++)
	{
		eigen_vecs[i] = Vector2d(pca_analysis.eigenvectors.at<double>(i, 0),
			pca_analysis.eigenvectors.at<double>(i, 1));
		eigen_val[i] = pca_analysis.eigenvalues.at<double>(i);
	}
	// Draw the principal components
	//Point p1 = cntr + 0.02 * Point(static_cast<int>(eigen_vecs[0].x * eigen_val[0]), static_cast<int>(eigen_vecs[0].y * eigen_val[0]));
	//Point p2 = cntr - 0.02 * Point(static_cast<int>(eigen_vecs[1].x * eigen_val[1]), static_cast<int>(eigen_vecs[1].y * eigen_val[1]));
	double angle = atan2(eigen_vecs[0].y(), eigen_vecs[0].x()); // orientation in radians
	return angle;
}
