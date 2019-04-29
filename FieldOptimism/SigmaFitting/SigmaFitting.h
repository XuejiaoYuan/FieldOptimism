#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"

class SigmaFitting {
public:
	void fromFluxPeak2Sigma(const string peak_file, vector<Heliostat*>& helios, const double DNI);
	void fitSigmaSurface();

	vector<Vector3d> gt_data;			// x, y, flux_peak_intensity
	vector<double> sigma_list;			// flux_peak to sigma
public:
	void readFluxPeakFile(const string peak_file);
	void fitSigma(double flux_peak, Heliostat* helio, const double DNI);
};