#pragma once
#include "../Common/CommonFunc.h"
#include "../DataStructure/Heliostat.h"
#include "../DataStructure/Receiver.h"

class SigmaFitting {
public:
	void fromFluxPeak2Sigma(const string peak_file, vector<Heliostat*>& helios, Receiver* recvs, const double DNI);

	vector<Vector3d> gt_data;			// x, y, flux_peak_intensity
public:
	void readFluxPeakFile(const string peak_file);
	void fitSigma(double flux_peak, Heliostat* helio, Receiver* recvs, const double DNI);
	float calcSigma(vector<float>& sigma_list) {
		return sqrt(pow(sigma_list[0], 2) * (pow(sigma_list[1], 2) + pow(sigma_list[2], 2) + pow(sigma_list[3], 2) + pow(sigma_list[4], 2))) / sqrt(sigma_list[5]);
	}
};