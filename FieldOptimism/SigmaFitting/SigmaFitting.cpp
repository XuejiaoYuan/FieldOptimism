#include "SigmaFitting.h"

void SigmaFitting::fromFluxPeak2Sigma(const string peak_file, vector<Heliostat*>& helios, const double DNI)
{
	readFluxPeakFile(peak_file);
	sigma_list.resize(helios.size());

#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i) {
		fitSigma(gt_data[i].z(), helios[i], DNI);
	}
}

void SigmaFitting::fitSigmaSurface()
{
}

void SigmaFitting::readFluxPeakFile(const string peak_file)
{
	fstream inFile(peak_file);;
	string line;
	while (getline(inFile, line)) {
		stringstream ss(line);
		Vector3d data;
		ss >> data.x() >> data.y() >> data.z();
		gt_data.push_back(data);
	}
	inFile.close();
}

void SigmaFitting::fitSigma(double flux_peak, Heliostat* helio, const double DNI)
{
	double low_range = 0;
	double high_range = 5;
	double percesion = 0.01;

	double min_err = INT_MAX;
	int fc_index = helio->focus_center_index;
	double helio_param = DNI * helio->flux_param * (1 - helio->sd_bk) * helio->cos_phi[fc_index];
	double sigma = low_range + (high_range - low_range) / 2.0;
	while (high_range - low_range > percesion) {
		double mid = low_range + (high_range - low_range) / 2.0;
		double dis = flux_peak - helio_param / pow(mid, 2);
	
		if (abs(dis) < min_err) {
			min_err = abs(dis);
			sigma = mid;
		}
		if (dis > Epsilon) high_range = mid;
		else low_range = mid;
	}
	

	helio->sigma = sigma;
	sigma_list[helio->helio_index] = sigma;
}
