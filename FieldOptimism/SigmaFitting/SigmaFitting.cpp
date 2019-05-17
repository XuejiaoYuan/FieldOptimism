#include "SigmaFitting.h"

void SigmaFitting::fromFluxPeak2Sigma(const string peak_file, vector<Heliostat*>& helios, Receiver* recvs, const double DNI)
{
	readFluxPeakFile(peak_file);

#pragma omp parallel for
	for (int i = 0; i < helios.size(); ++i) {
		fitSigma(gt_data[i].z(), helios[i], recvs, DNI);
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

void SigmaFitting::fitSigma(double flux_peak, Heliostat* helio, Receiver* recvs, const double DNI)
{
	int fc_index = helio->focus_center_index;
	Vector3d fc_pos = helio->focus_center;

	double low_range = 0;
	double high_range = 5.0;
	double percesion = 0.001;

	double min_err = INT_MAX;
	double helio_param = DNI * helio->flux_param * (1 - helio->sd_bk) *helio->cos_phi[fc_index];
	double sigma = low_range + (high_range - low_range) / 2.0;
	while (high_range - low_range > percesion) {
		double tmp_sigma = low_range + (high_range - low_range) / 2.0;
		double dis = flux_peak - helio_param / pow(tmp_sigma, 2);

		if (abs(dis) < min_err) {
			min_err = abs(dis);
			sigma = tmp_sigma;
		}
		if (dis > Epsilon) high_range = tmp_sigma;
		else low_range = tmp_sigma;
	}
	

	helio->sigma = sigma;
}
