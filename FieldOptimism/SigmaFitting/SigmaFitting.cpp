#include "SigmaFitting.h"

void SigmaFitting::fromFluxPeak2Sigma(const string peak_file, vector<Heliostat*>& helios, Receiver* recvs, const double DNI)
{
	readFluxPeakFile(peak_file);
	//sigma_list.resize(helios.size());

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
	//sqrt(pow(sigma_list[0], 2) * (pow(sigma_list[1], 2) + pow(sigma_list[2], 2) + pow(sigma_list[3], 2) + pow(sigma_list[4], 2))) / sqrt(sigma_list[5]);
	
	int fc_index = helio->focus_center_index;
	Vector3d fc_pos = recvs->focus_center[fc_index];
	vector<float> sig_list(6, 0);
	sig_list[0] = (fc_pos - helio->helio_pos).norm();		// dis
	sig_list[1] = SIGMA_SUN;								// sigma_sun
	sig_list[2] = 2 * SIGMA_S;								// sigma_bq
	sig_list[3] = 0;										// sigma_ast
	sig_list[4] = 0;										// sigma_t
	sig_list[5] = abs(helio->cos_phi[fc_index]);			// cos_rev

	double low_range = 0;
	double high_range = 5.0;
	double percesion = 0.001;

	double min_err = INT_MAX;
	//int fc_index = helio->focus_center_index;
	double helio_param = DNI * helio->flux_param * (1 - helio->sd_bk) *helio->cos_phi[fc_index];
	double sigma = low_range + (high_range - low_range) / 2.0;
	while (high_range - low_range > percesion) {
		double mid = low_range + (high_range - low_range) / 2.0;
		//sig_list[3] = mid;
		double dis = flux_peak - helio_param / pow(mid, 2);
		double tmp_sigma = mid;
		//double tmp_sigma = calcSigma(sig_list);
		//double dis = flux_peak - DNI*helio->flux_param*(1 - helio->sd_bk)*helio->cos_phi[fc_index] / pow(tmp_sigma, 2);

		if (abs(dis) < min_err) {
			min_err = abs(dis);
			sigma = tmp_sigma;
		}
		if (dis > Epsilon) high_range = mid;
		else low_range = mid;
		//if (dis > Epsilon) high_range = mid;
		//else low_range = mid;

	}
	

	helio->sigma = sigma;
	//sigma_list[helio->helio_index] = sigma;
}
