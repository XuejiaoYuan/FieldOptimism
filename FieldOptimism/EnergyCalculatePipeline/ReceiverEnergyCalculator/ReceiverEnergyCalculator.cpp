#include "ReceiverEnergyCalculator.h"

void ReceiverEnergyCalculator::calcRecvEnergySum()
{
	ReceiverType recv_type = solar_scene->recvs[0]->recv_type;
	
	// 1. Set heliostat device arguments
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);
	h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
	if (calcCenterMode) h_args.setHelioCenterBias(solar_scene->helios);

	// 2. Calculate field total energy
	// 2.1 Set device heliostat energy
	int helioNum = solar_scene->helios.size();
	setDeviceHelioEnergy();
	
	// 2.2 Set receiver arguments
	// 2.3 Calculate flux density integral
	switch (recv_type)
	{
	case RectangularRecvType:
	case PolyhedronRecvType:
		r_args = new ReceiverDeviceArgument();
		r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);
		calcRectRecvEnergySum(m, n, helioNum, h_args, *r_args, *gl_handler, d_helio_energy);
		break;
	case CylinderRecvType:
		r_args = new CylinderReceiverDeviceArgument();
		r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);
		calcCylinderRecvEnergySum(m, n, helioNum, h_args, *r_args, *gl_handler, d_helio_energy);
		break;
	default:
		throw runtime_error("[Error RecvEnergyCalculator] Wrong receiver type!!!\n");
	}

	// 2.4 Store heliostat energy
	storeHeliostatEnergy();

}

vector<float> ReceiverEnergyCalculator::calcRecvFluxDistribution()
{
	ReceiverType recv_type = solar_scene->recvs[0]->recv_type;
	Vector2i row_col = solar_scene->recvs[0]->rows_cols;

	// 1. Set heliostat device arguments
	h_args.setHelioDevicePos(solar_scene->helios);
	h_args.setHelioDeviceArguments(solar_scene->helios);
	h_args.setHelioRecvArguments(solar_scene->helios, *(solar_scene->recvs[0]));
	if (calcCenterMode) h_args.setHelioCenterBias(solar_scene->helios);

	// 2. Calculate field total energy
	// 2.1 Set device heliostat energy
	int helioNum = solar_scene->helios.size();
	setDeviceRecvFlux();

	// 2.2 Set receiver arguments
	// 2.3 Calculate flux density integral
	switch (recv_type)
	{
	case RectangularRecvType:
	case PolyhedronRecvType:
		r_args = new ReceiverDeviceArgument();
		r_args->setRecvDeviceArguments((*solar_scene->recvs[0]), solar_scene->helios);
		calcRectRecvFluxDistribution(h_args, *r_args, *gl_handler, d_recv_flux, row_col.x(), row_col.y());
		break;
	case CylinderRecvType:
		throw runtime_error("[Error] Not implement!!!\n");
	default:
		throw runtime_error("[Error RecvEnergyCalculator] Wrong receiver type!!!\n");
	}


	// 2.4 Store heliostat energy
	return storeDeviceRecvFlux();
}

void ReceiverEnergyCalculator::setDeviceHelioEnergy() {
	int helioNum = solar_scene->helios.size();
	h_helio_energy = new float[helioNum];
	for (int i = 0; i < helioNum; ++i)
		h_helio_energy[i] = 0;
	cudaMalloc((void**)&d_helio_energy, sizeof(float)*helioNum);
	cudaMemcpy(d_helio_energy, h_helio_energy, sizeof(float)*helioNum, cudaMemcpyHostToDevice);
}

void ReceiverEnergyCalculator::storeHeliostatEnergy() {
	int helioNum = solar_scene->helios.size();
	cudaMemcpy(h_helio_energy, d_helio_energy, sizeof(float)*helioNum, cudaMemcpyDeviceToHost);
	
	for (int i = 0; i < helioNum; ++i) {
		solar_scene->helios[i]->energy += h_helio_energy[i] * solar_scene->DNI;
	}
	delete[] h_helio_energy;
}

void ReceiverEnergyCalculator::setDeviceRecvFlux()
{
	int row = solar_scene->recvs[0]->rows_cols.x();
	int col = solar_scene->recvs[0]->rows_cols.y();
	int recv_num = solar_scene->recvs[0]->recv_face_num;
	h_recv_flux = new float[recv_num*row*col];
	for (int i = 0; i < recv_num*row*col; ++i)
		h_recv_flux[i] = 0;
	cudaMalloc((void**)&d_recv_flux, sizeof(float)*recv_num*row*col);
	cudaMemcpy(d_recv_flux, h_recv_flux, sizeof(float)*recv_num*row*col, cudaMemcpyHostToDevice);
}

vector<float> ReceiverEnergyCalculator::storeDeviceRecvFlux() {
	vector<float> flux;
	int row = solar_scene->recvs[0]->rows_cols.x();
	int col = solar_scene->recvs[0]->rows_cols.y();
	int recv_num = solar_scene->recvs[0]->recv_face_num;
	cudaMemcpy(h_recv_flux, d_recv_flux, sizeof(float)*recv_num*row*col, cudaMemcpyDeviceToHost);

	for (int i = 0; i < recv_num*row*col; ++i)
		flux.push_back(h_recv_flux[i]);

	delete[] h_recv_flux;
	return flux;
}