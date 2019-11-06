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
		if (h_helio_energy[i] < 0) {
			cout << "pause" << endl;
		}
	}
	delete[] h_helio_energy;
}
