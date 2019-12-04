#pragma once
#include "../../DataStructure/SolarScene.h"
#include "../../DataStructure/Receiver.h"
#include "../../GaussLegendre/GaussLegendre.cuh"
#include "../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../DeviceArgument/ReceiverDeviceArgument/CylinderReceiverDeviceArgument.h"
#include "RectangleReceiverFluxIntegral/RectRecvFluxIntegral.cuh"
#include "CylinderReceiverFluxIntegral/CylinderRecvFluxIntegral.cuh"
#include "RectangleReceiverDiscreteFlux/RectRecvDiscreteFlux.cuh"

class ReceiverEnergyCalculator
{
public:
	ReceiverEnergyCalculator(SolarScene *&solar_scene, GaussLegendre* gl_handler, int m, int n, bool calcCenterMode = false) :
		m(m), n(n), solar_scene(solar_scene), gl_handler(gl_handler), calcCenterMode(calcCenterMode), 
		h_helio_energy(nullptr), d_helio_energy(nullptr), h_recv_flux(nullptr), d_recv_flux(nullptr) {}
	void calcRecvEnergySum();
	vector<float> calcRecvFluxDistribution();
	~ReceiverEnergyCalculator() {
		solar_scene = nullptr;
		gl_handler = nullptr;

		h_args.clear();
		r_args->clear();
		//delete r_args;
		//delete[] h_helio_energy;
		cudaFree(d_helio_energy);
		cudaFree(d_recv_flux);
		h_helio_energy = nullptr;
		d_helio_energy = nullptr;
		h_recv_flux = nullptr;
		d_recv_flux = nullptr;
	}

private:
	int m, n;
	bool calcCenterMode;
	SolarScene* solar_scene;
	GaussLegendre* gl_handler;
	IntegralHelioDeviceArgumet h_args;
	ReceiverDeviceArgument* r_args;
	float *h_helio_energy, *d_helio_energy;
	float *h_recv_flux, *d_recv_flux;
	vector<float> helio_energy;
	void setDeviceHelioEnergy();
	void storeHeliostatEnergy();
	void setDeviceRecvFlux();
	vector<float> storeDeviceRecvFlux();
};

