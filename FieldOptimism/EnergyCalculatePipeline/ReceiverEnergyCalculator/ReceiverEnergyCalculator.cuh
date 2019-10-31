#pragma once
#include "../../DataStructure/SolarScene.h"
#include "../../DataStructure/Receiver.h"
#include "../../GaussLegendre/GaussLegendre.cuh"
#include "../../DeviceArgument/HeliostatDeviceArgument/HeliostatDeviceArgument.h"
#include "../../DeviceArgument/ReceiverDeviceArgument/ReceiverDeviceArgument.h"
#include "../../DeviceArgument/ReceiverDeviceArgument/CylinderReceiverDeviceArgument.h"
#include "RectangleReceiverFluxIntegral/RectRecvFluxIntegral.cuh"
#include "CylinderReceiverFluxIntegral/CylinderRecvFluxIntegral.cuh"

class ReceiverEnergyCalculator
{
public:
	ReceiverEnergyCalculator(SolarScene *solar_scene, GaussLegendre* gl_handler, int m, int n, bool calcCenterMode = false) :
		m(m), n(n), solar_scene(solar_scene), gl_handler(gl_handler), calcCenterMode(calcCenterMode) {}
	void calcRecvEnergySum(float* d_total_energy);
	~ReceiverEnergyCalculator() {
		solar_scene = nullptr;
		gl_handler = nullptr;
	}

private:
	int m, n;
	bool calcCenterMode;
	SolarScene* solar_scene;
	GaussLegendre* gl_handler;
	IntegralHelioDeviceArgumet h_args;
	ReceiverDeviceArgument* r_args;
	void calcRectRecvEnergySum(float* d_total_energy);
	void calcCylinderRecvEnergySum(float* d_total_energy);
};

