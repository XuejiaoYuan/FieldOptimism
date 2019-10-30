#pragma once
#include "ReceiverDeviceArgument.h"

class CylinderReceiverDeviceArgument : public ReceiverDeviceArgument {
	void setRecvDeviceArguments(Receiver& recv, vector<Heliostat*>& helios);
};