#include "FieldSimulationController/FieldSimulationController.h"

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << "[ERROR] Wrong parameters!" << endl;
		cout << "Usage: [InputJson] [SimulationType: -o -hf -ff]" << endl;
		return -1;
	}
	string str = string(argv[2]);
	ControllerMode controller_mode = FieldOptimismMode;
	if (str == "-hf")
		controller_mode = HelioFluxSimulationMode;
	else if (str == "-ff")
		controller_mode = FieldFluxSimulationMode;
	else if (str == "-e")
		controller_mode = EnergyCalculateMode;

	BaseController* controller = ControllerCreator::getController(controller_mode);
	controller->handler(argc, argv);

	return 1;
}