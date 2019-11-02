#include "FieldSimulationController/FieldSimulationController.h"

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << "[ERROR] Wrong parameters!" << endl;
		cout << "Usage: [InputJson] [SimulationType: -o -f]" << endl;
		return -1;
	}
	ControllerType controller_type = string(argv[2]) == "-f" ? FluxSimulationType : FieldOptimismType;
	BaseController* controller = ControllerCreator::getController(controller_type);
	controller->handler(argc, argv);


	return 1;
}