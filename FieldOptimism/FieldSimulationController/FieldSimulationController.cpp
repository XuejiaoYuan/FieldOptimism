#include "FieldSimulationController.h"

void FieldOptimismController::handler(int argc, char ** argv)
{
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);

	// 2. Start field optimization
	cout << "2. Start field optimism" << endl;
	DifferentialEvolution* de_handler = DECreator::getDE(argumenParser.getLayoutType());
	json field_args = de_handler->fieldOptimism(argumenParser);

	// 3. Save field result
	cout << "3. Saving result" << endl;
	fstream out(argumenParser.getOutputFn());
	out << field_args;
}

void FluxSimulationController::handler(int argc, char** argv) {
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);

	// 2. Calculate receiver flux distribution
	EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(FluxDensityMode, argumenParser);
	e_handler->handler(argumenParser.getConfig()["FieldArgs"].as<json>());
}
