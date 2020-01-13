#include "FieldSimulationController.h"

void EnergyCalculateController::handler(int argc, char** argv) {
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);
	
	// 2. Calculate filed energy
	EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, argumenParser);
	double res = e_handler->handler();
	cout << "\tTotal energy: " << res << endl;
}

void FieldOptimismController::handler(int argc, char ** argv)
{
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumentParser;
	argumentParser.parser(argc, argv);

	// 2. Start field optimization
	cout << "2. Start field optimism" << endl;
	DifferentialEvolution* de_handler = DECreator::getDE(argumentParser.getLayoutType());
	json field_args = de_handler->fieldOptimism(argumentParser);

		// 3. Save field result
	cout << "3. Saving result" << endl;
	EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, argumentParser);
	double new_cost = e_handler->handler(field_args);
	e_handler->saveSolarScene();
}

void HelioFluxSimulationController::handler(int argc, char** argv) {
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);

	// 2. Calculate receiver flux distribution
	EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(HelioFluxDensityMode, argumenParser);
	e_handler->handler(argumenParser.getConfig()["FieldArgs"].as<json>());
}

void FieldFluxSimulationController::handler(int argc, char ** argv)
{
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);

	// 2. Calculate receiver flux distribution
	EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(FieldFluxDensityMode, argumenParser);
	e_handler->handler(argumenParser.getConfig()["FieldArgs"].as<json>());
}
