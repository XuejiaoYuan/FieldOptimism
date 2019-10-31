#include "DifferentialEvolution.h"


json DifferentialEvolution::fieldOptimism(ArgumentParser & argumentParser)
{
	// 1. Initialize populations
	initialization();
	
	// 2. Iteration
	json config = argumentParser.getConfig();
	for (int i = 0; i < config["Iterations"].as<int>(); ++i) {
		// 2.1 Mutation

		// 2.2 Crossover

		// 2.3 Selection
	}
	
	// 3. Get best field params

	return json();
}
