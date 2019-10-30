#include "DifferentialEvolution.h"


json DifferentialEvolution::fieldOptimism(ArgumentParser & argumentParser)
{
	// 1. Initialize populations
	initialization();
	
	// 2. Iteration
	for (int i = 0; i < argumentParser.getIterations(); ++i) {
		// 2.1 Mutation

		// 2.2 Crossover

		// 2.3 Selection
	}
	
	// 3. Get best field params

	return json();
}
