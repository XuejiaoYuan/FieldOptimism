#include "FieldOptimismPipeline.h"
#include "../Tool/ArgumentParser/ArgumentParser.h"

void FieldOptimismPipeline::fieldOptimism(int argc, char* argv[]) {
	// 1. Initialize parameters
	cout << "1. Start loading arguments" << endl;
	ArgumentParser argumenParser;
	argumenParser.parser(argc, argv);
	
	// 2. Start field optimization
	DifferentialEvolution* de_handler = DECreator::getDE(argumenParser.getLayoutType());
	json field_args = de_handler->fieldOptimism(argumenParser);

	// 3. Save field result
	fstream out(argumenParser.getOutputFn());
	out << field_args;
}
