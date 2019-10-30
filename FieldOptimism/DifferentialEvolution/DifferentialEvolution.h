#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"


class DifferentialEvolution
{
public:
	json fieldOptimism(ArgumentParser& argumentParser);

protected:
	vector<json> fields;
	virtual void initialization() {};
	virtual void mutation() {};
	virtual void crossover() {};
	virtual void seletction() {};
};


class RectDifferentialEvolution : public DifferentialEvolution{
public:
	virtual void initialization() {};
	virtual void mutation() {};
	virtual void crossover() {};
	virtual void seletction() {};
};

class CrossRectDifferentialEvolution :public DifferentialEvolution {
public:
	virtual void initialization() {};
	virtual void mutation() {};
	virtual void crossover() {};
	virtual void seletction() {};
};

class FermatDifferentialEvolution:public DifferentialEvolution {
public:
	virtual void initialization() {};
	virtual void mutation() {};
	virtual void crossover() {};
	virtual void seletction() {};
};

class RadialDifferentialEvolution:public DifferentialEvolution {
public:
	virtual void initialization() {};
	virtual void mutation() {};
	virtual void crossover() {};
	virtual void seletction() {};
};

class DECreator {
public:
	static DifferentialEvolution* getDE(LayoutType layout_type) {
		switch (layout_type)
		{
		case RectLayoutType:
			return new RectDifferentialEvolution();
		case CrossRectLayoutType:
			return new CrossRectDifferentialEvolution();
		case FermatLayoutType:
			return new FermatDifferentialEvolution();
		case RadialLayoutType:
			return new RadialDifferentialEvolution();
		default:
			return new RectDifferentialEvolution();
		}
	}
};