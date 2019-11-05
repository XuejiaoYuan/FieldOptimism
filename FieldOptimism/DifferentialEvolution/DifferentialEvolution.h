#pragma once
#include "../Tool/ArgumentParser/ArgumentParser.h"
#include <random>
using namespace std;

class Constraint {
public:
	Constraint(string key="", double lower = 0.0, double upper = 0.0, bool is_constrained = false)
		:key(key), lower(lower), upper(upper), is_constrained(is_constrained) {};
	string key;
	double lower;
	double upper;
	bool is_constrained;
	bool check(double candidate) {
		if (is_constrained) {
			if (lower <= candidate && candidate <= upper)
				return true;
			else
				return false;
		}
		else 
			return true;
	}
};

class DifferentialEvolution
{
public:
	void initParams(json& constraints);
	json fieldOptimism(ArgumentParser& argumentParser);

protected:
	int randomSeed;
	int m_populationSize;
	int m_bestAgentIndex;
	double m_F;
	double m_CR;
	double m_maxCost;
	ArgumentParser* argumentParser;
	vector<string> m_paramKeys;
	unordered_map<string, Constraint> m_constraints;
	default_random_engine m_generator;
	vector<json> m_population;
	vector<double> m_maxCostPerAgent;

	void DEPipeline();
	virtual void initialization();
	virtual json mutation(uniform_int_distribution<int>& distribution, const int idx);
	virtual json crossover(json& z, const int idx);
	virtual bool seletction(json& newX, const int idx);
	vector<int> getRandomIndex(uniform_int_distribution<int>& distribution, const int idx);
};


class RectDifferentialEvolution : public DifferentialEvolution{
public:
	RectDifferentialEvolution() {
		m_paramKeys = { "z_start", "interval_ratio_x", "interval_ratio_z" };
	}
};

class CrossRectDifferentialEvolution :public RectDifferentialEvolution {
};

class FermatDifferentialEvolution:public DifferentialEvolution {
public:
	FermatDifferentialEvolution() {
		m_paramKeys = { "dsep", "gap", "helio_recv_dis" };
	}

protected:
	json mutation(uniform_int_distribution<int>& distribution, const int idx);
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