#include "DifferentialEvolution.h"
#include "../EnergyCalculatePipeline/EnergyCalculatePipeline.h"
#include "../Tool/Timer/Timer.h"

void DifferentialEvolution::initParams(json& de_config) {
	randomSeed = int(time(NULL));
	m_generator.seed(randomSeed);
	m_populationSize = de_config.get_with_default("populationSize", 50);
	m_population.resize(m_populationSize);
	m_maxCostPerAgent.resize(m_populationSize, 0);

	m_F = de_config.get_with_default("F", 0.8);
	m_CR = de_config.get_with_default("CR", 0.9);
	json constraints = de_config["constraints"].as<json>();
	for (string& param : m_paramKeys) {
		json con = constraints.get_with_default(param).as<json>();
		m_constraints[param] = Constraint(param, 
			con.get_with_default("lower", 0), con.get_with_default("upper", 1.0), con.get_with_default("is_constrained", false)
		);
	}
}

json DifferentialEvolution::fieldOptimism(ArgumentParser & _argumentParser)
{
	// 1. Preprocess 
	cout << "    2.1 Differential evolution initialize" << endl;
	argumentParser = &_argumentParser;
	json config = argumentParser->getConfig();
	initParams(config["DEParams"].as<json>());

	// 2. Initialize populations
	cout << "    2.2 Initialize population" << endl;
	initialization();
	
	// 3. Iteration
	cout << "    2.3 Start iteration" << endl;
	fstream out(argumentParser->getOutputPath() + "H" + to_string(argumentParser->getNumOfHelio()) + "de_process.txt", ios_base::out | ios_base::app);
	for (int i = 0; i < config["DEParams"]["Iterations"].as<int>(); ++i) {
		double old_maxCost = m_maxCostPerAgent[m_bestAgentIndex];
		out << "E: " << setprecision(12) << old_maxCost << "\nParams: " << m_population[m_bestAgentIndex] << endl;
		DEPipeline();
		double new_maxCost = m_maxCostPerAgent[m_bestAgentIndex];
		
	}
	out << m_maxCostPerAgent[m_bestAgentIndex] << endl;
	out.close();
	
	return m_population[m_bestAgentIndex];
}


void DifferentialEvolution::DEPipeline()
{
	uniform_int_distribution<int> distribution(0, m_populationSize - 1);

	for (int x = 0; x < m_populationSize; ++x) {
		json z = mutation(distribution, x);
		json newX = crossover(z, x);
		bool status = seletction(newX, x);
		if (!status) {
			--x;
			continue;
		}
		if (m_maxCostPerAgent[x] > m_maxCost) {
			m_maxCost = m_maxCostPerAgent[x];
			m_bestAgentIndex = x;
		}
	}
}

void DifferentialEvolution::initialization() {
	shared_ptr<uniform_real_distribution<double>> distribution;
	for (int i = 0; i < m_populationSize; ++i) {
		initializationHelper(distribution, i);
		cout << m_population[i] << endl;
		Timer::resetStart();
		EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, *argumentParser);
		m_maxCostPerAgent[i] = e_handler->handler(m_population[i]);
		Timer::printDuration("T");
		if (m_maxCostPerAgent[i] > m_maxCost) {
			m_maxCost = m_maxCostPerAgent[i];
			m_bestAgentIndex = i;
		}
		delete e_handler;
	}	
}

void DifferentialEvolution::initializationHelper(shared_ptr<uniform_real_distribution<double>>& distribution, const int idx)
{
	for (auto& constraint : m_constraints) {
		distribution = make_shared<uniform_real_distribution<double>>(uniform_real_distribution<double>(constraint.second.lower, constraint.second.upper));
		m_population[idx][constraint.second.key] = (*distribution)(m_generator);
	}
}

json DifferentialEvolution::mutation(uniform_int_distribution<int>& distribution, const int idx) {
	vector<int> ids = getRandomIndex(distribution, idx);

	// 变异操作
	json z;
	for (auto& param : m_paramKeys) {
		z[param] = m_population[ids[0]][param].as<double>() + m_F*(m_population[ids[1]][param].as<double>() - m_population[ids[2]][param].as<double>());
	}
	return z;
}

json DifferentialEvolution::crossover(json& z, const int idx) {
	// 3. 设置随机变量
	uniform_real_distribution<double> distributionParam(0, m_paramKeys.size());
	int R = distributionParam(m_generator);

	vector<double> r(m_paramKeys.size());
	uniform_real_distribution<double> distributionPerX(0, 1);
	for (auto& var : r)
		var = distributionPerX(m_generator);

	// 4. 交叉操作
	json newX;
	for (int i = 0; i < m_paramKeys.size(); ++i) {
		string key = m_paramKeys[i];
		if (r[i] < m_CR || i == R) 
			newX[key] = z[key].as<double>();
		else
			newX[key] = m_population[idx][key].as<json>();
	}
	return newX;
}

bool DifferentialEvolution::seletction(json& newX, const int idx) {
	// 1. 检查是否超过约束条件
	bool status = constraintCheck(newX);
	if (!status) {
		return false;
	}

	// 2. 计算新的cost并选择是否保留newX
	EnergyCalculatePipeline *e_handler = NULL;
	try
	{
		Timer::resetStart();
		e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, *argumentParser);
		double new_cost = e_handler->handler(newX);
		Timer::printDuration("T");
		if (new_cost > m_maxCostPerAgent[idx]) {
			m_population[idx] = newX;
			m_maxCostPerAgent[idx] = new_cost;
		}
	}
	catch (const std::exception&)
	{
		if(!e_handler) delete e_handler;
		return false;
	}

	return true;
}

bool DifferentialEvolution::constraintCheck(json& newX) {
	for (auto& key : m_paramKeys) {
		if (!m_constraints[key].check(newX[key].as<double>()))
			return false;
	}
	return true;
}

vector<int> DifferentialEvolution::getRandomIndex(uniform_int_distribution<int>& distribution, const int idx)
{
	// 1. 随机选取种群中的不同个体
	vector<int> ids(3, idx);
	while (ids[0] == idx || ids[1] == idx || ids[2] == idx || ids[0] == ids[1] || ids[0] == ids[2] || ids[1] == ids[2]) {
		ids[0] = distribution(m_generator);
		ids[1] = distribution(m_generator);
		ids[2] = distribution(m_generator);
	}
	return ids;
}

void FermatDifferentialEvolution::initializationHelper(shared_ptr<uniform_real_distribution<double>>& distribution, const int idx) {
	for (auto& constraint : m_constraints) {
		distribution = make_shared<uniform_real_distribution<double>>(uniform_real_distribution<double>(constraint.second.lower, constraint.second.upper));
		if (constraint.first == "gap")
			m_population[idx][constraint.first] = vector<double>({ (*distribution)(m_generator) });
		else
			m_population[idx][constraint.first] = (*distribution)(m_generator);
	}
}

json FermatDifferentialEvolution::mutation(uniform_int_distribution<int>& distribution, const int idx) {
	vector<int> ids = getRandomIndex(distribution, idx);

	// 变异操作
	json z;								
	for (auto& param : m_paramKeys) {
		if (param == "gap") {
			vector<double> a = m_population[ids[0]][param].as<vector<double>>();
			vector<double> b = m_population[ids[1]][param].as<vector<double>>();
			vector<double> c = m_population[ids[2]][param].as<vector<double>>();
			int maxIdx = 0;
			double maxCost = m_maxCostPerAgent[ids[0]];
			for (int i = 1; i < 3; ++i) {
				if (maxCost < m_maxCostPerAgent[ids[i]]) {
					maxCost = m_maxCostPerAgent[ids[i]];
					maxIdx = i;
				}
			}

			uniform_real_distribution<double> distributionParam(m_constraints[param].lower, m_constraints[param].upper);
			vector<double> d = m_population[ids[maxIdx]][param].as<vector<double>>();
			double v1, v2, v3;
			for (int i = 0; i < d.size(); ++i) {
				v1 = i < a.size() ? a[i] : distribution(m_generator);
				v2 = i < b.size() ? b[i] : distribution(m_generator);
				v3 = i < c.size() ? c[i] : distribution(m_generator);
				d[i] = v1 + m_F*(v2 - v3);
			}
			z[param] = d;
		}
		else{
			z[param] = m_population[ids[0]][param].as<double>() + m_F*(m_population[ids[1]][param].as<double>() - m_population[ids[2]][param].as<double>());
		}
	}
	return z;
}

json FermatDifferentialEvolution::crossover(json& z, const int idx) {
	// 3. 设置随机变量
	uniform_real_distribution<double> distributionParam(0, m_paramKeys.size());
	int R = distributionParam(m_generator);

	vector<double> r(m_paramKeys.size());
	uniform_real_distribution<double> distributionPerX(0, 1);
	for (auto& var : r)
		var = distributionPerX(m_generator);

	json newX;
	for (int i = 0; i < m_paramKeys.size(); ++i) {
		string key = m_paramKeys[i];
		if (key == "gap") {
			int sz_z = z[key].as<vector<double>>().size();
			int sz_m = m_population[idx][key].as<vector<double>>().size();
			int sz = max(sz_z, sz_m);
			vector<double> newGap;
			for (int j = 0; j < sz; ++j) {
				double r = distributionPerX(m_generator);
				double a = j < sz_z ? z[key].as<vector<double>>()[j] : m_population[idx][key].as<vector<double>>()[j];
				double b = j < sz_m ? m_population[idx][key].as<vector<double>>()[j] : z[key].as<vector<double>>()[j];
				if (r < m_CR || i == R)
					newGap.push_back(a);
				else
					newGap.push_back(b);
			}
			newX[key] = newGap;
		}
		else {
			if (r[i] < m_CR || i == R)
				newX[key] = z[key].as<double>();
			else
				newX[key] = m_population[idx][key].as<json>();
		}
	}
	return newX;
}

bool FermatDifferentialEvolution::constraintCheck(json& newX) {
	for (auto& key : m_paramKeys) {
		if (key == "gap") {
			for (auto& val : newX[key].as<vector<double>>())
				if (!m_constraints[key].check(val)) {
					//cout << key << ": " << val << endl;
					return false;
				}
		}
		else if (!m_constraints[key].check(newX[key].as<double>())) {
			//cout << key << ": " << newX[key].as<double>() << endl;
			return false;
		}
	}
	return true;
}