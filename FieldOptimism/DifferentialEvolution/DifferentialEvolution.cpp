#include "DifferentialEvolution.h"
#include "../EnergyCalculatePipeline/EnergyCalculatePipeline.h"
#include "../Tool/Timer/Timer.h"

// 
// [差分进化] 初始化问题参数
// 设置参数的上下界线，种群数量，F，CR参数
//
void DifferentialEvolution::initParams(json& de_config) {
	randomSeed = int(time(NULL));
	m_generator.seed(randomSeed);
	m_populationSize = de_config.get_with_default("populationSize", 10);
	m_population.resize(m_populationSize);
	m_maxCostPerAgent.resize(m_populationSize, 0);

	m_F = de_config.get_with_default("F", 0.8);
	m_CR = de_config.get_with_default("CR", 0.9);
	H = de_config.get_with_default("H", 10);

	adaptive = de_config.get_with_default("adaptive", true);

	M_CR.resize(H, 0.5);
	M_F.resize(H, 0.5);

	k = 0;

	json constraints = de_config["constraints"].as<json>();
	for (string& param : m_paramKeys) {
		json con = constraints.get_with_default(param).as<json>();

		bool is_constrained = con.get_with_default("is_constrained", false);
		if (is_constrained) {
			double lower = con.get_with_default("lower").as<double>();
			double upper = con.get_with_default("upper").as<double>();
			m_constraints[param] = Constraint(param, lower, upper, is_constrained);
		}
		else
			m_constraints[param] = Constraint(param, 0, 1, false);
		
	}
}

//
// [差分进化] 镜场优化主函数
//
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
	int cnt = 0;
	try {
		string str = "ad";
		if (!adaptive)
			str = "cl";
		fstream out(argumentParser->getOutputPath() + str+ "_small_H" + to_string(argumentParser->getNumOfHelio()) + "_L" + to_string(argumentParser->getLayoutType()) + "_de_process.txt", ios_base::out | ios_base::app);
		for (int i = 0; i < config["DEParams"]["Iterations"].as<int>(); ++i) {
			double old_maxCost = m_maxCostPerAgent[m_bestAgentIndex];
			out << "E: " << setprecision(12) << old_maxCost << "\nParams: " << m_population[m_bestAgentIndex] << endl;
			bool status = DEPipeline();
			double new_maxCost = m_maxCostPerAgent[m_bestAgentIndex];
			if (status) {
				++cnt;
				if (cnt == 10)
					break;
			}
		}
		out << m_maxCostPerAgent[m_bestAgentIndex] << endl;
		out.close();
	}
	catch(exception e){
		cerr << "[Warning] bad alloc!!!" << endl;
		return m_population[m_bestAgentIndex];
	}
	
	return m_population[m_bestAgentIndex];
}


//
// [差分进化] 单次差分进化操作，包括变异、交叉、选择
//
bool DifferentialEvolution::DEPipeline()
{
	uniform_int_distribution<int> distribution(0, m_populationSize - 1);
	S_CR.clear();
	S_F.clear();
	w_dis.clear();

	for (int x = 0; x < m_populationSize; ++x) {
		if(adaptive)
			getCRandF();
		json z = mutation(distribution, x);

		json newX = crossover(z, x);
		bool status = seletction(newX, x);
		if (!status) {
			--x;
			continue;
		}
		cout << x << ' ';
		if (m_maxCostPerAgent[x] > m_maxCost) {
			m_maxCost = m_maxCostPerAgent[x];
			m_bestAgentIndex = x;
		}
	}

	bool status = updateMCRandMF();
	return status;
}

//
// [差分进化] 初始化操作
//
void DifferentialEvolution::initialization() {
	shared_ptr<uniform_real_distribution<double>> distribution;
	for (int i = 0; i < m_populationSize; ++i) {
		cout << i << ' ';
		initializationHelper(distribution, i);
		EnergyCalculatePipeline* e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, *argumentParser);
		
		m_maxCostPerAgent[i] = e_handler->handler(m_population[i]);
		if (m_maxCostPerAgent[i] < 1) {
			--i;
			continue;
		}
		if (m_maxCostPerAgent[i] > m_maxCost) {
			m_maxCost = m_maxCostPerAgent[i];
			m_bestAgentIndex = i;
		}
		delete e_handler;
	}	
}

//
// [差分进化] 初始化操作核心函数
//
void DifferentialEvolution::initializationHelper(shared_ptr<uniform_real_distribution<double>>& distribution, const int idx)
{
	for (auto& constraint : m_constraints) {
		distribution = make_shared<uniform_real_distribution<double>>(uniform_real_distribution<double>(constraint.second.lower, constraint.second.upper));
		m_population[idx][constraint.second.key] = (*distribution)(m_generator);
	}
}


//
// [差分进化] 变异操作
//
json DifferentialEvolution::mutation(uniform_int_distribution<int>& distribution, const int idx) {
	vector<int> ids = getRandomIndex(distribution, idx);

	json z;
	for (auto& param : m_paramKeys) {
		if(!adaptive)
			z[param] = m_population[ids[0]][param].as<double>() + m_F*(m_population[ids[1]][param].as<double>() - m_population[ids[2]][param].as<double>());
		else {
			z[param] = m_population[idx][param].as<double>() + m_F*(m_population[m_bestAgentIndex][param].as<double>() - m_population[idx][param].as<double>()) +
				m_F *(m_population[ids[0]][param].as<double>() - m_population[ids[1]][param].as<double>());
			if (z[param] < m_constraints[param].lower)
				z[param] = (m_constraints[param].lower + m_population[idx][param].as<double>()) / 2.;
			else if (z[param] > m_constraints[param].upper)
				z[param] = (m_constraints[param].upper + m_population[idx][param].as<double>()) / 2.;
		}
	}
	return z;
}

// 
// [差分进化] 交叉操作
//
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

//
// [差分进化] 选择操作
//
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
		e_handler = EnergyCalculateCreator::getPipeline(SceneEnergyMode, *argumentParser);
		double new_cost = e_handler->handler(newX);
		if (new_cost > m_maxCostPerAgent[idx]) {
			w_dis.push_back(abs(new_cost - m_maxCostPerAgent[idx]));
			S_CR.push_back(m_CR);
			S_F.push_back(m_F);

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

//
// [差分进化] 判断生成的参数是否超过界线
//
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
	vector<int> ids;
	if (!adaptive) {
		ids.resize(3);
		while (ids[0] == idx || ids[1] == idx || ids[2] == idx || ids[0] == ids[1] || ids[0] == ids[2] || ids[1] == ids[2]) {
			ids[0] = distribution(m_generator);
			ids[1] = distribution(m_generator);
			ids[2] = distribution(m_generator);
		}
	}
	else {
		ids.resize(2);
		while (ids[0] == idx || ids[1] == idx || ids[0] == ids[1] || ids[0] == m_bestAgentIndex || ids[1] == m_bestAgentIndex) {
			ids[0] = distribution(m_generator);
			ids[1] = distribution(m_generator);
		}
	}
	
	return ids;
}

//
// [差分进化] 更新CR和F参数的辅助数组
//
vector<int> DifferentialEvolution::getCRandF()
{
	uniform_int_distribution<int> distribution(0, H - 1);
	int r = distribution(m_generator);
	normal_distribution<double> n_distribution(M_CR[r], 0.1); 
	cauchy_distribution<double> c_distribution(M_F[r], 0.1);
	uniform_real_distribution<double> p_distribution(2. / m_populationSize, 0.2);

	m_CR = n_distribution(m_generator);
	if (m_CR < 0)
		m_CR = 0;
	else if (m_CR > 1) 
		m_CR = 1;
	
	do {
		m_F = min(1.0, c_distribution(m_generator));
	} while (m_F <= 0);
	
	m_p = p_distribution(m_generator);

	return vector<int>();
}

//
// [差分进化] 更新CR和F参数
//
bool DifferentialEvolution::updateMCRandMF()
{
	if (w_dis.empty())
		return false;
	vector<double> w;
	double sum = 0;
	for (int i = 0; i < w_dis.size(); ++i)
		sum += w_dis[i];
	cout << sum << endl;
	double mean_cr = 0;
	double up = 0, down = 0;
	for (int i = 0; i < w_dis.size(); ++i) {
		mean_cr += w_dis[i] / sum*S_CR[i];
		up += w_dis[i] * S_F[i] * S_F[i];
		down += w_dis[i] * S_F[i];
	}
	M_CR[k] = mean_cr;
	M_F[k] = up / down;
	k = (k + 1) % H;
	return sum < 0.5;
}

void RadialStaggerDifferentialEvolution::initializationHelper(shared_ptr<uniform_real_distribution<double>>& distribution, const int idx) {
	for (auto& constraint : m_constraints) {
		distribution = make_shared<uniform_real_distribution<double>>(uniform_real_distribution<double>(constraint.second.lower, constraint.second.upper));
		if (adaptiveKey(constraint.first))
			m_population[idx][constraint.first] = vector<double>({ (*distribution)(m_generator) });
		else
			m_population[idx][constraint.first] = (*distribution)(m_generator);
	}
}

json RadialStaggerDifferentialEvolution::mutation(uniform_int_distribution<int>& distribution, const int idx) {
	vector<int> ids = getRandomIndex(distribution, idx);

	// 变异操作
	json z;								
	for (auto& param : m_paramKeys) {
		if (adaptiveKey(param)) {
			uniform_real_distribution<double> distributionParam(m_constraints[param].lower, m_constraints[param].upper);
			double v1, v2, v3, v4;
			vector<double> a, b, c, best;
			if (!adaptive) {
				a = m_population[ids[0]][param].as<vector<double>>();
				b = m_population[ids[1]][param].as<vector<double>>();
				c = m_population[ids[2]][param].as<vector<double>>();
				best = m_population[m_bestAgentIndex][param].as<vector<double>>();
			}
			else {
				a = m_population[idx][param].as<vector<double>>();
				best = m_population[m_bestAgentIndex][param].as<vector<double>>();
				b = m_population[ids[0]][param].as < vector<double>>();
				c = m_population[ids[1]][param].as < vector<double>>();
			}
			
			for (int i = 0; i < best.size(); ++i) {
				if (!adaptive) {
					v1 = i < a.size() ? a[i] : distribution(m_generator);
					v2 = i < b.size() ? b[i] : distribution(m_generator);
					v3 = i < c.size() ? c[i] : distribution(m_generator);
					best[i] = v1 + m_F*(v2 - v3);
				}
				else {
					v1 = i < a.size() ? a[i] : distribution(m_generator);
					v2 = i < best.size() ? best[i] : distribution(m_generator);
					v3 = i < b.size() ? b[i] : distribution(m_generator);
					v4 = i < c.size() ? c[i] : distribution(m_generator);
					best[i] = v1 + m_F*(v2 - v1 + v3 - v4);
				}
				// 随机选择进行修改
			}
			z[param] = best;
		}
		else{
			if(!adaptive)
				z[param] = m_population[ids[0]][param].as<double>() + m_F*(m_population[ids[1]][param].as<double>() - m_population[ids[2]][param].as<double>());
			else {
				z[param] = m_population[idx][param].as<double>() + m_F*(m_population[m_bestAgentIndex][param].as<double>() - m_population[idx][param].as<double>()) +
					m_F *(m_population[ids[0]][param].as<double>() - m_population[ids[1]][param].as<double>());
				if (z[param] < m_constraints[param].lower)
					z[param] = (m_constraints[param].lower + m_population[idx][param].as<double>()) / 2.;
				else if (z[param] > m_constraints[param].upper)
					z[param] = (m_constraints[param].upper + m_population[idx][param].as<double>()) / 2.;
			}
		}
	}
	return z;
}

json RadialStaggerDifferentialEvolution::crossover(json& z, const int idx) {
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
		if (adaptiveKey(key)) {
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

bool RadialStaggerDifferentialEvolution::constraintCheck(json& newX) {
	for (auto& key : m_paramKeys) {
		if (adaptiveKey(key)) {
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