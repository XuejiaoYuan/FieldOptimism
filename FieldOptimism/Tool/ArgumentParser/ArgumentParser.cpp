#include "ArgumentParser.h"

bool ArgumentParser::parser(int argc, char** argv) {
	if (argc != 2) {
		cout << "The number of arguments is wrong!" << endl;
		return false;
	}

	config_path = string(argv[1]);

	fstream in_file(config_path, ios_base::in);
	in_file >> config;
	in_file.close();
	
	try
	{
		auto path = config.get_with_default("Path");
		input_path = path["InputPath"].as_string();
		output_path = path["OutputPath"].as_string();
		sunray_fn = path["SunrayFile"].as_string();
		output_fn = path["OutputFileName"].as_string();
		sunray.calcSunRay(sunray_fn);
		
		auto scene = config.get_with_default("Scene");
		layout_type = (LayoutType)(scene["LayoutType"].as<int>());
		num_of_helio = scene["NumOfHelio"].as<int>();
		model_type = (ModelType)(scene["ModelType"].as<int>());
		calc_sigma = scene["CalcSigma"].as<bool>();

		auto helio = config.get_with_default("Heliostat");
		heliostat = HeliostatCreator::getHeliostat((HelioType)helio["type"].as<int>());
		heliostat->initHelio(helio);

		auto recv = config.get_with_default("Receiver");
		Receiver *receiver = ReceiverCreator::getReceiver((ReceiverType)recv["type"].as<int>());		// Only handle one receiver
		receiver->initRecv(recv);
		recvs.push_back(receiver);

		getTimeParams();
		
		auto de = config.get_with_default("DE");
		iterations = de["iterations"].as<int>();
	}
	catch (const std::exception&)
	{
		cout << "Wrong configuration!" << endl;
		return false;
	}
	return true;
	
}

vector<int> ArgumentParser::getTimeParams(TimeType type) {
	switch (type) {
	case MonthType:
		return months;
	case DayType:
		return days;
	case HourType:
		return hours;
	default:
		return{};
	}
}

void ArgumentParser::getTimeParams() {
	auto t = config.get_with_default("Time");
	unordered_map<string, vector<int>&> hash({{"month", months}, {"day", days}, {"hour", hours}});
	for(auto& item: hash){
		string key = item.first;
		for (int m = t.get_with_default(key)["start"].as<int>(); m <= t.get_with_default(key)["end"].as<int>(); m += t.get_with_default(key)["gap"].as<int>())
			item.second.push_back(m);
	}
	minute_gap = t.get_with_default("minute")["gap"].as<int>();
}