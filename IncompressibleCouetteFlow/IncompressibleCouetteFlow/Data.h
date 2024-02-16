#pragma once
#include <vector>

class Data {
public:
	std::vector<double> x_u;
	std::vector<double> y_u;
	std::vector<double> x_v;
	std::vector<double> y_v;
	std::vector<double> x_p;
	std::vector<double> y_p;
	std::vector<double> t;
	std::vector<std::vector<std::vector<double>>> u;
	std::vector<std::vector<std::vector<double>>> v;
	std::vector<std::vector<std::vector<double>>> p;
};
