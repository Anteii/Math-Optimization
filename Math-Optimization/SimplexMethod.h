#pragma once
#include <vector>
#include <string>
#include <regex>
#include <iostream>
#include "CppConsoleTable.hpp"

struct Expression {
	std::vector<double> c;
	std::string op;
	double b;
};

struct SimplexTable {
	std::vector<std::string> basis;
	std::vector<std::string> vars;
	std::vector<std::vector<double>> rows;
	std::vector<double> delta;
	std::vector<double> b;
};

std::ostream& operator<< (std::ostream& os, const Expression& exp);
std::ostream& operator<< (std::ostream& os, const SimplexTable& exp);

std::vector<double> parseExpression(std::string input);
Expression parseLimitation(std::string input);
Expression standLimitation(Expression& exp);

std::vector<double> doSimplex(std::vector<double>&, std::vector<Expression>&);