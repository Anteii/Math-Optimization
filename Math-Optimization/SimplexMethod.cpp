#include "SimplexMethod.h"
#include <algorithm>
double eps = 1e-10;

std::vector<double> parseExpression(std::string input) {
	std::regex var_regex("[+-][\\d\\.]+x_\\d+");
	std::vector<std::string> vars;
	std::vector<double> c;

	auto input_begin =
		std::sregex_iterator(input.begin(), input.end(), var_regex);
	auto input_end = std::sregex_iterator();

	for (std::sregex_iterator i = input_begin; i != input_end; ++i) {
		std::smatch match = *i;
		std::string match_str = match.str();
		//std::cout << match_str << std::endl;
		vars.push_back(match_str);
	}
	std::string last_var = vars[vars.size() - 1];
	int n = std::atoi(last_var.substr(last_var.find("x_") + 2, last_var.size()).c_str()) + 1;

	for (size_t i = 0; i < n; i++) c.push_back(0.0);

	for (size_t i = 0; i < vars.size(); i++) {
		int index = std::atoi(vars[i].substr(vars[i].find("x_") + 2, vars[i].size()).c_str());
		c[index] = std::atof(vars[i].substr(0, vars[i].find("x_")).c_str());
		//std::cout << std::atof(vars[i].substr(0, vars[i].find("x_")).c_str()) << std::endl;
	}

	return c;
}

Expression parseLimitation(std::string input) {
	std::vector<double> c = parseExpression(input);
	Expression a;
	a.c = c;
	std::string ops[] = { "<=", ">=", "=", ">", "<" };
	for (std::string op : ops) {
		if (input.find(op) != std::string::npos) {
			a.op = op;
			break;
		}
	}
	a.b = std::atof(input.substr(input.find(a.op) + (a.op == ">=" || a.op == "<=" ? 2 : 1), input.size()).c_str());
	return a;
}

Expression standLimitation(Expression& exp) {
	if (exp.op == ">=" || exp.op == ">") {
		for (size_t i = 0; i < exp.c.size(); i++)
		{
			exp.c[i] = -exp.c[i];
		}
		exp.op = "<=";
		exp.b = -exp.b;
	}
	return exp;
}

SimplexTable generateSimplexTable(std::vector<double>& fc, std::vector<Expression>& limitations) {
	int n = limitations[0].c.size();
	SimplexTable table;
	for (int i = n - limitations.size(); i < n; ++i)
	{
		table.basis.push_back(std::string("x_") + std::to_string(i));
	}
	for (int i = 0; i < n; ++i)
	{
		table.vars.push_back(std::string("x_") + std::to_string(i));
		//std::cout << std::string("x_") + std::to_string(i) << std::endl;
	}
	table.vars.push_back("b");
	for (auto l : limitations) {
		table.rows.push_back(l.c);
		table.b.push_back(l.b);
	}

	table.delta = fc;
	for (size_t i = 0; i < table.delta.size(); i++)
	{
		table.delta[i] = -table.delta[i];
	}
	// целевая равна нулюб
	table.b.push_back(0);
	return table;
}

bool isOptimal(SimplexTable& table) {
	bool flag = true;
	for (auto el : table.b) {
		flag &= el > -eps;
	}
	return flag;
}


std::vector<double> saxpy(const std::vector<double>& a, const std::vector<double>& b, double alpha) {
	std::vector<double> c(a.size());
	for (size_t i = 0; i < a.size(); i++)
	{
		c[i] = (a[i] + alpha * b[i]);
	}

	return c;
}

std::vector<double> addVec(std::vector<double>& a, std::vector<double>& b) {
	return saxpy(a, b, 1);
}

std::vector<double> subVec(std::vector<double>& a, std::vector<double>& b) {
	return saxpy(a, b, -1);
}

std::vector<double> zeros(int n) {
	std::vector<double> c(n);
	return c;
}

std::vector<double> restoreOptimal(std::vector<int>& depVar, SimplexTable& table) {
	std::vector<double> c;
	for (auto i : depVar) {
		std::string varName = std::string("x_") + std::to_string(i);
		auto it = std::find(table.basis.begin(), table.basis.end(), varName);
		if (it != table.basis.end()) {
			int index = it - table.basis.begin();
			c.push_back(table.b[index]);
		}
		else {
			c.push_back(0);
		}
	}
	return c;
}

void doSimplexStep(SimplexTable& table) {
	double max = std::abs(table.b[0]);

	int index1 = 0;
	for (int i = 0; i < table.b.size(); ++i) {
		if (table.b[i] < 0 && std::abs(table.b[i]) > max) {
			max = std::abs(table.b[i]);
			index1 = i;
		}
	}
	double min = std::numeric_limits<double>::max();
	int index2 = -1;
	for (size_t i = 0; i < table.delta.size(); i++)
	{
		double t = table.delta[i] / table.rows[index1][i];
		if (!(std::isnan(t) || std::isinf(t) || std::abs(t) < eps) && std::abs(t) < min) {
			min = t;
			index2 = i;
		}
	}
	if (index2 == -1) {
		//
		exit(1);
	}


	for (size_t i = 0; i < table.rows.size(); i++)
	{
		if (i == index1) {
			table.b[i] = table.b[i] / table.rows[i][index2];
			table.rows[i] = saxpy(zeros(table.rows[i].size()), table.rows[i], 1 / table.rows[i][index2]);
		}
		else {
			table.b[i] = table.b[i] - table.rows[i][index2] / table.rows[index1][index2] * table.b[index1];
			table.rows[i] = saxpy(table.rows[i], table.rows[index1], -table.rows[i][index2] / table.rows[index1][index2]);
		}
		//std::cout << table << std::endl;
	}
	table.b[table.b.size() - 1] = table.b[table.b.size() - 1] - table.delta[index2] / table.rows[index1][index2] * table.b[index1];
	table.delta = saxpy(table.delta, table.rows[index1], -table.delta[index2] / table.rows[index1][index2]);
	table.basis[index1] = table.vars[index2];
}

std::vector<double> doSimplex(std::vector<double>& fc, std::vector<Expression>& limitations) {
	std::cout << "start simplex" << std::endl;
	int n = fc.size();
	std::vector<int> depVar;
	for (size_t i = 0; i < fc.size(); i++)
	{
		if (std::abs(fc[i]) > eps)
			depVar.push_back(i);
	}
	for (auto& l : limitations) {
		standLimitation(l);
		if (l.c.size() > n) n = l.c.size();
	}
	--n;
	for (size_t i = 0; i < limitations.size(); i++)
	{
		for (size_t j = 0; j < i; j++) limitations[i].c.push_back(0.0);
		limitations[i].c.push_back(1);
		for (size_t j = i; j < limitations.size() - 1; j++) limitations[i].c.push_back(0.0);
		limitations[i].op = "=";
	}
	for (size_t i = fc.size(); i < limitations[0].c.size(); i++)
	{
		fc.push_back(0.0);
	}
	for (auto l : limitations) {
	}

	SimplexTable table = generateSimplexTable(fc, limitations);
	std::cout << table << std::endl;
	while (!isOptimal(table)) {
		doSimplexStep(table);
		std::cout << table << std::endl;
	}
	auto optimal = restoreOptimal(depVar, table);
	for (size_t i = 0; i < depVar.size(); i++)
	{
		std::cout << "x_" + std::to_string(depVar[i]) << " = " << optimal[i] << std::endl;
	}
	return optimal;
}

template <typename T>
std::vector<T> mergeVectors(const std::vector<T>& v1, const std::vector<T>& v2) {
	std::vector<T> res;
	res.reserve(v1.size() + v2.size()); // preallocate memory
	res.insert(res.end(), v1.begin(), v1.end());
	res.insert(res.end(), v2.begin(), v2.end());
	return res;
}

std::ostream& operator<< (std::ostream& os, const Expression& exp)
{
	for (size_t i = 0; i < exp.c.size(); i++)
	{
		os << exp.c[i] << " ";
	}
	os << exp.op << " " << exp.b;
	return os;
}

std::ostream& operator<< (std::ostream& os, const SimplexTable& table)
{
	samilton::ConsoleTable tableB;
	std::vector<std::string> row;

	row.push_back("Basis");
	for (auto h : table.vars) {
		row.push_back(h);
	}
	tableB.addRow(std::vector<std::string>(row));
	row.clear();
	for (size_t i = 0; i < table.basis.size(); i++)
	{
		row.push_back(table.basis[i]);
		for (auto h : table.rows[i]) {
			row.push_back(std::to_string(h));
		}
		row.push_back(std::to_string(table.b[i]));
		tableB.addRow(std::vector<std::string>(row));
		row.clear();
	}
	row.push_back("delta");
	for (auto h : table.delta) {
		row.push_back(std::to_string(h));
	}
	row.push_back(std::to_string(table.b[table.b.size() - 1]));
	tableB.addRow(std::vector<std::string>(row));
	os << tableB;
	return os;
}