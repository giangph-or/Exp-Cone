#include <iostream>
#include <vector>
#include <fstream>
#include "ProblemReader.h"

using namespace std;
class Report
{
public:
	int nProducts;
	string sol_file;
	string report_file;
	void create_report(string instance, vector<string> nGroup, vector<string> capacity, vector<string> noPay);
	void create_report_Sen(string instance, vector<string> noPay, vector<string> capacity);
};
#pragma once
