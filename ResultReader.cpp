#include "ResultReader.h"
#include <iomanip>

void Report::create_report(string instance, vector<string> nGroup, vector<string> capacity, vector<string> noPay) {
	report_file = instance + ".csv";
	ofstream reportFile(report_file, ofstream::out);
	reportFile << "noPay,nGroup,capacity,obj,runtime\n";
	for (int n = 0; n < noPay.size(); ++n) {
		for (int b = 0; b < nGroup.size(); ++b) {
			for (int c = 0; c < capacity.size(); ++c) {
				ifstream resultFile;
				string name = instance + "_" + nGroup[b] + "_" + capacity[c] + "_" + noPay[n] + ".txt";
				cout << name << endl;
				resultFile.open(name, ofstream::in);
				double tmp, obj, time;
				resultFile >> tmp >> obj >> time;
				if (!resultFile) {
					cout << "Cannot open file " << name << " ..." << endl;
				}
				resultFile.close();
				reportFile << setprecision(5) << fixed << noPay[n] << "," << nGroup[b] << "," << capacity[c] << "," << obj << "," << time << endl;
			}
		}
	}
	reportFile.close();
}

void Report::create_report_Sen(string instance, vector<string> noPay, vector<string> capacity) {
	report_file = instance + ".csv";
	vector<string> set = { "set1", "set2", "set3", "set4", "set5"};
	ofstream reportFile(report_file, ofstream::out);
	reportFile << "noPay,capacity,set,obj,runtime\n";
	for (int n = 0; n < noPay.size(); ++n)
		for (int c = 0; c < capacity.size(); ++c)
			for (int s = 0; s < set.size(); ++s) {
				ifstream resultFile;
				string name = instance + "_" + set[s] + "_" + noPay[n] + "_" + capacity[c] + ".txt";
				cout << name << endl;
				resultFile.open(name, ofstream::in);
				double tmp, obj, time;
				resultFile >> tmp >> obj >> time;
				if (!resultFile) {
					cout << "Cannot open file " << name << " ..." << endl;
				}
				resultFile.close();
				reportFile << setprecision(5) << fixed << noPay[n] << "," << capacity[c] << "," << set[s] << "," << obj << "," << time << endl;
			}

	reportFile.close();
}