//
// Created by zyn66 on 6/15/2024.
//

#ifndef PLOT_H
#define PLOT_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

void generatePairplot(const vector<vector<double>>& data, const string& outputFolder);
void generateGnuplotScript(std::vector<std::vector<double>>& data);
void generateDensityPlot(const string& outputFolder, const vector<vector<double>>& densityData, int i, int j);
void generatePairplotWithDensity(vector<vector<double>>& data, const vector<string>& hueData, const string& outputFolder, double R, int max_x, int max_y);
void LinearRegPlot(const std::string& outputFolder, int i, int j, double a0, double b, const std::vector<std::pair<double, double>>& data);
void BoxPlot(const std::string& outputFolder, const std::vector<std::vector<std::string>>& data,
             const std::string& title, const std::string& xLabel, const std::string& yLabel);
#endif //PLOT_H
