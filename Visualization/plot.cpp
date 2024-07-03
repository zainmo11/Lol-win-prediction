#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <cmath>
#include <map>
#include <boost/algorithm/string.hpp>
#include "plot.h"
using namespace std;

mutex mtx;

void generatePlot(const string& dataFilePath, const string& outputFolder, int i, int j) {
    string scriptFilename = outputFolder + "/pairplot_" + to_string(i+1) + "_vs_" + to_string(j+1) + ".gp";
    string outputFilename = outputFolder + "/pairplot_" + to_string(i+1) + "_vs_" + to_string(j+1) + ".png";

    ofstream scriptFile(scriptFilename);
    scriptFile << "set terminal pngcairo size 800,600\n";
    scriptFile << "set output '" << outputFilename << "'\n";
    scriptFile << "set title 'Variable " << i + 1 << " vs Variable " << j + 1 << "'\n";
    scriptFile << "plot '" << dataFilePath << "' using " << j + 1 << ":" << i + 1 << " with points title ''\n";
    scriptFile.close();

    // Run Gnuplot script
    string command = "gnuplot " + scriptFilename;
    system(command.c_str());

    // Lock mutex before printing to avoid interleaved output
    mtx.lock();
    cout << "Generated plot: " << outputFilename << endl;
    cout << "Generated script: " << scriptFilename << endl;
    mtx.unlock();
}

void generatePairplot(const vector<vector<double>>& data, const string& outputFolder) {
    int numVariables = data[0].size();

    // Create the output folder if it doesn't exist
    string mkdirCommand = "mkdir -p " + outputFolder;
    system(mkdirCommand.c_str());

    // Print the output folder path
    cout << "Output folder: " << outputFolder << endl;

    // Write data to file
    string dataFilePath = "pairplot_data.dat";
    ofstream dataFile(dataFilePath);
    for (const auto& row : data) {
        for (const auto& val : row) {
            dataFile << val << " ";
        }
        dataFile << "\n";
    }
    dataFile.close();

    // Print the data file path
    cout << "Data file path: " << dataFilePath << endl;

    // Vector to hold threads
    vector<thread> threads;

    // Generate Gnuplot script for each pair of variables
    for (int i = 0; i < numVariables; ++i) {
        for (int j = i + 1; j < numVariables; ++j) {
            // Create a thread for each plot generation task
            threads.emplace_back(generatePlot, dataFilePath, outputFolder, i, j);
        }
    }

    // Wait for all threads to complete
    for (auto& th : threads) {
        th.join();
    }

    cout << "All pairplots generated in folder: " << outputFolder << endl;
}


void generateGnuplotScript(std::vector<std::vector<double>>& data) {
    std::ofstream scriptFile("heatmap_script.gp");
    if (!scriptFile.is_open()) {
        std::cerr << "Error: Unable to create Gnuplot script file." << std::endl;
        return;
    }

    // Write Gnuplot commands to the script file
    scriptFile << "set terminal png\n"; // Output to PNG file
    scriptFile << "set output 'heatmap.png'\n";
    scriptFile << "set xrange [0:" << data[0].size() - 1 << "]\n";
    scriptFile << "set yrange [0:" << data.size() - 1 << "]\n";
    scriptFile << "set pm3d map\n"; // Heatmap style plot
    scriptFile << "unset key\n";
    scriptFile << "set title 'Heat Map'\n";
    scriptFile << "plot '-' matrix with image\n";

    // Write data to plot into the script file
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            scriptFile << data[i][j] << " ";
        }
        scriptFile << "\n";
    }
    scriptFile << "e\n";

    scriptFile.close();
}

// Function to calculate the density for each data point
void calculateDensity(const vector<vector<double>>& data, double R, vector<vector<double>>& densityData) {
    int numDataPoints = data.size();
    for (int i = 0; i < numDataPoints; ++i) {
        double x0 = data[i][0];
        double y0 = data[i][1];
        int c = 0;
        for (int j = 0; j < numDataPoints; ++j) {
            if (i != j) {
                double x1 = data[j][0];
                double y1 = data[j][1];
                double distance = sqrt((x1 - x0) * (x1 - y0) + (y1 - y0) * (y1 - y0));
                if (distance <= R) {
                    ++c;
                }
            }
        }
        double density = c / (M_PI * R * R);
        densityData.push_back({x0, y0, density});
    }
}

// Function to generate the Gnuplot script and run it
void generateDensityPlot(const string& outputFolder, const vector<vector<double>>& densityData, int i, int j, int max_x , int max_y )  {
    string scriptFilename = outputFolder + "/densityplot_" + to_string(i+1) + "_vs_" + to_string(j+1) + ".gp";
    string outputFilename = outputFolder + "/densityplot_" + to_string(i+1) + "_vs_" + to_string(j+1) + ".png";

    ofstream scriptFile(scriptFilename);
    scriptFile << "reset session\n";
    scriptFile << "set size ratio -1\n";
    scriptFile << "set palette rgb 33,13,10\n";
    scriptFile << "set terminal pngcairo size 800,600\n";
    scriptFile << "set output '" << outputFilename << "'\n";
    scriptFile << "# Set the axis ranges to zoom out\n";
    scriptFile << "set xrange [0:"<< max_x <<"]\n";
    scriptFile << "set yrange [0:"<< max_y <<"]\n";
    scriptFile << "plot '-' using 1:2:3 with points pt 7 lc palette z notitle\n";

    for (const auto& row : densityData) {
        scriptFile << row[0] << " " << row[1] << " " << row[2] << "\n";
    }

    scriptFile << "e\n";
    scriptFile.close();

    // Run Gnuplot script
    string command = "gnuplot " + scriptFilename;
    system(command.c_str());

    cout << "Generated density plot: " << outputFilename << endl;
    cout << "Generated script: " << scriptFilename << endl;
}

// Function to generate pair plots with density
void generatePairplotWithDensity(vector<vector<double>>& data, const vector<string>& hueData, const string& outputFolder, double R, int max_x, int max_y) {
    int numVariables = data[0].size();

    // Create the output folder if it doesn't exist
    string mkdirCommand = "mkdir -p " + outputFolder;
    system(mkdirCommand.c_str());

    // Print the output folder path
    cout << "Output folder: " << outputFolder << endl;

    // Categorize data based on hue
    map<string, vector<vector<double>>> categorizedData;
    for (size_t i = 0; i < data.size(); ++i) {
        categorizedData[hueData[i]].push_back(data[i]);
    }

    // Print the categorized data file paths
    for (const auto& category : categorizedData) {
        cout << "Data file path for category '" << category.first << "': " << outputFolder + "/data_" + category.first + ".dat" << endl;
    }

    // Generate Gnuplot script for each pair of variables
    for (int i = 0; i < numVariables; ++i) {
        for (int j = i + 1; j < numVariables; ++j) {
            // Calculate density for each category
            for (const auto& category : categorizedData) {
                vector<vector<double>> densityData;
                calculateDensity(category.second, R, densityData);
                generateDensityPlot(outputFolder, densityData, i, j , max_x, max_y);
            }
        }
    }

    cout << "All density plots generated in folder: " << outputFolder << endl;
}

void LinearRegPlot(const std::string& outputFolder, int i, int j, double a0, double b, const std::vector<std::pair<double, double>>& data) {
    std::string scriptFilename = outputFolder + "/pairplot_" + std::to_string(i + 1) + "_vs_" + std::to_string(j + 1) + ".gp";
    std::string outputFilename = outputFolder + "/pairplot_" + std::to_string(i + 1) + "_vs_" + std::to_string(j + 1) + ".png";

    std::ofstream scriptFile(scriptFilename);
    scriptFile << "reset session\n";

    scriptFile << "$Data <<EOD\n";
    for (const auto& point : data) {
        scriptFile << point.first << " " << point.second << "\n";
    }
    scriptFile << "EOD\n";

    scriptFile << "b = " << b << "\n";
    scriptFile << "f(x) = " << a0 << "*x + b\n";

    // Least square fit
    scriptFile << "set fit nolog quiet\n";
    scriptFile << "fit f(x) $Data u 1:2 via a0\n";

    scriptFile << "set terminal pngcairo size 800,600\n";
    scriptFile << "set output '" << outputFilename << "'\n";
    scriptFile << "set key top left\n";
    scriptFile << "set title 'Variable " << i + 1 << " vs Variable " << j + 1 << "'\n";

    scriptFile << "plot $Data u 1:2 w p pt 7, \\\n";
    scriptFile << "     f(x) w l ti sprintf(\"Least square fit: f(x)="<<a0<<"*x+ "<<b<<", a0=%%g\",a0)\n";
    scriptFile << "     $Data using 1:2:(0):(f($1)-$2) with vectors linecolor rgb \"red\" nohead title \"Differences to f(x)\"\n";
    scriptFile.close();

    // Run Gnuplot script
    std::string command = "gnuplot " + scriptFilename;
    system(command.c_str());

    // Lock mutex before printing to avoid interleaved output

    std::cout << "Generated plot: " << outputFilename << std::endl;
    std::cout << "Generated script: " << scriptFilename << std::endl;

}

void BoxPlot(const std::string& outputFolder, const std::vector<std::vector<std::string>>& data,
             const std::string& title, const std::string& xLabel, const std::string& yLabel) {
    if (data.empty() || data[0].empty()) {
        std::cerr << "Empty data provided." << std::endl;
        return;
    }

    // Extract feature names from the first row
    std::vector<std::string> featureNames = data[0];
    std::string scriptFilename = outputFolder + "/boxplot.gp";
    std::string dataFilename = outputFolder + "/boxplot_data.dat";
    std::string outputFilename = outputFolder + "/boxplot.png";

    // Write numerical data to file, starting from the second row
    std::ofstream dataFile(dataFilename);
    for (size_t i = 1; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            dataFile << data[i][j];
            if (j < data[i].size() - 1) {
                dataFile << " ";  // Space delimiter
            }
        }
        dataFile << "\n";
    }
    dataFile.close();

    // Write Gnuplot script
    std::ofstream scriptFile(scriptFilename);
    scriptFile << "reset session\n";
    scriptFile << "set border 2 front lt black linewidth 1.000 dashtype solid\n";
    scriptFile << "set boxwidth 0.5 absolute\n";
    scriptFile << "set style fill solid 0.50 border lt -1\n";
    scriptFile << "unset key\n";
    scriptFile << "set pointsize 0.5\n";
    scriptFile << "set style data boxplot\n";
    scriptFile << "set xtics border in scale 0,0 nomirror norotate autojustify\n";
    scriptFile << "set xtics norangelimit\n";

    // Construct the xtics string from feature names
    scriptFile << "set xtics (";
    for (size_t i = 0; i < featureNames.size(); ++i) {
        scriptFile << "\"" << featureNames[i] << "\" " << (i + 1);
        if (i != featureNames.size() - 1) {
            scriptFile << ", ";
        }
    }
    scriptFile << ")\n";

    scriptFile << "set ytics border in scale 1,0.5 nomirror norotate autojustify\n";
    scriptFile << "set title \"" << title << "\\n\"\n";
    scriptFile << "set xlabel \"" << xLabel << "\"\n";
    scriptFile << "set ylabel \"" << yLabel << "\"\n";
    scriptFile << "set xrange [ * : * ] noreverse writeback\n";
    scriptFile << "set x2range [ * : * ] noreverse writeback\n";
    scriptFile << "set yrange [ * : * ] noreverse writeback\n";
    scriptFile << "set y2range [ * : * ] noreverse writeback\n";
    scriptFile << "set zrange [ * : * ] noreverse writeback\n";
    scriptFile << "set cbrange [ * : * ] noreverse writeback\n";
    scriptFile << "set rrange [ * : * ] noreverse writeback\n";
    scriptFile << "set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front noinvert bdefault\n";
    scriptFile << "set style boxplot candles range 1.50 nooutliers pt 7 separation 1 labels auto sorted\n";
    scriptFile << "t(x) = x/1.e6\n";
    scriptFile << "filter(col, factor_col, level) = (strcol(factor_col) eq word(factors, level)) ? t(column(col)) : 1/0\n";
    scriptFile << "NO_ANIMATION = 1\n";
    scriptFile << "NF = " << featureNames.size() << "\n";
    scriptFile << "Last datafile plotted: \"" << dataFilename << "\"\n";
    scriptFile << "plot for [i=1:NF] '" << dataFilename << "' using (i):(filter(8, 4, i))\n";

    scriptFile << "set terminal pngcairo size 800,600\n";
    scriptFile << "set output '" << outputFilename << "'\n";
    scriptFile << "replot\n";
    scriptFile.close();

    // Run Gnuplot script
    std::string command = "gnuplot " + scriptFilename;
    system(command.c_str());

    // Lock mutex before printing to avoid interleaved output
    std::cout << "Generated plot: " << outputFilename << std::endl;
    std::cout << "Generated script: " << scriptFilename << std::endl;
}
