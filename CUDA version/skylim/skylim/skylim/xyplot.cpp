// xyplot.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>  // For std::pair
#include <cstdlib>  // For system()
#include <string>

using namespace std;

// Function to plot data using Gnuplot

void plotData(const string& filename, char Mtitle[],char Xtitle[],char Ytitle[],char lg1[],char lg2[],char lg3[],char lg4[]) {
    FILE* gnuplotPipe = _popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        cerr << "Error: Could not open Gnuplot." << endl;
        return;
    }

    // Gnuplot settings
    fprintf(gnuplotPipe, "set title '%s'\n",Mtitle);
    fprintf(gnuplotPipe, "set xlabel '%s'\n", Xtitle);
    fprintf(gnuplotPipe, "set ylabel '%s'\n", Ytitle);
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "plot '%s' using 1:2 with lines linecolor 'green' linewidth 2 title '%s', \\\n", filename.c_str(), lg1);
    fprintf(gnuplotPipe, "     '%s' using 1:3 with lines linecolor 'red' linewidth 2 title '%s', \\\n", filename.c_str(), lg2);
    fprintf(gnuplotPipe, "     '%s' using 1:4 with lines linecolor 'blue' linewidth 2 title '%s', \\\n", filename.c_str(), lg3);
    fprintf(gnuplotPipe, "     '%s' using 1:5 with lines linecolor 'black' linewidth 2 title '%s' \n", filename.c_str(), lg4);
    fflush(gnuplotPipe);
    _pclose(gnuplotPipe);
}

void thexyplot(int mpl, double* timpl, double* vpl1, double* vpl2, double* vpl3, double* vpl4, char Mtitle[], char Xtitle[], char Ytitle[], char lg1[], char lg2[], char lg3[], char lg4[])
//void thexyplot() 
{
    int k;
    FILE* mydata;
    mydata = fopen("datatoplot.txt", "w");
    for (k = 1; k <= mpl; k++)
    {
        fprintf(mydata, "%.5e %.5e %.5e %.5e %.5e \n", timpl[k], vpl1[k], vpl2[k], vpl3[k], vpl4[k]);
    }
    fclose(mydata);

    string filename = "datatoplot.txt";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    vector<vector<double>> data;
    double x, y1, y2, y3, y4;
    // Read data from file
    while (file >> x >> y1 >> y2 >> y3 >> y4) {
        data.push_back({ x, y1, y2, y3, y4 });
    }
    file.close();

    // Check if data is available
    if (data.empty()) {
        cerr << "Error: No valid data found in the file!" << endl;
        return;
    }

//    cout << "Data loaded successfully from " << filename << endl;
    // Plot the data using Gnuplot
    plotData(filename, Mtitle, Xtitle, Ytitle, lg1, lg2, lg3, lg4);
      return;
}
