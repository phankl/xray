#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <numbers>

#include <complex>
#include <fftw3.h>

#include "analysis.h"
#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure{
  
  public:
   
    Structure(int, vector<int>, vector<double>, XYZ, XYZ, XYZ, XYZ, int, bool);
    
    vector<vector<vector<double>>> density;
    vector<vector<vector<double>>> intensity;

    vector<vector<double>> density2D;
    vector<vector<double>> intensity2D;
    
    void savePoints(string);
    void saveDensity(string);
    void saveIntensity(string);
    void ewaldSphere(XYZ, string);

    vector<double> odf(int);
    vector<double> intensityDistributionTube(double, int);
    vector<double> intensityDistribution2D(double, int);

  private:

    int nTube;
    double rTube;
    double lTube;

    bool saxs;

    XYZ boxSize;
    XYZ voxelSize;
    XYZ voxelFFTSize;
    XYZ nVoxel;

    XYZ mean;
    XYZ std;

    int zPoints;
    int phiPoints;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;
    vector<XYZ> points;
     
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate(int);
    void voxelise();
    void fft();
};

#endif
