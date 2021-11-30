#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <complex>
#include <fftw3.h>

#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure{
  
  public:
   
    Structure(int, vector<int>, vector<double>, XYZ, XYZ, XYZ, XYZ);
    void savePoints(string);
    void saveDensity(string);
    void saveIntensity(string);

  private:

    int nTube;
    double rTube;
    double lTube;

    XYZ boxSize;
    XYZ voxelSize;
    XYZ nVoxel;

    XYZ mean;
    XYZ std;

    int zPoints;
    int phiPoints;

    vector<Tube> tubes;
    vector<XYZ> points;
    vector<vector<vector<double>>> density;
    vector<vector<vector<double>>> intensity;

    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate();
    void voxelise();
    void fft();
};

#endif
