#include "structure.h"

Structure::Structure(int nTubeNew, vector<int> rasterProperties, vector<double> tubeProperties, XYZ boxSizeNew, XYZ nVoxelNew, XYZ meanNew, XYZ stdNew, int seed, bool saxsNew):
  nTube(nTubeNew),
  rTube(tubeProperties[0]),
  lTube(tubeProperties[1]),
  phiPoints(rasterProperties[0]),
  zPoints(rasterProperties[1]),
  boxSize(boxSizeNew),
  nVoxel(nVoxelNew),
  mean(meanNew),
  std(stdNew),
  saxs(saxsNew)
{
  generate(seed);
  voxelise();
  fft();
}

void Structure::savePoints(string fileName) {
  ofstream file(fileName);

  file << points.size() << endl;
  file << "Atoms. Timestep: 0" << endl;

  for (int i = 0; i < points.size(); i++) {
    XYZ point = points[i];
    file << 1 << " " << point.x << " " << point.y << " " << point.z << endl;
  }
  file.close();
}

void Structure::saveDensity(string fileName) {
  ofstream file(fileName);

  if (saxs) {
    
    // VTK file header

    file << "# vtk DataFile Version 4.2" << endl;
    file << "Charge density voxel data" << endl;
    file << "ASCII" << endl;

    // define data structure

    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << nVoxel.x << " 1 " << nVoxel.z << endl;
    file << "ORIGIN 0.0 0.0 0.0" << endl;
    file << "SPACING " << voxelSize.x << " 0.0 " << voxelSize.z << endl;

    // data block

    file << "POINT_DATA " << int(nVoxel.x * nVoxel.z) << endl;
    file << "SCALARS DENSITY FLOAT 1" << endl;
    file << "LOOKUP_TABLE DENSITY" << endl;

    for (int j = 0; j < nVoxel.z; j++)
      for (int i = 0; i < nVoxel.x; i++)
        file << density2D[i][j] << endl;
  }
  else {

    // VTK file header

    file << "# vtk DataFile Version 4.2" << endl;

    file << "ASCII" << endl;

    // define data structure

    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << nVoxel.x << " " << nVoxel.y << " " << nVoxel.z << endl;
    file << "ORIGIN 0.0 0.0 0.0" << endl;
    file << "SPACING " << voxelSize.x << " " << voxelSize.y << " " << voxelSize.z << endl;

    // data block

    file << "POINT_DATA " << int(nVoxel.x * nVoxel.y * nVoxel.z) << endl;
    file << "SCALARS DENSITY FLOAT 1" << endl;
    file << "LOOKUP_TABLE DENSITY" << endl;

    for (int k = 0; k < nVoxel.z; k++)
      for (int j = 0; j < nVoxel.y; j++)
        for (int i = 0; i < nVoxel.z; i++)
          file << density[i][j][k] << endl;
  }

  file.close();
}

void Structure::saveIntensity(string fileName) {
  ofstream file(fileName);

  if (saxs) {
    
    // VTK file header

    file << "# vtk DataFile Version 4.2" << endl;
    file << "X-ray intensity voxel data" << endl;
    file << "ASCII" << endl;

    // define data structure

    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << nVoxel.x << " 1 " << nVoxel.z << endl;
    file << "ORIGIN " << -0.5*voxelFFTSize.x*nVoxel.x << " 0.0 " << -0.5*voxelFFTSize.z*nVoxel.z <<  endl;
    file << "SPACING " << voxelFFTSize.x << " 0.0 " << voxelFFTSize.z << endl;

    // data block

    file << "POINT_DATA " << int(nVoxel.x * nVoxel.z) << endl;
    file << "SCALARS INTENSITY FLOAT 1" << endl;
    file << "LOOKUP_TABLE INTENSITY" << endl;

    for (int j = 0; j < nVoxel.z; j++)
      for (int i = 0; i < nVoxel.z; i++)
        file << intensity2D[i][j] << endl;

  }
  else {

    // VTK file header

    file << "# vtk DataFile Version 4.2" << endl;
    file << "X-ray intensity voxel data" << endl;
    file << "ASCII" << endl;

    // define data structure

    file << "DATASET STRUCTURED_POINTS" << endl;
    file << "DIMENSIONS " << nVoxel.x << " " << nVoxel.y << " " << nVoxel.z << endl;
    file << "ORIGIN " << -0.5*voxelFFTSize.x*nVoxel.x << " " << -0.5*voxelFFTSize.y*nVoxel.y << " " << -0.5*voxelFFTSize.z*nVoxel.z <<  endl;
    file << "SPACING " << voxelFFTSize.x << " " << voxelFFTSize.y << " " << voxelFFTSize.z << endl;

    // data block

    file << "POINT_DATA " << int(nVoxel.x * nVoxel.y * nVoxel.z) << endl;
    file << "SCALARS INTENSITY FLOAT 1" << endl;
    file << "LOOKUP_TABLE INTENSITY" << endl;

    for (int k = 0; k < nVoxel.z; k++)
      for (int j = 0; j < nVoxel.y; j++)
        for (int i = 0; i < nVoxel.z; i++)
          file << intensity[i][j][k] << endl;
  }

  file.close();
}

void Structure::ewaldSphere(XYZ kIn, string fileName) {
  
  if (saxs) {
    cout << "Ewald sphere not used for SAXS simulation!" << endl;
    return;
  }

  double r = kIn.length();

  vector<XYZ> voxels;
  vector<double> voxelData;

  // find voxelised sphere

  for (int i = 0; i < nVoxel.x; i++)
    for (int j = 0; j < nVoxel.y; j++)
      for (int k = 0; k < nVoxel.z; k++) {
        
        double x = (i - 0.5*nVoxel.x)*voxelFFTSize.x;
        double y = (j - 0.5*nVoxel.y)*voxelFFTSize.y;
        double z = (k - 0.5*nVoxel.z)*voxelFFTSize.z;
        XYZ pos(x, y, z);
        
        double d = (pos - kIn).length();
        double dMin = d;
        double dMax = d;

        // check all vertices of each voxel
        
        double dTemp;

        XYZ posZ = pos + XYZ(0.0, 0.0, voxelFFTSize.z);
        dTemp = (posZ - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
        
        XYZ posY = pos + XYZ(0.0, voxelFFTSize.y, 0.0);
        dTemp = (posY - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
 
        XYZ posYZ = pos + XYZ(0.0, voxelFFTSize.y, voxelFFTSize.z);
        dTemp = (posYZ - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
        
        XYZ posX = pos + XYZ(voxelFFTSize.x, 0.0, 0.0);
        dTemp = (posX - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
        
        XYZ posXZ = pos + XYZ(voxelFFTSize.x, 0.0, voxelFFTSize.z);
        dTemp = (posXZ - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
        
        XYZ posXY = pos + XYZ(voxelFFTSize.x, voxelFFTSize.y, 0.0);
        dTemp = (posXY - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;
        
        XYZ posXYZ = pos + XYZ(voxelFFTSize.x, voxelFFTSize.y, voxelFFTSize.z);
        dTemp = (posXYZ - kIn).length();
        if (dTemp < dMin) dMin = dTemp;
        else if (dTemp > dMax) dMax = dTemp;

        // if voxel intersects Ewald sphere, add voxel to data in correct order

        if (dMin <= r && r <= dMax) {
          voxels.push_back(pos);
          voxels.push_back(posX);
          voxels.push_back(posY);
          voxels.push_back(posXY);
          voxels.push_back(posZ);
          voxels.push_back(posXZ);
          voxels.push_back(posYZ);
          voxels.push_back(posXYZ);

          voxelData.push_back(intensity[i][j][k]);
        }
      }

  // output Ewald sphere to file
  
  ofstream file(fileName);

  // VTK file header

  file << "# vtk DataFile Version 4.2" << endl;
  file << "Ewald sphere voxel data" << endl;
  file << "ASCII" << endl;

  // define data structure

  file << "DATASET UNSTRUCTURED_GRID" << endl;

  // cell points

  file << "POINTS " << voxels.size() << " FLOAT" << endl;

  for (int i = 0; i < voxels.size(); i++)
    file << voxels[i].x << " " << voxels[i].y << " " << voxels[i].z << endl;

  // cell definitions

  file << "CELLS " << voxels.size()/8 << " " << 9*voxels.size()/8 << endl;

  for (int i = 0; i < voxels.size() / 8; i++) {
    file << "8 ";
    for (int j = 0; j < 8; j++)
      file << 8*i + j << " ";
    file << endl;
  }

  file << "CELL_TYPES " << voxels.size() / 8 << endl;

  for (int i = 0; i < voxels.size() / 8; i++)
    file << 11 << endl;

  // data block

  file << "CELL_DATA " << voxelData.size() << endl;
  file << "SCALARS INTENSITY FLOAT 1" << endl;
  file << "LOOKUP_TABLE INTENSITY" << endl;

  for (int i = 0; i < voxelData.size(); i++)
    file << voxelData[i] << endl;

  file.close();

}

vector<double> Structure::odf(int nBins) {
  
  vector<double> result(nBins, 0.0);
  double binSize = numbers::pi / nBins;

  for (int i = 0; i < nTube; i++) {
    XYZ t = tubes[i].t;
    double angle = acos(t.z / t.length());
    int bin = floor(angle / binSize);
    result[bin] += 1.0;
  }

  double sum = 0.0;
  for (int i = 0; i < nBins; i++)
    sum += result[i] * binSize;

  for (int i = 0; i < nBins; i++)
    result[i] /= sin((i+0.5) * binSize) * sum;

  return result;
}

vector<double> Structure::intensityDistribution2D(double radius, int nBins) {
  
  vector<double> result(nBins, 0.0);
  vector<int> binCount(nBins, 0);
  double binSize = 2.0 * numbers::pi / nBins;

  for (int i = 0; i < nVoxel.x; i++)
    for (int j = 0; j < nVoxel.z; j++) {
      double x = (i - 0.5*nVoxel.x)*voxelFFTSize.x;
      double z = (j - 0.5*nVoxel.z)*voxelFFTSize.z;
    
      double d = x*x + z*z;
      double dMin = d;
      double dMax = d;

      // check all vertices of pixel

      d = x*x + pow(z+voxelFFTSize.z, 2);
      if (d < dMin) dMin = d;
      if (d > dMax) dMax = d;
      
      d = pow(x+voxelFFTSize.x, 2) + z*z;
      if (d < dMin) dMin = d;
      if (d > dMax) dMax = d;
      
      d = pow(x+voxelFFTSize.x, 2) + pow(z+voxelFFTSize.z, 2);
      if (d < dMin) dMin = d;
      if (d > dMax) dMax = d;

      // add intensity to bin if pixel contains desired contour

      if (dMin < radius*radius && radius*radius < dMax) {

        double xCentre = x + 0.5*voxelFFTSize.x;
        double zCentre = z + 0.5*voxelFFTSize.z;

        double angle = -atan2(zCentre, xCentre) + 1.5*numbers::pi;
        if (angle > 2.0 * numbers::pi) angle -= 2.0 * numbers::pi;
        int bin = floor(angle / binSize);
        result[bin] += intensity2D[i][j];
        binCount[bin]++;
      }
    }

  for (int i = 0; i < nBins; i++)
    result[i] /= binCount[i];

  return result;
}


double Structure::distance(Tube tube1, Tube tube2, double lambda1, double lambda2) {
  return (tube1.s - tube2.s + lambda1*tube1.t - lambda2*tube2.t).length();
}

bool Structure::checkOverlap(Tube tube1) {
  
  double r1 = tube1.r;
  double l1 = tube1.l;

  XYZ t1 = tube1.t;
  XYZ s1 = tube1.s;

  for (int i = 0; i < ghostTubes.size(); i++) {   
    
    Tube tube2 = ghostTubes[i];
    
    double r2 = tube2.r;
    double l2 = tube2.l;

    XYZ t2 = tube2.t;
    XYZ s2 = tube2.s;

    // check for minimum distance of lines

    XYZ n = cross(t1, t2);
    n.normalise();

    XYZ delta = s2 - s1;

    double d = fabs(delta * n);

    if (d > r1+r2) continue;

    // check if minimum distance is within line segments

    double c = t1 * t2;
    double frac = 1.0 / (1.0 - c*c);
    
    double lambda1 = (t1 - c*t2) * delta * frac;
    double lambda2 = (c*t1 - t2) * delta * frac;
 
    if (lambda1 >= 0.0 && lambda1 <= l1 && lambda2 >= 0.0 && lambda2 <= l2) {
      if (d < r1+r2) return true;
      else continue;
    }

    // check cases where minimum distance is outside the line segments
    // not fully exploiting convexity on subdomains yet, kinda brute force

    // domain edges
    
    lambda1 = 0.0;
    lambda2 = -(delta * t2);
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = c*l1 - delta*t2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = delta * t1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    lambda1 = c*l2 + delta*t1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    // domain vertices

    lambda1 = 0.0;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = 0.0;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
    
    lambda1 = l1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
  }

  return false;
}

void Structure::generate(int seed) {

  mt19937_64 generator(seed);
  generator.discard(10000);

  uniform_real_distribution<double> sDistribution(0.0, 1.0);
  normal_distribution<double> tDistributionX(mean.x, std.x);
  normal_distribution<double> tDistributionY(mean.y, std.y);
  normal_distribution<double> tDistributionZ(mean.z, std.z);

  int attempts = 0;
  int tubeIndex = 1;

  while (tubes.size() < nTube) {
    
    // generate new tube

    double sX = sDistribution(generator) * boxSize.x;
    double sY = sDistribution(generator) * boxSize.y;
    double sZ = sDistribution(generator) * boxSize.z;
  
    double tX = tDistributionX(generator);
    double tY = tDistributionY(generator);
    double tZ = tDistributionZ(generator);
 
    // if (sDistribution(generator) < 0.5) tZ *= -1.0;

    XYZ s(sX, sY, sZ);
    XYZ t(tX, tY, tZ);

    t.normalise();

    Tube tube(rTube, lTube, s, t);
    
    // check periodic boundary conditions and create secondary ghost tube for trial

    XYZ eTemp = s + lTube*t;
    XYZ sTemp = s;
    
    bool ghost = false;

    if (eTemp.x < 0.0) {
      sTemp.x += boxSize.x;
      ghost = true;
    }
    if (eTemp.y < 0.0) {
      sTemp.y += boxSize.y; 
      ghost = true;
    }
    if (eTemp.z < 0.0) {
      sTemp.z += boxSize.z; 
      ghost = true;
    }
    if (eTemp.x > boxSize.x) {
      sTemp.x -= boxSize.x;
      ghost = true;
    }
    if (eTemp.y > boxSize.y) {
      sTemp.y -= boxSize.y;
      ghost = true;
    }
    if (eTemp.z > boxSize.z) {
      sTemp.z -= boxSize.z;
      ghost = true;
    }

    attempts++;

    // add new tube if no overlap is detected

    bool ghostOverlap = (ghost) ? checkOverlap(Tube(rTube, lTube, sTemp, t)) : false;

    if (!checkOverlap(tube) && !ghostOverlap) {
      
      // cout << "Tube " << tubeIndex << " generated after " << attempts << " attempts." << endl;
      attempts = 0;
      tubes.push_back(tube);
      ghostTubes.push_back(tube);
      tubeIndex++;

      // check periodic boundary conditions and add ghost tube if boundary is crossed
      
      XYZ eTemp = s + lTube*t;
      XYZ sTemp = s;
      
      bool ghost = false;

      if (eTemp.x < 0.0) {
        sTemp.x += boxSize.x;
        ghost = true;
      }
      if (eTemp.y < 0.0) {
        sTemp.y += boxSize.y; 
        ghost = true;
      }
      if (eTemp.z < 0.0) {
        sTemp.z += boxSize.z; 
        ghost = true;
      }
      if (eTemp.x > boxSize.x) {
        sTemp.x -= boxSize.x;
        ghost = true;
      }
      if (eTemp.y > boxSize.y) {
        sTemp.y -= boxSize.y;
        ghost = true;
      }
      if (eTemp.z > boxSize.z) {
        sTemp.z -= boxSize.z;
        ghost = true;
      }

      if (ghost) ghostTubes.push_back(Tube(rTube, lTube, sTemp, t));
    }
  }
}

void Structure::voxelise() {
    
  double phiSpacing = 6.28318530718 / phiPoints;
  
  // generate points for charge density

  for (int i = 0; i < nTube; i++) {
    
    Tube tube = tubes[i];
    XYZ s = tube.s;
    XYZ t = tube.t;
    double r = tube.r;
    double l = tube.l;

    double zSpacing = l / (zPoints - 1);

    // generate points along surface of tube

    XYZ n;
    if (t.x != 0.0) n = XYZ(-(t.y + t.z)/t.x, 1.0, 1.0);
    else if (t.y != 0.0) n = XYZ(1.0, -(t.x + t.z)/t.y, 1.0);
    else n = XYZ(1.0, 1.0, -(t.x + t.y)/t.z);
    n.normalise();

    for (int j = 0; j < phiPoints; j++) {

      XYZ nRotated = rotate(n, t, j*phiSpacing);
      nRotated.normalise();
      
      for (int k = 0; k < zPoints; k++) {
        XYZ point = s + r*nRotated + k*zSpacing*t;
        
        // apply periodic boundary conditions
        
        if (point.x < 0.0) point.x += boxSize.x;
        if (point.y < 0.0) point.y += boxSize.y;
        if (point.z < 0.0) point.z += boxSize.z;
        if (point.x > boxSize.x) point.x -= boxSize.x;
        if (point.y > boxSize.y) point.y -= boxSize.y;
        if (point.z > boxSize.z) point.z -= boxSize.z;

        points.push_back(point);
      } 
    }
  }

  // calculate voxel sizes

  double voxelSizeX = boxSize.x / nVoxel.x;
  double voxelSizeY = boxSize.y / nVoxel.y;
  double voxelSizeZ = boxSize.z / nVoxel.z;

  voxelSize = XYZ(voxelSizeX, voxelSizeY, voxelSizeZ);
  
  // voxelise density

  if (saxs) {
    density2D = vector<vector<double>>(nVoxel.x, vector<double>(nVoxel.z, 0.0));

    for (int i = 0; i < points.size(); i++) {
      XYZ point = points[i];
      
      int indexX = floor(point.x / voxelSizeX);
      int indexZ = floor(point.z / voxelSizeZ);

      density2D[indexX][indexZ] += 1.0;
    }
  }
  else {

    density = vector<vector<vector<double>>>(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));

    for (int i = 0; i < points.size(); i++) {
      XYZ point = points[i];

      int indexX = floor(point.x / voxelSizeX);
      int indexY = floor(point.y / voxelSizeY);
      int indexZ = floor(point.z / voxelSizeZ);
    
      density[indexX][indexY][indexZ] += 1.0;
    }
  }
}

void Structure::fft() {

  // calculate reciprocal voxel sizes;
  
  double voxelFFTSizeX = 2.0 * numbers::pi / boxSize.x;
  double voxelFFTSizeY = 2.0 * numbers::pi / boxSize.y;
  double voxelFFTSizeZ = 2.0 * numbers::pi / boxSize.z;

  voxelFFTSize = XYZ(voxelFFTSizeX, voxelFFTSizeY, voxelFFTSizeZ);

  if (saxs) {
    
    // setup fft procedure
    
    vector<double> densityUnwrapped(int(nVoxel.x*nVoxel.z), 0.0);
    vector<complex<double>> amplitudeUnwrapped(int(nVoxel.x)*(int(nVoxel.z)/2 + 1), {0.0, 0.0});
    fftw_complex *fftwAmplitudeUnwrapped = reinterpret_cast<fftw_complex*>(&amplitudeUnwrapped[0]);

    fftw_plan plan = fftw_plan_dft_r2c_2d(int(nVoxel.x), int(nVoxel.z), &densityUnwrapped[0], fftwAmplitudeUnwrapped, FFTW_ESTIMATE);

    // unwrap density data to 1D vector

    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.z; j++)
        densityUnwrapped[i*int(nVoxel.z) + j] = density2D[i][j];

    // run fft and clean up

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // calculate intensity and wrap up into 2d vector

    intensity2D = vector<vector<double>>(nVoxel.x, vector<double>(nVoxel.z, 0.0));
    
    int nVoxelZTemp = int(nVoxel.z/2) + 1;

    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.z; j++) {
        int jTemp = j;
        if (j > int(nVoxel.z)/2) jTemp = int(nVoxel.z) - j - 1;
        int index = i*nVoxelZTemp + jTemp;
        intensity2D[i][j] = norm(amplitudeUnwrapped[index]);
      } 

    // wrap fft around and center at origin
    
    vector<vector<double>> intensityCentred(nVoxel.x, vector<double>(nVoxel.z, 0.0));
    
    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.z; j++) {
        int iTemp, jTemp, kTemp;
        if (i > nVoxel.x/2) iTemp = i - int(nVoxel.x)/2;
        else iTemp = i + int(nVoxel.x)/2 - 1;
        if (j > nVoxel.z/2) jTemp = j - int(nVoxel.z)/2;
        else jTemp = j + int(nVoxel.z)/2 - 1;
      
        intensityCentred[i][j] = intensity2D[iTemp][jTemp];
      }

    intensity2D = intensityCentred;

  }
  else {

    // setup fft procedure
    
    vector<double> densityUnwrapped(int(nVoxel.x*nVoxel.y*nVoxel.z), 0.0);
    vector<complex<double>> amplitudeUnwrapped(int(nVoxel.x*nVoxel.y)*(int(nVoxel.z)/2 + 1), {0.0, 0.0});
    fftw_complex *fftwAmplitudeUnwrapped = reinterpret_cast<fftw_complex*>(&amplitudeUnwrapped[0]);

    fftw_plan plan = fftw_plan_dft_r2c_3d(int(nVoxel.x), int(nVoxel.y), int(nVoxel.z), &densityUnwrapped[0], fftwAmplitudeUnwrapped, FFTW_ESTIMATE);

    // unwrap density data to 1D vector

    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.y; j++)
        for (int k = 0; k < nVoxel.z; k++)
          densityUnwrapped[i*int(nVoxel.y)*int(nVoxel.z) + j*int(nVoxel.z) + k] = density[i][j][k];

    // run fft and clean up

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // calculate intensity and wrap up into 3d vector

    intensity = vector<vector<vector<double>>>(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));
    
    int nVoxelZTemp = int(nVoxel.z/2) + 1;

    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.y; j++)
        for (int k = 0; k < nVoxel.z; k++) {
          int kTemp = k;
          if (k > int(nVoxel.z)/2) kTemp = int(nVoxel.z) - k - 1;
          int index = i*int(nVoxel.y)*nVoxelZTemp + j*nVoxelZTemp + kTemp;
          intensity[i][j][k] = norm(amplitudeUnwrapped[index]);
        } 

    // wrap fft around and center at origin
    
    vector<vector<vector<double>>> intensityCentred(nVoxel.x, vector<vector<double>>(nVoxel.y, vector<double>(nVoxel.z, 0.0)));
    
    for (int i = 0; i < nVoxel.x; i++)
      for (int j = 0; j < nVoxel.y; j++)
        for (int k = 0; k < nVoxel.z; k++) {
          int iTemp, jTemp, kTemp;
          if (i > nVoxel.x/2) iTemp = i - int(nVoxel.x)/2;
          else iTemp = i + int(nVoxel.x)/2 - 1;
          if (j > nVoxel.y/2) jTemp = j - int(nVoxel.y)/2;
          else jTemp = j + int(nVoxel.y)/2 - 1;
          if (k > nVoxel.z/2) kTemp = k - int(nVoxel.z)/2;
          else kTemp = k + int(nVoxel.z)/2 - 1;
        
          intensityCentred[i][j][k] = intensity[iTemp][jTemp][kTemp];
        }

    intensity = intensityCentred;
  }
}
