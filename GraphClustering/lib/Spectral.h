#ifndef AROYA_SPECTRAL_H
#define AROYA_SPECTRAL_H


#include"eigen3/Eigen/Dense"
using namespace Eigen;
#include"armadillo/include/armadillo"
using namespace arma;
#include<vector>
#include"KMeans.h"
using namespace std;

enum SpectralType{RatioCut,NormalizedCut};

class Spectral {
public:
	Spectral();
	void setSpectralType(const SpectralType &);
	void setPointData(const vector<vector<double>> &newData);
	void setAdjData(const vector<vector<double>> &newData);
	void run();
	void writeFile(const char *fileName);
	void setKEigenVectors(const int &);
	vector<int>getFlag();
	//clustering method
	AroyaKMeans myans;
private:
	//store origin data
	MatrixXd data;
	MatrixXd D, W, L;
	//store clustering ans
	vector<vector<double>>origin;
	vector<int>flag;
	//setting of Spectral Type
	SpectralType spectraltype;
	//anti-span the dimension
	int dimension;
	//record
	double maxDistance;
	void RatioCut();
	void NormalizedCut();
};

#endif
#include"Spectral.cpp"