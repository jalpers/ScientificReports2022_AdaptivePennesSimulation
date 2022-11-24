#ifndef PENNESHEATTRANSFER_H
#define PENNESHEATTRANSFER_H
#include<Matrix3D.h>
#include<vtkImageData.h>
#include<vtkSmartPointer.h>
#include<vector>
#include<array>
#include<SimulationRessources.h>


//test stuff
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <cmath>

using namespace std::chrono;
using namespace std;
class PennesHeatTransfer
{
public:
	PennesHeatTransfer();
	PennesHeatTransfer(bool _multisteps);
	void setInitValues(Matrix3D<PennesEquationParameter>* _initValues);
	void setInitValues(vtkSmartPointer<vtkImageData> _initValues);
	Matrix3D<PennesEquationParameter>* getNewValues();
	Matrix3D<PennesEquationParameter>* getInitialValues();
	void modifyInitValues(const double* x, const int* _param, int n, float angle = 1000);

	~PennesHeatTransfer();
	void setUIConnection(vtkSmartPointer < vtkImageData> _data);
	void finiteStep(int _scheme = 0);
	void finiteStep_opti(int* _T = nullptr, int* _needleAxis = nullptr);
	void updateDataUi();

	void setTimeSteps(double _dt, int _N = 1);
	void saveTimeStep();
	Matrix3D<PennesEquationParameter>* getSavedTimeStep() { return P_saved; }
	void setVoxelConstant();
	void setVesselVoxels(Matrix3D<bool>& _map);
	void updateNeedleT(const vector<int> _T);
	void updateNeedleQ(const vector<double> _Q);

	Matrix3D<PennesEquationParameter>* P_n;
	Matrix3D<PennesEquationParameter>* P_new;
	Matrix3D<PennesEquationParameter>* P_saved;
private:


	void ADI_heatEquation();
	double error();
	bool is_onborder(int i, int j, int l);
	void thomas(vector<double>& r, double lamda);
	void printMatrix();
	void copyMatrix();
	void readData();
	vector<vector<vector<double> > > matrix3D, previousMatrix3D;
	//ifstream fin("ADI_input.txt");
	int x, y, z;
	int u0, k, M, Q;
	double precision;
	double lamdaX, lamdaY, lamdaZ;
	double deltaX, deltaY, deltaZ;




	void stepImplicitADI();
	void stepExplicitEuler();
	void thomas(double* a, double* b, double* c, double* d, int n);
	void debugViewer(Matrix3D<PennesEquationParameter>* _b);
	void debugViewer(Matrix3D<double>* _b);
	
	bool setCylinderCoord = false;

	vtkSmartPointer<vtkImageData> dataUI;


	double* dx_2;
	double dt;
	int N;

	bool multiStep;

	double avagHeat;

};
#endif // PENNESHEATTRANSFER_H