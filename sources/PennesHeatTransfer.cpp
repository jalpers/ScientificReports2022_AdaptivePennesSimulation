#include "..\include\PennesHeatTransfer.h"
#include <math.h>
#include<opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <omp.h>
#include<ParameterOptimization.h>
using namespace cv;
PennesHeatTransfer::PennesHeatTransfer()
{
	dt = 0.1;
	N = 1;
	P_n = new Matrix3D<PennesEquationParameter>();
	P_new = new Matrix3D<PennesEquationParameter>();
	P_saved = new Matrix3D<PennesEquationParameter>();
	dx_2 = new double[3];
	multiStep = false;
}
PennesHeatTransfer::PennesHeatTransfer(bool _ms)
{
	dt = 0.1;
	N = 1;
	P_n = new Matrix3D<PennesEquationParameter>();
	P_new = new Matrix3D<PennesEquationParameter>();
	P_saved = new Matrix3D<PennesEquationParameter>();
	dx_2 = new double[3];
	multiStep = _ms;
}
PennesHeatTransfer::~PennesHeatTransfer()
{
	delete P_n;
	delete P_new;
	delete P_saved;
	delete dx_2;
}
void PennesHeatTransfer::setTimeSteps(double _dt, int _N)
{
	dt = _dt;
	if (_N > 1)
	{
		//std::cout << "multistep" << std::endl;
		multiStep = true;
	}

	N = _N;
}
void PennesHeatTransfer::setInitValues(Matrix3D<PennesEquationParameter>* _initValues)
{
	
	int* dim = _initValues->getDimensions();
	double* spacing = _initValues->getSpacing();
	dx_2 = new double[3];
	for (int i = 0; i < 3; i++)
	{
		dx_2[i] = spacing[i] * spacing[i] * 0.001 * 0.001 ;/* 0.001 * 0.001*/
	}

	P_n = new Matrix3D<PennesEquationParameter>(*_initValues);
	P_new = new Matrix3D<PennesEquationParameter>(dim[0], dim[1], dim[2], setCylinderCoord);
	P_new->setSpacing(spacing);
	PennesEquationParameter pn;
	if (GlobalVariables::boundaryCondition == 1)
	{
		P_n->setAmbientValue(pn);
		P_new->setAmbientValue(pn);
		//std::cout << "Default: " << (*P_n)[{-1, 0, 0}].T << std::endl;
	}
	//Debug-----------------------------
	/*long double Heatsum = 0;
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
				Heatsum += (*P_n)[{i, j, k}].T;
			}
		}
	}
	std::cout << "StartHeat: " << Heatsum << " " << Heatsum / (dim[0] * dim[1] * dim[2]) << std::endl;*/
	//std::cout << "DeltaX_2: " << dx_2[0] <<" "<< dx_2[1]<< " "<<dx_2[2]<< std::endl;
	//std::cout << "Dimension: " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
	//Debug-----------------------------
	

}
void PennesHeatTransfer::setInitValues(vtkSmartPointer<vtkImageData> _initialTValues)
{
	int* dim = new int[3];
	for (int i = 0; i < 3; i++)
	{
		dim[i] = _initialTValues->GetDimensions()[i];
	}

	//dim[0] = dim[1] = dim[2] = 60;
	//dim[1] = 1;


	Matrix3D<PennesEquationParameter>* initValues = new Matrix3D<PennesEquationParameter>(dim[0], dim[1], dim[2], setCylinderCoord);

	initValues->setSpacing(_initialTValues->GetSpacing());
	


	double Heatsum = 0;
	PennesEquationParameter p;
	for( int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
				p = PennesEquationParameter();

				(*initValues)[{i, j, k}] = p;
				//(*initValues)[{(dim[0] / 2), j, (dim[2] / 2)}].T =21.0;
				
				/*if(j<25)
					(*initValues)[{( dim[0] / 2), j,(dim[2] / 2)}].T =exp(.1 * j)*10.0;
				if(j>=25 && j<35)
					(*initValues)[{(dim[0] / 2), j, (dim[2] / 2)}].T =  10.0;*/


				Heatsum += (*initValues)[{i, j, k}].T;
			}
		}
	}

	//(*P_n)[{size_t(dim[0] / 2), size_t(dim[1] / 2), size_t(dim[2] / 2)}].Q_r = 100000.0;

	std::cout << "StartHeat: " << Heatsum<< " "<<Heatsum / (dim[0] * dim[1] * dim[2])<<std::endl;
	std::cout << "Pnnesconstr-done" << endl;

	setInitValues(initValues);

}



void PennesHeatTransfer::setUIConnection(vtkSmartPointer<vtkImageData> _data)
{
	dataUI = _data;
}
void PennesHeatTransfer::finiteStep_opti(int* _T, int* _needleAxis)
{
	std::cout << "Start finiteStep_opti" << std::endl;
	float devHeat;
	for (int i = 0; i < N; i++)
	{
		
		if (_T != nullptr)
		{

			for (int j = 0; j < 60; j++)
			{

				devHeat = (*P_new)[{ j, _needleAxis[1], 0}].T - (float)_T[j];

				//std::cout << "oldT: " << (int)(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " ";

				(*P_new)[{ j, _needleAxis[1], 0}].T -= (int)devHeat;

				devHeat = (*P_n)[{ j, _needleAxis[1], 0}].T - (float)_T[j];

				//std::cout << "oldT: " << (int)(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " ";

				(*P_n)[{ j, _needleAxis[1], 0}].T -= (int)devHeat;
				/*(*P_n)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].Q_r += (int)(devHeat *
					((*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].rho*
						(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].c));*/

						//	std::cout <<"newT: " <<(int)(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " devheat " << devHeat << " talong " << (float)TalongNeedle[y] << std::endl;
			}
		}
		stepImplicitADI();
	}
	
}
void PennesHeatTransfer::finiteStep(int _scheme)
{
	//std::cout << "finiteStep" << std::endl;
	//stepExplicitEuler();
	double q,t;
	int* dim = P_n->getDimensions();
	double Qsum = 0;
	long double heatSumIs = 0;
	double heatSumNew = 0;
	//for (int i = 0; i < dim[0]; i++)
	//{
	//	for (int j = 0; j < dim[1]; j++)
	//	{
	//		for (int k = 0; k < dim[2]; k++)
	//		{
	//			heatSumIs += (*P_n)[{i, j, k}].T;
	//		}
	//	}
	//}
	//std::cout << "HeatsumIs: " << heatSumIs / (dim[0] * dim[1] * dim[2])<< std::endl;

	//std::cout << "Heating Temperatur: " << std::endl;
	//
	//for (int y = 0; y < P_n->getDimensions()[1]; y++)
	//{
	//	PennesEquationParameter* v_n = &(*P_n)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}];
	//	v_n->T = GlobalVariables::baseLineT + v_n->Q_r_rel * (v_n->T_ind_max - GlobalVariables::baseLineT);
	//	std::cout << (*P_n)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "scheme: " << _scheme << std::endl;
	for (int i = 0; i < N; i++)
	{
		if (_scheme == 0)
		{
			std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
			stepImplicitADI();
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << std::endl;

			//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds > (end - begin).count() << "[ms]" << std::endl;
			
		}	
		if (_scheme == 1)
		{
			
			stepExplicitEuler();
			
		}
		
	}
}
void PennesHeatTransfer::updateNeedleQ(const vector<double> _Q)
{
	std::cout << "updateNeedleQ" << std::endl;
	ParameterOptimization::adjustParameterNeedleAxis(_Q, P_n);
}
void PennesHeatTransfer::updateNeedleT(const vector<int> _T)
{
	std::cout << "updateNeedleT" << std::endl;
	ParameterOptimization::adjustParameterNeedleAxis(_T, P_n);
}
void PennesHeatTransfer::setVoxelConstant()
{
	for (int y = 0; y < P_n->getDimensions()[1]; y++)
	{
		int pON[3] = { GlobalVariables::pointOnNeedle_volume[0], y ,GlobalVariables::pointOnNeedle_volume[2] };
		std::cout << "Pon" << pON[0] << " " << pON[1] << " " << pON[2] << std::endl;
		(*P_n)[pON].keepTconstant = true;
		(*P_new)[pON].keepTconstant = true;
	}
}
void PennesHeatTransfer::setVesselVoxels(Matrix3D<bool>& _maps)
{
	int* dim = _maps.getDimensions();
	double T_a = 25;
	double w_b = 800 / (1e6 * 60) * 997 / (0.005 * 0.005 * 1 * 3.1415);
	double c_b = 4182;//
	for (int x = 0; x < dim[0]; x++)
	{
		for (int y = 0; y < dim[1]; y++)
		{
			for (int z = 0; z < dim[2]; z++)
			{
				
				if (_maps[{x, y, z}])
				{
					(*P_n)[{x, y, z}].T_a = T_a;
					(*P_new)[{x, y, z}].T_a = T_a;

					(*P_n)[{x, y, z}].w_b = w_b;
					(*P_new)[{x, y, z}].w_b = w_b;

					(*P_n)[{x, y, z}].c_b= c_b;
					(*P_new)[{x, y, z}].c_b = c_b;

					(*P_n)[{x, y, z}].isVessel = true;
					(*P_new)[{x, y, z}].isVessel = true;
				}

			}
		}
	}
}




void PennesHeatTransfer::updateDataUi()
{
	//TODO Exchange data, das beides auf dem gleiche stand ist
	std::cout << "updateDataUi" << std::endl;
	double sumHeat = 0;
	int* dim = P_new->getDimensions();
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
				

				unsigned char* v = static_cast<unsigned char*>(dataUI->GetScalarPointer(i, j, k));
				unsigned char l = (unsigned char)(int)std::round((*P_new)[{i, j, k}].T);
				v[0] = l;
				/*if(i==GlobalVariables::pointOnNeedle_volume[0] && k == GlobalVariables::pointOnNeedle_volume[2])
					std::cout <<" nadle: " <<(*P_new)[{i, j, k}].T<<" "<< (*P_n)[{i, j, k}].T ;*/
				//sumHeat += (*P_new)[{i, j, k}].T;
		

			}
			/*if (i == (int)dim[1] / 2)
				std::cout << std::endl;*/
		}
		
	}
	dataUI->Modified();
	//std::cout << std::endl;
	//avagHeat = sumHeat / (dim[0] * dim[1] * dim[2]);
	//std::cout <<"sumHeat<: "<< sumHeat<<" avagHeat: " << avagHeat << std::endl;
}
void PennesHeatTransfer::stepExplicitEuler()
{
	//std::cout << "StepExplicitEuler" << std::endl;
	long double Heatsum = 0;
	int* dim = P_n->getDimensions();
	Matrix3D<double>* T = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	Matrix3D<double>* T_new = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	double res;
	double heatflux = 0;
	double td;
	
	double defaultValue = 21;
	if (GlobalVariables::boundaryCondition == 1)
	{
		T->setAmbientValue(defaultValue);
		T_new->setAmbientValue(defaultValue);
	}

#pragma omp parallel for collapse(3)
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
				
				PennesEquationParameter v_n = (*P_n)[{i, j, k }];
				(*T)[{i, j, k}] = v_n.T;
				(*P_new)[{i, j, k}] = v_n;
				double hf_b = v_n.w_b * v_n.c_b * (v_n.T_a - v_n.T);
				//std::cout <<"hf_b "<< hf_b;
				//double hf_t = dt * (hf_b + v_n.Q_m + v_n.Q_r) / (v_n.rho * v_n.c) ;
				double hf_t = dt * (hf_b + v_n.Q_m + v_n.Q_r) / (v_n.M);
				if (v_n.keepTconstant == true)
				{
					(*P_new)[{i, j, k}].T = v_n.T;
					(*T_new)[{i, j, k}] = v_n.T;
				}
				else
				{
					td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
					heatflux = 
						(dt * td * ((*P_n)[{i-1, j, k }].T + (*P_n)[{i + 1, j, k }].T - 2 * v_n.T))/dx_2[0]
						+ (dt * td * ((*P_n)[{i, j - 1, k }].T + (*P_n)[{i, j + 1, k }].T - 2 * v_n.T)) / dx_2[1]
						+ (dt * td * ((*P_n)[{i, j, k - 1 }].T + (*P_n)[{i, j, k + 1 }].T - 2 * v_n.T)) / dx_2[2];

					res = v_n.T + heatflux;
					
					double plu = (
						(*P_n)[{i - 1, j, k }].T + (*P_n)[{i + 1, j, k }].T
						+ (*P_n)[{i, j - 1, k }].T + (*P_n)[{i, j + 1, k }].T
						+ (*P_n)[{i, j, k - 1 }].T + (*P_n)[{i, j, k + 1 }].T
						) / dx_2[0];
					
					res = (1 - (6 * dt * td / dx_2[0])) * v_n.T +
						dt * td * plu;
					//cout << plu << " " << res << std::endl;

					

					//std::cout << " res: " << res << " heatflux " << heatflux;
					(*P_new)[{i, j, k}].T = res;
					(*T_new)[{i, j, k}] = res;
				}
				
			}
		}
	}

	if (multiStep)
	{
		for (int i = 0; i < dim[0]; i++)
		{
			for (int j = 0; j < dim[1]; j++)
			{
				for (int k = 0; k < dim[2]; k++)
				{

 					(*P_n)[{i, j, k}].T = (*P_new)[{i, j, k}].T;
					Heatsum += abs(defaultValue - (*P_new)[{i, j, k}].T);
					/*if ((*P_new)[{i, j, k}].T <= 0.1)
					{
						std::cout <<"null: " << i << " " << " " << j << " " << k << std::endl;
					}*/

				}
			}
		}
	}
	//std::cout << "-------T.............." << std::endl;
	//ParameterOptimization::debugViewer(T);
	//std::cout << "-------T_new.............." << std::endl;
	//ParameterOptimization::debugViewer(T_new);
	delete T;
	delete T_new;
	//std::cout << "StepExplicitEuler - Heatsum: " << Heatsum << endl;// / (P_n->getDimensions()[0] * P_n->getDimensions()[1] * P_n->getDimensions()[2]) << std::endl;
}
void PennesHeatTransfer::stepImplicitADI()
{
	//std::cout << "StepImlicitADI-Start" << std::endl;

	//crank-Nicolson - von nuemann adiabatisch randbedingung
	
	int* dim = P_n->getDimensions();

	//std::cout << "T:mitte: " << (*P_n)[{dim[0] / 2, dim[1] / 2, dim[2] / 2}].T << std::endl;

	//Matrix3D<double>* T = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	Matrix3D<double>* T_new = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	Matrix3D<double>* T_s = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	Matrix3D<double>* T_ss = new Matrix3D<double>(dim[0], dim[1], dim[2]);
	
	double defaultValue = 21;
	if (GlobalVariables::boundaryCondition == 1)
	{
	//	T->setAmbientValue(defaultValue);
		T_s->setAmbientValue(defaultValue);
		T_ss->setAmbientValue(defaultValue);
		T_new->setAmbientValue(defaultValue);
	}

	//Heating
	//std::cout << "Heating Temperatur: " << std::endl;
	#pragma omp parallel for 
	for (int y = 0; y < P_n->getDimensions()[1]; y++)
	{
		PennesEquationParameter* v_n = &(*P_n)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}];
		if (v_n->Q_r_rel != 0)
		{
			v_n->T = GlobalVariables::baseLineT + v_n->Q_r_rel * (v_n->T_ind_max - GlobalVariables::baseLineT);
		}
	
		//std::cout << (*P_n)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " ";
	}
	//std::cout << std::endl;
	//#pragma omp parallel for collapse(3)
	//for (int i = 0; i < dim[0]; i++)
	//{
	//	for (int j = 0; j < dim[1]; j++)
	//	{
	//		for (int k = 0; k < dim[2]; k++)
	//		{

	//			(*T)[{i, j, k}] = (*P_n)[{i, j, k}].T;
	//			/*if ((*T)[{i, j, k}] > 22)
	//				std::cout << i << j << k << std::endl;*/

	//		}
	//	}
	//}
	/*std::cout << (*P_n)[{-1, 0, 0}].T << " " << (*P_n)[{0, -1, 0}].T << " " << (*P_n)[{0, 0, -1}].T
		<< " " << (*P_n)[{60, 0, 0}].T << " " << (*P_n)[{0, dim[1], 0}].T << " " << (*P_n)[{0, 0, 60}].T << std::endl;
	std::cout<<"TEst: "<< (*P_n)[{0, dim[1], 0}].T << " " << (*P_n)[{0, 0, 60}].T << std::endl;*/
	//x-direction	

	#pragma omp parallel for collapse(2)
	for (int j = 0; j < dim[1]; j++)
	{
		for (int k = 0; k < dim[2]; k++)
		{


			double *a = new double[dim[0]]();
			double *b = new double[dim[0]]();
			double* c = new double[dim[0]]();
			double *d = new double[dim[0]]();

			double r1, r2, r3;
			double td; //thermal diffusivity (m^2/s)
			
			for (int i = 0; i < dim[0]; i++)
			{



				//TODO k flux conservative 
					//TODO boundary conditition
				//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
				//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
				td = (*P_n)[{i, j, k}].getTd();
				r1 = r2 = r3 = td * dt / dx_2[0];

				a[i] = -r1 / 2;
				b[i] = 1 + r1;
				c[i] = -r1 / 2;


				d[i] = r1 / 2 * (*P_n)[{i - 1, j, k}].T + r1 / 2 * (*P_n)[{i + 1, j, k}].T
					+ r2 * (*P_n)[{i, j - 1, k}].T + r2 * (*P_n)[{i, j + 1, k}].T
					+ r3 * (*P_n)[{i, j, k - 1}].T + r3 * (*P_n)[{i, j, k + 1}].T
					+ (1 - r1 - 2 * r2 - 2 * r3) * (*P_n)[{i, j, k}].T;

			}
			a[0] = 0;
			c[dim[0] - 1] = 0;
			if (GlobalVariables::boundaryCondition == 0)//neumann -gleiche t_-1 = t_0
			{
				b[0] -= r1 / 2;
				b[dim[0] - 1] -=  r1 / 2;
			}
			else
			{
				
				d[0] += r1 / 2 * defaultValue;
				
				d[dim[0] - 1] += r1 / 2 * defaultValue;
				/*d[0] += r3 / 2 * defaultValue + r2 / 2 * defaultValue +  r1 / 2 * defaultValue;
				d[dim[0] - 1] += r3 / 2 * defaultValue + r2 / 2 * defaultValue +  r1 / 2 * defaultValue;*/
			}
			//d[0] = d[0] + a[0] * (*P_n)[{0, j, k}].T;
			//d[dim[0] - 1] = d[dim[0] - 1] + c[0] * (*P_n)[{dim[0] - 1, j, k}].T;
			thomas(a, b, c, d, dim[0]);
			
			for (int l = 0; l < dim[0]; l++)
			{
				/*if (d[l] < 20)
					std::cout << "dx :" << l << " " << j << " " << k << " : " << d[l] << std::endl;*/
				(*T_s)[{l, j, k}] = d[l];
			}
			delete[] a;
			delete[] b;
			delete[] c;
			delete[] d;
		}
	}




	
	//y-Direction
	
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < dim[0]; i++)
	{
		for (int k = 0; k < dim[2]; k++)
		{
			double* a = new double[dim[1]]();
			double* b = new double[dim[1]]();
			double* c = new double[dim[1]]();
			double* d = new double[dim[1]]();
			double r1, r2, r3;
			double td; //thermal diffusivity (m^2/s)
			
			for (int j = 0; j < dim[1]; j++) 
			{
				
				//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
				//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
				td = (*P_n)[{i, j, k}].getTd();
				r1 = r2 = r3 = td * dt / dx_2[1];

				a[j] = -r2 / 2;
				b[j] = 1 + r2;
				c[j] = -r2 / 2;

				d[j] = r1 / 2 * (*P_n)[{i - 1, j, k}].T + r1 / 2 * (*P_n)[{i + 1, j, k}].T 
					+ r1 / 2 * (*T_s)[{i - 1, j, k}] + r1 / 2 * (*T_s)[{i + 1, j, k}] 
					+ r2 / 2 * (*P_n)[{i, j - 1, k}].T + r2 / 2 * (*P_n)[{i, j + 1, k}].T
					+ r3 * (*P_n)[{i, j, k - 1}].T + r3 * (*P_n)[{i, j, k + 1}].T
					- r1 * (*T_s)[{i, j, k}] 
					+ (1 - r1 - r2 - 2 * r3) * (*P_n)[{i, j, k}].T;
			}
			a[0] = 0;
			c[dim[1] - 1] = 0;
			//cout << "d[0]_y " << d[0] << std::endl << std::endl;
			if (GlobalVariables::boundaryCondition == 0 )//neumann -gleiche t_-1 = t_0
			{
				b[0] -= r2 / 2;
				b[dim[1] - 1] -= r2 / 2;
			}
			else
			{
				d[0] += r2 / 2 * defaultValue;
				d[dim[1] - 1] += r2 / 2 * defaultValue;
				/*d[0] += r3 / 2 * defaultValue + r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;
				d[dim[1] - 1] += r3 / 2 * defaultValue + r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;*/
			}
			
			
			thomas(a, b, c, d, dim[1]);
			//cout << "d[0]_y " << d[0] << std::endl << std::endl;
			
			for (int l = 0; l < dim[1]; l++)
			{
				/*if (d[l] < 20)
					std::cout << "dy :" << i << " " << l << " " << k << " : " << d[l] << std::endl;*/
				/*if ((*P_n)[{i, l, k}].keepTconstant)
					(*T_ss)[{i, l, k}] = (*T_s)[{i, l, k}];
				else*/
					(*T_ss)[{i, l, k}] = d[l];
			}
			delete[] a;
			delete[] b;
			delete[] c;
			delete[] d;
		}
	}



	//z.direction
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			double *a = new double[dim[2]]();
			double *b = new double[dim[2]]();
			double *c = new double[dim[2]]();
			double *d = new double[dim[2]]();
			double r1, r2, r3;
			double td; //thermal diffusivity (m^2/s)
			double res;
			
			for (int k = 0; k < dim[2]; k++)
			{
				
				//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
				//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
				td = (*P_n)[{i, j, k}].getTd();
				r1 = r2 = r3 = td * dt / dx_2[2];
				
				a[k] = -r3 / 2;
				b[k] = 1 + r3;
				c[k] = -r3 / 2;
				
				d[k] = r1 / 2 * (*P_n)[{i - 1, j, k}].T + r1 / 2 * (*P_n)[{i + 1, j, k}].T
					+ r1 / 2 * (*T_s)[{i - 1, j, k}] + r1 / 2 * (*T_s)[{i + 1, j, k}]
					+ r2 / 2 * (*P_n)[{i, j - 1, k}].T + r2 / 2 * (*P_n)[{i, j + 1, k}].T
					+ r2 / 2 * (*T_ss)[{i, j - 1, k}] + r2 / 2 * (*T_ss)[{i, j + 1, k}] 
					+ r3 / 2 * (*P_n)[{i, j, k - 1}].T + r3 / 2 * (*P_n)[{i, j, k + 1}].T
					- r1 * (*T_s)[{i, j, k}] - r2 * (*T_ss)[{i, j, k}] 
					+ (1 - r1 - r2 - r3) * (*P_n)[{i, j, k}].T;
		/*		d[k] = r1 / 2 * (*T)[{i - 1, j, k}] + r1 / 2 * (*T)[{i + 1, j, k}]
					+ r1 / 2 * (*T_s)[{i - 1, j, k}] + r1 / 2 * (*T_s)[{i + 1, j, k}]
					+ r2 / 2 * (*T)[{i, j - 1, k}] + r2 / 2 * (*T)[{i, j + 1, k}]
					+ r2 / 2 * (*T_ss)[{i, j - 1, k}] + r2 / 2 * (*T_ss)[{i, j + 1, k}] 
					+ r3 / 2 * (*T)[{i, j, k - 1}] + r3 / 2 * (*T)[{i, j, k + 1}]
					- r1 * (*T_s)[{i, j, k}] - r2 * (*T_ss)[{i, j, k}] 
					+ (1 - r1 - r2 - r3) * (*T)[{i, j, k}];*/

			}
			//cout <<"d[0]_z "<< d[0]<<std::endl;
			a[0] = 0;
			c[dim[2] - 1] = 0;
		
			if (GlobalVariables::boundaryCondition == 0)//neumann -gleiche t_-1 = t_0
			{

				b[0] -= r3 / 2;
				b[dim[2] - 1] -= r3 / 2;
			}
			else
			{
				d[0] += r3 / 2 * defaultValue;
				//cout << "d[0] " << d[0] << std::endl;
				d[dim[2] - 1] += r3 / 2 * defaultValue;
				/*d[0] += r3 / 2 * defaultValue + 2 * r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;
				d[dim[2]-1] += r3 / 2 * defaultValue + 2 * r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;*/
			}
			
			
			//cout << "d[0] " << d[0] << std::endl;
			thomas(a, b, c, d, dim[2]);
			//cout << "d[0]_z " << d[0] << std::endl << std::endl;
			
			for (int l = 0; l < dim[2]; l++)
			{
				PennesEquationParameter v_n = (*P_n)[{i, j, l}];
				(*P_new)[{i, j, l}] = v_n;
				double hf_b = v_n.w_b * v_n.c_b * (v_n.T_a - v_n.T);
				//std::cout <<"hf_b "<< hf_b;
				//double hf_t = dt * (hf_b + v_n.Q_m + v_n.Q_r) / (v_n.rho * v_n.c) ;
				double hf_t = dt * (hf_b + v_n.Q_m + v_n.Q_r) / (v_n.M);
				if (d[i] > 100 || d[i] < 0)
					{
					std::cout << "Falsch d pen " << d[i];
					}
				/*if(hf_t > 0.1|| hf_t < -0.1)
					std::cout  <<"hft_t: " <<hf_t << std::endl;*/
				/*if (d[l] < 20)
					std::cout<<"d :" << i << " " << j << " " << l <<" : "<<d[l]<<  std::endl;*/
				/*if (hf_t> 0)
					std::cout <<"hft: "<<hf_t<<" " << i << " " << j << " " << l << std::endl;*/
				/*if (d[l] > 100)
					std::cout << i<< " " << j << " " <<l << std::endl;*/
				//int v;
				//if (v_n.isVessel)
				//{
				//	v = 1;
				//	//std::cout << "isVessel" << std::endl;
				//}
				//	
				//else
				//	v = 0;
				//Pt = rho_b * w * Cb / (rho * Ct) * (T_amb - T_new)
				double P_t = v_n.w_b * v_n.c_b * (v_n.T_a - d[l]) * (v_n.getTd() / v_n.k);
				/*if (P_t > 0)
					std::cout <<"P_t" << P_t << " Diff "<< (v_n.T_a - d[l])<<std::endl;*/
				res = d[l] + dt * P_t ;//+ hf_t   /((v_n.getTd()/v_n.k))
				//+dt * (v_n.w_b * v_n.c_b * (v_n.T_a - v_n.T))
				/*if (v_n.keepTconstant == true)
				{
					(*P_new)[{i, j, l}].T = v_n.T;
					(*T_new)[{i, j, l}] = v_n.T;
					
				}
				else*/
				{
					(*P_new)[{i, j, l}].T = res;
					(*T_new)[{i, j, l}] = res;
				}
				
			}
			delete[] a;
			delete[] b;
			delete[] c;
			delete[] d;
		}
	}


	long double Heatsum = 0;

	if (multiStep)
	{
		std::cout << "multi" << std::endl;
		#pragma omp parallel for collapse(3)
		for (int i = 0; i < dim[0]; i++)
		{
			for (int j = 0; j < dim[1]; j++)
			{
				for (int k = 0; k < dim[2]; k++)
				{ 

					(*P_n)[{i, j, k}].T = (*P_new)[{i, j, k}].T;
					#pragma omp atomic
					Heatsum += (*P_new)[{i, j, k}].T;
					//if ((*P_new)[{i, j, k}].T <= 20 || (*P_new)[{i, j, k}].T >= 22)
					//{
					//	std::cout <<"null: " << i << " " << " " << j << " " << k <<" "<< (*P_new)[{i, j, k}].T << std::endl;
					//}

				}
			}
		}
	}
	
	//debugViewer(T_s);
	//ParameterOptimization::debugViewer(T_ss);
	//ParameterOptimization::debugViewer(P_new);


	//delete T;
	delete T_s;
	delete T_ss;
	delete T_new;
	//std::cout << "stepImplicitADI - Heatsum: "<<Heatsum<<" " << (double)(Heatsum / (P_n->getDimensions()[0] * P_n->getDimensions()[1] * P_n->getDimensions()[2])) << std::endl;
	//debugViewer(P_new);
	//std::cout << "stepImplicitADI-End" << std::endl;



}




void PennesHeatTransfer::thomas(double* a, double* b, double* c, double* d, int n) {
	/*
	// n is the number of unknowns

	|b0 c0 0 ||x0| |d0|
	|a1 b1 c1||x1|=|d1|
	|0  a2 b2||x2| |d2|

	1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

		x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

	2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
		from 1st it.: -| a1x0 + a1g0x1        = a1r0
					-----------------------------
						  (b1 - a1g0)x1 + c1x2 = d1 - a1r0

		x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

	3rd iteration:      | a2x1 + b2x2   = d2
		from 2nd it. : -| a2x1 + a2g1x2 = a2r2
					   -----------------------
					   (b2 - a2g1)x2 = d2 - a2r2
		x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
	Finally we have a triangular matrix:
	|1  g0 0 ||x0| |r0|
	|0  1  g1||x1|=|r1|
	|0  0  1 ||x2| |r2|

	Condition: ||bi|| > ||ai|| + ||ci||

	in this version the c matrix reused instead of g
	and             the d matrix reused instead of r and x matrices to report results
	Written by Keivan Moradi, 2014
	*/
	n--; // since we start from x0 (not x1)
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {
		
		d[i] -= c[i] * d[i + 1];
	}
}

Matrix3D<PennesEquationParameter>* PennesHeatTransfer::getNewValues()
{
	return P_new;
}
Matrix3D<PennesEquationParameter>* PennesHeatTransfer::getInitialValues()
{
	return P_n;
}
void PennesHeatTransfer::modifyInitValues(const double* x, const int* _param, int _Np, float angle)
{
	int* dim = P_n->getDimensions();
	//std::cout << dim[0] << " " << dim[1] << " " << dim[2] <<" x0 "<< x[0]<< std::endl;
	int index = 0;
	long double Tdsum = 0;
	long double heatSumIs = 0;
	float diff;
	#pragma omp parallel for collapse(4)
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
		
				for (int p = 0; p < _Np; p++)
				{
					//if (!((*P_n)[{i, j, k}]).keepTconstant)
					{
						switch (_param[p])
						{

						case ParameterFlags::td:
							(*P_n)[{i, j, k}].td = x[p];//TODO noch abhöngig von T machen
							break;

						case ParameterFlags::k_a:
							(*P_n)[{i, j, k}].k_a = x[p];
							break;
						case ParameterFlags::k_b:
							(*P_n)[{i, j, k}].k_b = x[p];
							break;
						case ParameterFlags::k:
							(*P_n)[{i, j, k}].k = x[p];//TODO noch abhöngig von T machen
							break;
						case ParameterFlags::td_a:
						{
							//P_n->getAngle({ i, j, k });
							diff = std::abs(angle - P_n->getAngle({ i, j, k }));
							/*if (angle != 1000 && diff <= 22.5)
							{
								(*P_n)[{i, j, k}].td_a = (1 - diff / 22.5) * (*P_n)[{i, j, k}].td_a + (diff / 22.5) * x[p];
							}
							else*/
								(*P_n)[{i, j, k}].td_a = x[p];

							break;
						}
						case ParameterFlags::td_b:
							(*P_n)[{i, j, k}].td_b = x[p];
							break;
						case ParameterFlags::td_c:
							(*P_n)[{i, j, k}].td_c = x[p];
							break;
						case ParameterFlags::T_ind_max:
						{
							(*P_n)[{i, j, k}].T_ind_max = x[p];
							break;
						}
						case ParameterFlags::c_b:
						{
							if((*P_n)[{i, j, k}].isVessel)
								(*P_n)[{i, j, k}].c_b = x[p] * 1e3;
							break;

						}
						default:
							break;
						}
					}

					
				}
				/*heatSumIs += (*P_n)[{i, j, k}].T;
				Tdsum += (*P_n)[{i, j, k}].getTd();*/
			}
		}
	}
	//std::cout <<" Td_avag: " << Tdsum/ (dim[0] * dim[1] * dim[2]) << std::endl;
}

void PennesHeatTransfer::saveTimeStep()
{
	P_saved = new Matrix3D<PennesEquationParameter>(*P_n);
}