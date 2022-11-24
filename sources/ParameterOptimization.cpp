#include "ParameterOptimization.h"

//#include<opencv2/opencv.hpp>
#include <thread>
#include <mutex>
#include<Debugger.h>
//PennesHeatTransfer* ParameterOptimization::pht = new PennesHeatTransfer(false);
double ParameterOptimization::duration = 0;
int  ParameterOptimization::T[60] = { 0 };
Matrix3D<double>* ParameterOptimization::mTSoll = new Matrix3D<double>();


Matrix3D<PennesEquationParameter>* ParameterOptimization::mIs = new Matrix3D<PennesEquationParameter>();
double ParameterOptimization::squaredError[GlobalVariables::numberOfSlices] = { 0 };
std::vector<double>* ParameterOptimization::measValues = new std::vector<double>();
std::vector<std::array<int, 3>>*ParameterOptimization::measPositionInMatrix = new std::vector<std::array<int, 3>>();
std::vector<std::array<int, 3>>* ParameterOptimization::slicePositionInMatrix = new std::vector<std::array<int, 3>>();
Parameters ParameterOptimization::simPar = Parameters();
Parameter2Optimize ParameterOptimization::par2opti = Parameter2Optimize();
PennesHeatTransfer_Gpu* ParameterOptimization::pht = new PennesHeatTransfer_Gpu();
//float* ParameterOptimization::wb = new float();
// Namespace nullifies the use of cv::function(); 
using namespace std;


ParameterOptimization::ParameterOptimization()
{
	
	PennesEquationParameter pn;
	double s[GlobalVariables::numberOfOptiParam];
	s[0] = 1;
	s[1] = 1;
//	s[2] = 1;
	//s[2] = 1;
	//s[3] = 1;
	/*s[2] = 1;
	s[3] = 1;*/
	double x[GlobalVariables::numberOfOptiParam];
	x[0] = GlobalVariables::baseLineT;
	x[1] = 1.5;
	//x[2] = 4182;
	//x[2] = 0.1;
	//x[3] = 0.01;
//	x[2] = pn.c_b / 1e3;
	/*x[2] = pn.td_b;
	x[3] = pn.td_c;*/
	//cout << "PN.td " << pn.td <<" " << pn.td_a<<" " <<pn.td_b<< std::endl;
	for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
	{
		x_is[i].setcontent(GlobalVariables::numberOfOptiParam, x);
		s_is[i].setcontent(GlobalVariables::numberOfOptiParam, s);
	}
	//std::cout << "initalParamToOpti: " << x_is[0][0] << std::endl;
	bestSlice = 0;
    
	
}
ParameterOptimization::~ParameterOptimization()
{
	//delete pht;
	delete mTSoll;
	delete mIs;
	deleteStuff();
}
ParameterOptimization::ParameterOptimization(PennesHeatTransfer_Gpu& _pht)
{
	
	double s[GlobalVariables::numberOfOptiParam];
	s[0] = 1;
	s[1] = 1;
	//s[2] = 1;
	//s[3] = 1;
	/*s[2] = 1;
	s[3] = 1;*/
	double x[GlobalVariables::numberOfOptiParam];
	x[0] = GlobalVariables::baseLineT;
	x[1] = 1.5;
//	x[1] = 4182;
	//exp
	//x[2] = 0.1;
	//x[3] = 0.01;
	//	x[2] = pn.c_b / 1e3;
		/*x[2] = pn.td_b;
		x[3] = pn.td_c;*/
		//cout << "PN.td " << pn.td <<" " << pn.td_a<<" " <<pn.td_b<< std::endl;
	for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
	{
		x_is[i].setcontent(GlobalVariables::numberOfOptiParam, x);
		s_is[i].setcontent(GlobalVariables::numberOfOptiParam, s);
	}
	//std::cout << "initalParamToOpti: " << x_is[0][0] << std::endl;
	bestSlice = 0;

	pht = new PennesHeatTransfer_Gpu(_pht);

}
void fctOptParam_jacc(const real_1d_array& x, real_1d_array& fi, real_2d_array& jac, void* ptr)
{
	std::cout << "fctOptParam" << std::endl;

	double x_[GlobalVariables::numberOfOptiParam];
	for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
	{
		x_[i] = x[i];
	}
	//std::cout << std::endl;
	//ParameterOptimization::debugViewer(ParameterOptimization::mTSoll);
	//ParameterOptimization::debugViewer(ParameterOptimization::mIs);
	PennesHeatTransfer pht = PennesHeatTransfer(true);
	double dt = 1;
	int Nt = (int)ParameterOptimization::duration / dt;
	pht.setTimeSteps(dt, Nt);

	pht.setInitValues(ParameterOptimization::mIs);
	/*this->pht->finiteStep();
	std::cout << heatDeviation() << std::endl;;*/

	int j[GlobalVariables::numberOfOptiParam];
	//j[0] = ParameterFlags::k;
	j[0] = ParameterFlags::k_a;
	j[1] = ParameterFlags::k_b;
	pht.modifyInitValues(x_, j, GlobalVariables::numberOfOptiParam);


	pht.finiteStep();//ParameterOptimization::T, l

	Matrix3D<double> m_soll = Matrix3D<double>(60, 60, 1);
	Matrix3D<double> m_is = Matrix3D<double>(60, 60, 1);
	Matrix3D<double> m_slice = Matrix3D<double>(60, 60, 1);
	int index = 0;
	for (int x = 0; x < 60; x++)
	{
		for (int y = 0; y < 60; y++)
		{
			(m_soll)[{x, y, 0}] = ParameterOptimization::measValues->at(index);
			(m_slice)[{x, y, 0}] = (*(pht.getNewValues()))[{x, 30, y}].T;
			//std::cout << ParameterOptimization::positionInMatrix->at(index)[0] << " "<<ParameterOptimization::positionInMatrix->at(index)[1] << " " << ParameterOptimization::positionInMatrix->at(index)[2] << std::endl;
			(m_is)[{x, y, 0}] = (*(pht.getNewValues()))[ParameterOptimization::measPositionInMatrix->at(index)].T;//ParameterOptimization::positionInMatrix->at(index)
			index++;
		}
	}
	//ParameterOptimization::debugViewer(&m_slice);
	//ParameterOptimization::debugViewer(&m_soll);
	//ParameterOptimization::debugViewer(&m_is);
	double sumSquared = 0;
	for (int i = 0; i < ParameterOptimization::measPositionInMatrix->size(); i++)
	{
		double t_soll = ParameterOptimization::measValues->at(i);
		double t_new = (*(pht.getNewValues()))[ParameterOptimization::measPositionInMatrix->at(i)].T;
		sumSquared += (t_soll - t_new) * (t_soll - t_new);
		fi[i] = t_soll - t_new;


		//jaccobi:
		/*double  t_is = (*ParameterOptimization::mIs)[ParameterOptimization::positionInMatrix->at(i)].T;
		jac[i][0] = (t_new - t_is) - [0][i][0]*/



		//std::cout << fi[i] << " ";
	}
	std::cout << "SumSquared: " << sqrt(sumSquared / ParameterOptimization::measPositionInMatrix->size());
	std::cout << std::endl;

	printf("%s\n", x.tostring(1).c_str());
	

}
void fctOptParam(const real_1d_array& x, real_1d_array& fi, void* ptr)
{
	//std::cout << "fctOptParam" << std::endl;

	float x_[GlobalVariables::numberOfOptiParam];
	//std::cout << "derzeitiger wert: ";
	for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
	{
		x_[i] = x[i];
		//cout << x_[i] << " ";
	}
	//ParameterOptimization::debugViewer(ParameterOptimization::mTSoll);
	//ParameterOptimization::debugViewer(ParameterOptimization::mIs);
	unsigned int dim[3] = { ParameterOptimization::simPar.dim.x,ParameterOptimization::simPar.dim.y,ParameterOptimization::simPar.dim.z };
	//PennesHeatTransfer_Gpu pht = PennesHeatTransfer_Gpu(dim);

	
	
	//pht.setTimeSteps(dt, Nt);
	
	//pht.setInitValues(ParameterOptimization::mIs);
#if GPU
	ParameterOptimization::pht->setT_d(ParameterOptimization::pht->getTsaved_d());
#else
	ParameterOptimization::pht->setT_h(ParameterOptimization::pht->getTsaved_h());
#endif
	
	//(*ParameterOptimization::pht->par_h) = Parameters((ParameterOptimization::simPar));

	int Nt = GlobalVariables::timeBetweenMeas / ParameterOptimization::pht->par_h->dt;
	/*this->pht->finiteStep();
	std::cout << heatDeviation() << std::endl;;*/
	
	/*int j[GlobalVariables::numberOfOptiParam];*/
	//j[0] = ParameterFlags::k;
	/*for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
		j[i] = GlobalVariables::paramToOpti[i];*/
	
	ParameterOptimization::modifyParameter(ParameterOptimization::pht->par_h, x_);
	//pht.modifyInitValues(x_, j, GlobalVariables::numberOfOptiParam);//,GlobalVariables::slices_time_f[GlobalVariables::slice_index]//todo


	ParameterOptimization::pht->finiteStep(Nt);//ParameterOptimization::T, l

	Matrix3D<double> m_soll = Matrix3D<double>(60, 60, 1);
	Matrix3D<double> m_is = Matrix3D<double>(60, 60, 1);
	Matrix3D<double> m_slice= Matrix3D<double>(60, 60, 1);

	
	Matrix3D<float> *Tnew_h = ParameterOptimization::pht->getTnew_h();
	//ParameterOptimization::debugViewer(&m_slice);
	//ParameterOptimization::debugViewer(&m_soll);
	//Debugger::debugViewer(&m_is,"mIs");
	//std::cout << "size soll " << ParameterOptimization::measPositionInMatrix->size() << std::endl;
	double sumSquared = 0;
	//std::cout <<"dsaf " << ParameterOptimization::measPositionInMatrix->size() << std::endl;
	for (int i = 0; i < ParameterOptimization::measPositionInMatrix->size(); i++)
	{
			double t_soll = ParameterOptimization::measValues->at(i);
			double t_new = (*Tnew_h)[ParameterOptimization::measPositionInMatrix->at(i)];
			
			//(*(pht.getNewValues()))[ParameterOptimization::measPositionInMatrix->at(i)].T = ParameterOptimization::measValues[i];
			sumSquared += (t_soll - t_new) * (t_soll - t_new);//abs(t_soll - t_new);
			fi[i] = t_soll - t_new;		

			
			//std::cout << fi[i] << " ";
	}
	ParameterOptimization::squaredError[GlobalVariables::slice_index] = sqrt(sumSquared / ParameterOptimization::measPositionInMatrix->size());
	//std::cout << "---------------------SumSquared: "<<sqrt(sumSquared / ParameterOptimization::measPositionInMatrix->size())<<std::endl;
	
	/*int index = 0;
	for (int x = 0; x < 60; x++)
	{
		for (int y = 0; y < 60; y++)
		{
			(m_slice)[{x, y, 0}] = (*(pht.getNewValues()))[ParameterOptimization::slicePositionInMatrix->at(index)].T;
			index++;
		}
	}
	Debugger::debugViewer(&m_slice, "mslice");*/

	///printf("%s\n", x.tostring(1).c_str());
}
void ParameterOptimization::optimize()
{
	std::cout << "---------------optimize------------------" << std::endl;

	int nbOfMeasPoints = measPositionInMatrix->size();
//	nbOfMeasPoints = 3600;
	if (nbOfMeasPoints == 0)
	{
		std::cout << "nbOfMeasPoints+ " << nbOfMeasPoints << std::endl;
		return;
	}

	//welche Schicht
	slice_index = 0;
	/*std::cout << "x_is: ";
	for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
		std::cout << x_is[GlobalVariables::slice_index][i]<<" ";*/
	//std::cout << std::endl;
	//double heatsumIs = 0;
	//double heatsumSoll = 0;

	////x0.assign(0, 0.59);
	//for (int i = 0; i < dim[0]; i++)
	//{
	//	for (int j = 0; j < dim[1]; j++)
	//	{
	//		for (int k = 0; k < dim[2]; k++)
	//		{
	//			{
	//				heatsumIs += (*mIs)[{i, j, k}].T;
	//				heatsumSoll += (*_mTSoll)[{i, j, k}];
	//				//x0.assign(index, (*mIs)[{i, j, k}].Q_r);//
	//				//index++;
	//				//std::cout << x0[index] << std::endl;
	//			}


	//		}
	//	}
	//}

	/*std::cout << "initial Heat: " << heatsumIs / (dim[0] * dim[1] * dim[2]) << std::endl;
	std::cout << "Soll Heat: " << heatsumSoll / (dim[0] * dim[1] * dim[2]) << std::endl;*/


	double epsx = 0.00000001;// 0.0000000001     0.00001;
	ae_int_t maxits = 10;
	minlmstate state;
	minlmreport rep;
	//double bndl_ = x_is[GlobalVariables::slice_index][0] - x_is[GlobalVariables::slice_index][0] * 0.1 > 0 ? x_is[GlobalVariables::slice_index][0] - x_is[GlobalVariables::slice_index][0] * 0.1 : 0;
	//double bndu_ = x_is[GlobalVariables::slice_index][0] - x_is[GlobalVariables::slice_index][0] * 0.1;
	double bndl_[GlobalVariables::numberOfOptiParam];
	double bndu_[GlobalVariables::numberOfOptiParam];
	bndl_[0] = 80;
	bndl_[1] = 0.1;
	//bndl_[2] = 0;
//	bndl_[3] = 0.0;
	/*bndl_[2] = 3;
	bndl_[3] = 1;*/
	//bndl_[2] = 2;

	bndu_[0] = 300;
	bndu_[1] = 5;
//	bndu_[2] = 10000;
	//expD
	//bndu_[2] = 1;
	//bndu_[3] = 0.1;
	/*bndu_[2] = 8;
	bndu_[3] = 5;*/
	//bndu_[2] = 6;

	real_1d_array bndl;
	real_1d_array bndu;
	bndl.setcontent(GlobalVariables::numberOfOptiParam, bndl_);
	bndu.setcontent(GlobalVariables::numberOfOptiParam, bndu_);
	//
	// Create optimizer, tell it to:
	// * use numerical differentiation with step equal to 0.0001
	// * use unit scale for all variables (s is a unit vector)
	// * stop after short enough step (less than epsx)
	//
	//real_1d_array x_ = x_is[GlobalVariables::slice_index];
	
	real_1d_array x;
	real_1d_array s;

	x= x_is[GlobalVariables::slice_index];
	s = s_is[GlobalVariables::slice_index];
	//int Nx = 60;
	//double x_[Nx];
	//double s_[Nx];
	//for (int i = 0; i < 60; i++)
	//{
	//	x_[i] = 50;
	//	s_[i] = 1;
	//}
	//std::cout << std::endl;
	//x.setcontent(Nx, x_);
	//s.setcontent(Nx, s_);
	/*for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;*/
//	std::cout << "xxxx: " << x[0] <<x[1]<<x[2]<< std::endl;
	minlmcreatev(GlobalVariables::numberOfOptiParam, nbOfMeasPoints, x, 0.1, state);//0.0001
	minlmsetcond(state, epsx, maxits);
	minlmsetbc(state, bndl, bndu);
	minlmsetscale(state, s_is[GlobalVariables::slice_index]);

	//
	// Optimize
	//
	alglib::minlmoptimize(state, fctOptParam);

	//
	// Test optimization results
	//
	// NOTE: because we use numerical differentiation, we do not
	//       verify Jacobian correctness - it is always "correct".
	//       However, if you switch to analytic gradient, consider
	//       checking it with OptGuard (see other examples).
	//
	minlmresults(state, x, rep);
	//std::cout << "coreect x " << std::endl;
	//printf("%s\n", x.tostring(1).c_str()); // EXPECTED: [-3,+3]
	//std::cout <<"td_0 "<< x[0] * 0 + x[1]<<" td_21 " << x[0] * 21 + x[1] << " td_100 " << x[0] * 100 + x[1] << std::endl;
	//update x 
	
	for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
	{
		x_is[GlobalVariables::slice_index][i] = x[i];
	}


	/*std::cout << "alle Optiparam" << std::endl;
	for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
	{
		for (int j = 0; j < GlobalVariables::numberOfOptiParam; j++)
		{
			std::cout << x_is[i][j] << " ";
		}
		std::cout << squaredError[i];
		std::cout << std::endl;
		
	}*/
	
	//DEBUG
	//PennesHeatTransfer pht = PennesHeatTransfer(true);
	//pht.setInitValues(mIs);
	//std::cout << "Mis size " << mIs->getDimensions()[0] << " " << mIs->getDimensions()[1] << " "<<mIs->getDimensions()[2] << std::endl;
	//double dt = GlobalVariables::dt;
	//int Nt = (int)duration / dt;
	//pht.setTimeSteps(dt, Nt);

	//
	//int j[GlobalVariables::numberOfOptiParam];
	//for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
	//	j[i] = GlobalVariables::paramToOpti[i];

	//pht.modifyInitValues(&x_is[0][0], j, GlobalVariables::numberOfOptiParam);//, GlobalVariables::slices_time_f[GlobalVariables::GlobalVariables::slice_index]

	//pht.finiteStep();




	//for (int i = 0; i < ParameterOptimization::measPositionInMatrix->size(); i++)
	//{
	//	(*(pht.getNewValues()))[ParameterOptimization::measPositionInMatrix->at(i)].T = ParameterOptimization::measValues->at(i);
	////	std::cout << " " << pht.getNewValues()->getAngle(ParameterOptimization::measPositionInMatrix->at(i)) << " ";
	//	
	//}
	//std::cout << std::endl;


	//int index = 0;
	//Matrix3D<double> m_slice = Matrix3D<double>(60, 60, 1);
	//for (int x = 0; x < 60; x++)
	//{
	//	for (int y = 0; y < 60; y++)
	//	{
	//		//int y_ = dim[1] - y - 1;
	//		(m_slice)[{x, y, 0}] = (*(pht.getNewValues()))[(*ParameterOptimization::slicePositionInMatrix)[index]].T;
	//		//std::cout << slicePositionInMatrix->at(index)[0] << " " << slicePositionInMatrix->at(index)[1] << " " << slicePositionInMatrix->at(index)[2]<<std::endl;
	//		index++;
	//	}
	//}

	//Debugger::debugViewer(&m_slice,"mSlice");
	
}


void ParameterOptimization::modifyParameter(Parameters* _p,  float* optiParam)
{
	for (int p = 0; p < GlobalVariables::numberOfOptiParam; p++)
	{
		switch (GlobalVariables::paramToOpti[p])
		{
		case ParameterFlags::td_a:
		{
			//P_n->getAngle({ i, j, k });
			//diff = std::abs(angle - P_n->getAngle({ i, j, k }));
			/*if (angle != 1000 && diff <= 22.5)
			{
				(*P_n)[{i, j, k}].td_a = (1 - diff / 22.5) * (*P_n)[{i, j, k}].td_a + (diff / 22.5) * x[p];
			}
			else*/
			_p->td_a = optiParam[p];

			break;
		}
		case ParameterFlags::td_b:
			_p->td_b = optiParam[p];
			break;
		case ParameterFlags::td_c:
			_p->td_c = optiParam[p];
			break;
		case ParameterFlags::T_ind_max:
		{

			ParameterOptimization::par2opti.T_ind_max = optiParam[p];
			for (int i = 0; i < _p->N_heat; i++)
			{
				_p->T_heat[i] = GlobalVariables::baseLineT + ParameterOptimization::par2opti.Q_rel[i] * (ParameterOptimization::par2opti.T_ind_max - GlobalVariables::baseLineT);
				//std::cout << _p->T_heat[i]<< " ";
			}
			//std::cout << std::endl;
			break;
		}
		case ParameterFlags::c_b:
		{
			_p->cb = optiParam[p];
			break;
		}
		default:
			break;
		}
	}
}
void ParameterOptimization::adjustParameterNeedleAxis(const vector<double> _QalongNeedle, Matrix3D<PennesEquationParameter>* _mIs)
{
	//TODO
	for (int y = 0; y < _mIs->getDimensions()[1]; y++)
	{
		(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].Q_r_rel = _QalongNeedle[y];
		
		std::cout << _QalongNeedle[y] << " ";
		
	}
	std::cout << std::endl;
}

void ParameterOptimization::adjustParameterNeedleAxis(const vector<int> TalongNeedle, Matrix3D<PennesEquationParameter>* _mIs)
{
	//TODO
	float devHeat;
	for (int y = 0; y < _mIs->getDimensions()[1]; y++)
	{
		//devHeat = (*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T - (float)TalongNeedle[y];
		
		//std::cout << "oldT: " << (int)(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " ";
		
		//(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T -= (int)devHeat;
		
		
		(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T = TalongNeedle[y];
		(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].k_a = 30;
		(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].k_b = 0;
		/*for (int i = -1; i < 2; i++)
		{
			for (int j = -1; j < 1; j++)
			{
				(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0] + i, y, GlobalVariables::pointOnNeedle_volume[2]+j}].k_a = 30;
			}
			
		}*/
		
		
		//(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].k_b
		/*(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].Q_r += (int)(devHeat *
			((*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].rho*
				(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].c));*/
		std::cout << TalongNeedle[y] << " ";
		//std::cout <<"newT: " <<(int)(*_mIs)[{GlobalVariables::pointOnNeedle_volume[0], y, GlobalVariables::pointOnNeedle_volume[2]}].T << " devheat " << devHeat << " talong " << (float)TalongNeedle[y] << std::endl;
	}
	std::cout << std::endl;
}
