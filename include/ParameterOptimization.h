#ifndef PARAMETEROPTIMIZATION_H
#define PARAMETEROPTIMIZATION_H

#include <vector>
#include<Matrix3D.h>
#include <SimulationRessources.h>
#include<PennesHeatTransfer.h>
#include<PennesHeatTransfer_Gpu.h>

/*
*	Include VTK header.
*/
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include<qstringlist.h>
//#include <PennesHeatTransfer.h>
//#include<Debugger.h>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;

struct Parameter2Optimize
{
    
    float td_a;
    float td_b;
    float td_c;
    float T_ind_max;
    float Q_rel[128];
    Parameter2Optimize() {};
};
class ParameterOptimization{
public:
    ParameterOptimization();
    ParameterOptimization(PennesHeatTransfer_Gpu& _pth);
    ~ParameterOptimization();

    //make is const
    void optimize(Matrix3D<PennesEquationParameter>* _mIs, Matrix3D<double>* _mTSoll, HasToBeOptimized _htbOpti, double _time);

    void setParam(Parameters* _p) { simPar = Parameters(*_p); }
    void optimize();
    void setIsParam(Matrix3D<PennesEquationParameter>* _mIs) { 
        mIs = _mIs; }
    void setMeasurement(std::vector<double>* _measValues, std::vector<std::array<int, 3>>* _measPositionInMatrix, std::vector<std::array<int, 3>>* _slicePositionInMatrix)
    {
        measValues = _measValues; measPositionInMatrix = _measPositionInMatrix; slicePositionInMatrix = _slicePositionInMatrix;
    }
    float* getOptiParamter(int _slice_index)
    {
        float* res = new float[GlobalVariables::numberOfOptiParam];
        int minElementIndex = std::min_element(squaredError, squaredError + GlobalVariables::numberOfSlices) - squaredError;
        std::cout << minElementIndex << std::endl;
        for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
        {
            //res[i] = (x_is[_slice_index][i] + x_is[minElementIndex][i])/2;
            //ohne zeitliche Optimuerng:
            res[i] = x_is[_slice_index][i];
        }
        
        return res;
    }
    double getSquaredError(int _slice_index)
    {
        return squaredError[_slice_index];
    }
    void optimizeTest();
    void optimize_LM(Matrix3D<PennesEquationParameter>* _mIs, Matrix3D<double>* _mTSoll, HasToBeOptimized _htbOpti, double _time);

    double heatDeviation();
    void deleteStuff()
    {
        delete measPositionInMatrix;
        delete slicePositionInMatrix;
        delete measValues;
    }
    //static void debugViewer(Matrix3D<PennesEquationParameter>* _b);
    //static void debugViewer(Matrix3D<double>* _b);
    //static void debugViewer(const std::vector<std::vector<int>>& _b);
    static int T[60];

    static Matrix3D<double>* mTSoll;


    static Matrix3D<PennesEquationParameter>* mIs;
    static std::vector<double>* measValues;
    static std::vector<std::array<int, 3>>* measPositionInMatrix;
    static std::vector<std::array<int, 3>>* slicePositionInMatrix;
    static double duration;
    static void adjustParameterNeedleAxis(const vector<int> TalongNeedle, Matrix3D<PennesEquationParameter>* _mIs);
    static void adjustParameterNeedleAxis(const vector<double> _QalongNeedle, Matrix3D<PennesEquationParameter>* _mIs);
    static double sumSquaredDistance(Matrix3D<PennesEquationParameter>* _b);
    static double squaredError[GlobalVariables::numberOfSlices];
    //static PennesHeatTransfer_Gpu* pht;
    static Parameter2Optimize par2opti;
    static Parameters simPar;
    static PennesHeatTransfer_Gpu* pht;


    static void modifyParameter(Parameters* _P, float* x);
    

   
    
private:
    
    //static FitFunc optfct;
    //static FitFunc fsphere;
    int bestSlice;
    real_1d_array x_is[GlobalVariables::numberOfSlices];
    real_1d_array s_is[GlobalVariables::numberOfSlices];
    real_2d_array jac_is[GlobalVariables::numberOfSlices];
    
    int slice_index = 0;

    double* setUpParamArray();

    
    
    

};

#endif