//----------------------------------------------------------------------------------
//! Class NecrosisMapComputation.
/*!
// \file
// \author  Julian Alpers
// \date    2020-12-10
//
// Basic class for computing the necrosis map from 2D phase images.
*/
//----------------------------------------------------------------------------------

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageShiftScale.h>
#include <vtkImageMathematics.h>
#include<vector>
#define _USE_MATH_DEFINES
#include <math.h>

class NecrosisMapComputation
{
public:
	NecrosisMapComputation();
	~NecrosisMapComputation();

	/*float _magneticFieldStrength;
	float _echoTime;*/
	void computeHeatMap(
		vtkSmartPointer<vtkImageData> p_vtk_phaseReferece, 
		vtkSmartPointer<vtkImageData> p_vtk_phase, 
		vtkSmartPointer<vtkImageData> p_vtk_heatMap, 
		float f_magneticFieldStrength,
		float f_echoTime);
	void computeReference(std::vector< vtkSmartPointer<vtkImageData>>& _phaseRefereceImages);
};