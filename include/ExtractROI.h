#include<DicomHandler.h>
#include<vtkImageData.h>
#include<SimulationRessources.h>
#include <vtkExtractVOI.h>
class ExtractROI
{
public:
	ExtractROI();
	~ExtractROI();
	void crop(vtkSmartPointer<vtkImageData>& _data, int _sliceNr);
	void defineROI(vtkSmartPointer<vtkImageData> _referenceData, int _sliceNr);
private:
	int ROI[GlobalVariables::numberOfSlices][4];
	float CalcMedian(int* arr, int N);
};

