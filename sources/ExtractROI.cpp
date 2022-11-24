#include<ExtractROI.h>
#define mytpye double
ExtractROI::ExtractROI()
{
	
}

ExtractROI::~ExtractROI()
{
}
void ExtractROI::crop(vtkSmartPointer<vtkImageData>& _data, int _sliceNr)
{
	vtkNew<vtkExtractVOI> extractROI;
	
	extractROI->SetInputData(_data);
	extractROI->SetVOI(ROI[_sliceNr][2], ROI[_sliceNr][3], ROI[_sliceNr][0], ROI[_sliceNr][1], 0, 0);
	extractROI->SetInformation(_data->GetInformation());
	extractROI->UpdateInformation();
	extractROI->Update();
	
	//extractROI->PropagateUpdateExtent();

	_data = extractROI->GetOutput();
//	_data->Modified();
	int* extractedDims = _data->GetDimensions();
	std::cout << "Dims: "
		<< " x: " << extractedDims[0] << " y: " << extractedDims[1]
		<< " z: " << extractedDims[2] << std::endl;

}
void  ExtractROI::defineROI(vtkSmartPointer<vtkImageData> _referenceData, int _sliceNr)
{
	int* dim = _referenceData->GetDimensions();
	int midpoint[2];

	midpoint[0] = dim[0] / 2;
	midpoint[1] = dim[1] / 2;

	const int width =10 ;

	int window = 20;


	double * grad;
	double* grad_acc;
	double grad_sum;

	double* intensityChangeOfGrad;





	int x_l[width];
	int x_r[width];
	int y_l[width];
	int y_r[width];

	for (int x = -width/2; x < width / 2; x++)
	{
		int size = dim[1] - 1;
		grad = new double[size]();
		grad_acc = new double[size]();
		grad_sum = 0;

		intensityChangeOfGrad = new double[size]{ 0 };

		for (int y = 0; y < size; y++)
		{
			double a = _referenceData->GetScalarComponentAsDouble(midpoint[0] + x, y,0,0);
			double b = _referenceData->GetScalarComponentAsDouble(midpoint[0] + x, y +1, 0, 0);
			grad[y] = abs(a - b);
			//std::cout <<std::round( grad[y]) << " ";
			grad_sum = grad_sum + (grad[y]);
			grad_acc[y] = grad_sum;
		}
		//std::cout<<std::endl;
		//std::cout << grad_sum << std::endl;

		for (int i = window-1; i < size - window-1; i++)
		{
			double left = std::max(grad_acc[i], grad_acc[i + window - 1]) - std::min(grad_acc[i], grad_acc[i + window - 1]);
			double right = std::max(grad_acc[i], grad_acc[i - window + 1]) - std::min(grad_acc[i], grad_acc[i - window + 1]);

			if (left > right)
				intensityChangeOfGrad[i] = (left / right);
			else
				intensityChangeOfGrad[i] = (right / left);
			//std::cout << intensityChangeOfGrad[i] << " ";
		}
		//std::cout << std::endl;
	
		int pl = std::distance(intensityChangeOfGrad,std::max_element(intensityChangeOfGrad, intensityChangeOfGrad + midpoint[1]-1));
		int pr = std::distance(intensityChangeOfGrad, std::max_element(intensityChangeOfGrad + midpoint[1]-1, intensityChangeOfGrad + dim[1]-1));
		//std::cout <<"p " << pl <<" "<<pr <<std::endl;
		x_l[x + width / 2] = pl;
		x_r[x + width / 2] = pr;


		delete[] grad;
		delete[] grad_acc;
		delete[] intensityChangeOfGrad;
	}

	//y
	//std::cout << "y" << std::endl;
	for (int y = -width / 2; y < width / 2; y++)
	{
		int size = dim[0] - 1;
		//std::cout << "size" << size << std::endl;
		grad = new double[size]();
		grad_acc = new double[size]();
		grad_sum = 0;

		intensityChangeOfGrad = new double[size] { 0 };

		for (int x = 0; x < size; x++)
		{
			double a = _referenceData->GetScalarComponentAsDouble(x, midpoint[1] + y, 0, 0);
			double b = _referenceData->GetScalarComponentAsDouble(x +1, midpoint[1] + y, 0, 0);
			grad[x] = abs(a - b);
			//std::cout <<std::round( grad[y]) << " ";
			grad_sum = grad_sum + (grad[x]);
			grad_acc[x] = grad_sum;
		}
		//std::cout << grad_sum << std::endl;

		for (int i = window-1; i < size - window-1; i++)
		{
			//double left = std::max(grad_acc[i], grad_acc[i + window - 1]) - std::min(grad_acc[i], grad_acc[i + window - 1]);
			//double right = std::max(grad_acc[i], grad_acc[i - window + 1]) - std::min(grad_acc[i], grad_acc[i - window + 1]);
			double left = abs(grad_acc[i] - grad_acc[i - window + 1]);
			double right = abs(grad_acc[i] - grad_acc[i + window - 1]);
			if (left > right)
				intensityChangeOfGrad[i] = (left / right);
			else
				intensityChangeOfGrad[i] = (right / left);
			//std::cout << intensityChangeOfGrad[i] << " ";
		}
		//std::cout << std::endl;
			
		int pl = std::distance(intensityChangeOfGrad, std::max_element(intensityChangeOfGrad, intensityChangeOfGrad + midpoint[0]-1));
		int pr = std::distance(intensityChangeOfGrad, std::max_element(intensityChangeOfGrad + midpoint[0]-1, intensityChangeOfGrad + dim[0]-1));
		std::cout <<"p " << pl <<" "<<pr <<std::endl;
		y_l[y + width / 2] = pl;
		y_r[y + width / 2] = pr;


		delete[] grad;
		delete[] grad_acc;
		delete[] intensityChangeOfGrad;
	}
	

	
	//ROI[_sliceNr][0] = std::round(CalcMedian(x_l, width));
	//ROI[_sliceNr][1] = std::round(CalcMedian(x_r, width));
	ROI[_sliceNr][0] = std::round(std::accumulate(x_l, x_l + width, 0.0f)/width);
	ROI[_sliceNr][1] = std::round(std::accumulate(x_r, x_r + width, 0.0f)/width);

	ROI[_sliceNr][2] = std::round(std::accumulate(y_l, y_l + width, 0.0f) / width);
	ROI[_sliceNr][3] = std::round(std::accumulate(y_r, y_r + width, 0.0f) / width);

	std::cout << "ROI[_sliceNr][ " << _sliceNr<<"] "<< ROI[_sliceNr][0] << " " << ROI[_sliceNr][1] << " " << ROI[_sliceNr][2] << " " << ROI[_sliceNr][3] << std::endl;



}
float ExtractROI::CalcMedian(int* arr, int N)
{
	std::sort(arr, arr + N);
	if (N % 2 == 0)
	{
		return (arr[N / 2 - 1] + arr[N / 2]) / 2;
	}
	else
	{
		return arr[N / 2];
	}
	
}