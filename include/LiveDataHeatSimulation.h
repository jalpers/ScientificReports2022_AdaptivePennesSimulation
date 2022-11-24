// MinGW und Libcmaes
#ifndef LIVEDATAHEATSIMULATION_H
#define LIVEDATAHEATSIMULATION_H
/*
*	Include QT header.
*/
#include <QTimer>
#include <QDebug>
#include <QMainWindow>
#include <QFileDialog>
//#include <QVTKWidget.h>
#include <QMainWindow>
#include <QFileDialog>
/*
*	Include VTK header.
*/
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkDICOMImageReader.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>

#include <vtkImageShiftScale.h>
#include<vtkResliceImageViewer.h>
#include<vtkImageActor.h>


#include <vtkNamedColors.h>
#include <vtkContourValues.h>
#include<vtkImageViewer.h>
/*
*	Own includes.
*/
#include <DicomHandler.h>
#include<DataVolume.h>
#include<PennesHeatTransfer.h>
#include<ParameterOptimization.h>
#include<VesselMap.h>
#include "Test_cuda.h"
#include<PennesHeatTransfer_gpu.h>
#include<ThermalDose.h>
#include<ExtractROI.h>
#include<NecrosisMapComputation.h>
#include<Matrix3D.h>

#include "vtkBoundedPlanePointPlacer.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkDICOMImageReader.h"
#include "vtkDistanceRepresentation.h"
#include "vtkDistanceRepresentation2D.h"
#include "vtkDistanceWidget.h"
//#include <vtkGenericOpenGLRenderWindow.h>
#include "vtkHandleRepresentation.h"
#include "vtkImageData.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkImageSlabReslice.h"
#include "vtkInteractorStyleImage.h"
#include "vtkLookupTable.h"
#include "vtkPlane.h"
#include "vtkPlaneSource.h"
#include "vtkPointHandleRepresentation2D.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkProperty.h"
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkRenderWindowInteractor.h"
#include "vtkResliceImageViewer.h"
#include "vtkResliceCursorLineRepresentation.h"
#include "vtkResliceCursorThickLineRepresentation.h"
#include "vtkResliceCursorWidget.h"
#include "vtkResliceCursorActor.h"
#include "vtkResliceCursorPolyDataAlgorithm.h"
#include "vtkResliceCursor.h"
#include "vtkResliceImageViewerMeasurements.h"


#include "vtkSmartPointer.h"
#include "vtkResliceImageViewer.h"
#include "vtkImagePlaneWidget.h"
#include "vtkDistanceWidget.h"
#include "vtkResliceImageViewerMeasurements.h"
#include <QMainWindow>

//#include<QVTKInteractor.h>
#include<vtkImageActor.h>
#include<vtkImageMapToColors.h>
#include<vtkImageMapper3D.h>


QT_BEGIN_NAMESPACE
namespace Ui { class LiveDataHeatSimulation; }
QT_END_NAMESPACE

class LiveDataHeatSimulation : public QMainWindow
{
    Q_OBJECT

public:
	LiveDataHeatSimulation(QWidget *parent = nullptr);
    ~LiveDataHeatSimulation();


private slots:
	
	void updateImage();
	/*!
	*	\brief Load specified dicom file.
	*	\param q_c_dicomFilePath Path to the dicom file.
	* 	\param q_c_viewer Specification of the desired viewer to change its content.
	*	\return none
	*/

	void slotButtonBrowse();
	/*!
	*	\brief Slot function connected to the "Compute heat map" button.
	*	\param none
	*	\return none
	*/
	void slotButtonStartSimulation();


	
	

private:
	void startSimulation();
	void Test();
	void startSimulationUI(DataVolume* _volume);

	void init3DVolumeViewer(vtkSmartPointer<vtkImageData> _volume);
	void initHeatMapViewer(vtkSmartPointer<vtkImageData> _currentSlice);
	void initSlicerVolumeViewer(vtkSmartPointer<vtkImageData> _volume);

	void updateUI();
	vtkSmartPointer<vtkRenderWindow> qt_vtk_renderWindow_3Dvolume;
	vtkSmartPointer<vtkRenderWindow> qt_vtk_renderWindow_HeatMap;
	vtkSmartPointer<vtkRenderWindow> qt_vtk_renderWindow_SlicerVolume;

	vtkSmartPointer<vtkImageViewer2> imageViewer;
	

	bool isCropped;
    Ui::LiveDataHeatSimulation *ui;
    QFileDialog* m_q_dialog ;			                    //!< Pointer to the dialog to browse files and directories.
	QTimer* timer;

	QString q_c_dicomFilePath ;
	std::string vesselFilePath;
	//
	QStringList phantoms = { "PerfusionPhantom_1","PerfusionPhantom_2","PerfusionPhantom_3","PerfusionPhantom_4","PerfusionPhantom_5", "PerfusionPhantom_6","Phantom_1","Phantom_2","Phantom_3","Phantom_4","Phantom_5","TempPhantom_1","TempPhantom_2"
	};
//perfu5 fehlt
	QString root_direc;
	// 
	int nrOfTimeSteps_list[13] = { 14,14,14,14,14,13,14,13,13,13,13,13,13 };//14,14,14,14,14,13,14,13,13,13,13
	//QStringList phantoms = { "Phantom_2" };


	QString tubeFilePath;
	int fileID;
	int timestep;
	QStringList angleList = { "0","22_5","45","67_5","90","112_5","135","157_5" };//"112_5"
	//QStringList slices_time = { "0","90","45","135","22_5","112_5","67_5","157_5" };
	QStringList slices_time = { "0","22_5","45","67_5","90","112_5","135","157_5" };
	DataVolume volume;
	DataVolume UIVolume;
	
	int nrOfTimeSteps;
	int currentSliceIndex;
	int currentPhantomIndex;
	PennesHeatTransfer_Gpu* pht;
	ParameterOptimization* parOpt;
	ThermalDose* thermD;

	std::vector<std::array<double, GlobalVariables::numberOfOptiParam>> optiParamCache;
protected:
	std::vector<double> squaredErrorCache;
	std::vector<std::chrono::milliseconds> optiTimeCache;
	vtkSmartPointer< vtkResliceImageViewer > riw[2];
	vtkSmartPointer< vtkImagePlaneWidget > planeWidget[2];
	vtkSmartPointer< vtkDistanceWidget > DistanceWidget[2];
	vtkSmartPointer< vtkResliceImageViewerMeasurements > ResliceMeasurements;

	std::chrono::steady_clock::time_point begin_simul;
	std::chrono::steady_clock::time_point end_simul;
	//double counter()
	//{
	//	return 0.001 * GetTickCount();
	//}

};
#endif // LIVEDATAHEATSIMULATION_H