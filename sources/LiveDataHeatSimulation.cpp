#include "LiveDataHeatSimulation.h"
#include "ui_LiveDataHeatSimulation.h"
#include<vtkWin32OutputWindow.h>
//#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImagePlaneWidget.h>
#include<Pathfinder.h>
#include <thread>
#include <mutex>
#include<Debugger.h>

#define  _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING 1
#define VISUAL 0
//#include <filesystem>
#include<experimental/filesystem>
//#include<cuda_util.h>

//std::mutex m;

//https://discourse.vtk.org/t/help-i-need-compiled-qt5-15-vtk8-2-debug-means-the-compiled-version-is-vtk8-2-debug-and-support-qt5-15-development/4965
class vtkImageInteractionCallback : public vtkCommand
{
public:
	static vtkImageInteractionCallback* New() { return new vtkImageInteractionCallback; }

	vtkImageInteractionCallback()
	{
		this->Slicing = 0;
		this->sliceNr = 30;
		this->ImageReslice = nullptr;
		this->Interactor = nullptr;
	}

	void SetImageReslice(vtkImageReslice* reslice) { this->ImageReslice = reslice; }

	vtkImageReslice* GetImageReslice() { return this->ImageReslice; }

	void SetInteractor(vtkRenderWindowInteractor* interactor) { this->Interactor = interactor; }

	vtkRenderWindowInteractor* GetInteractor() { return this->Interactor; }

	void Execute(vtkObject*, unsigned long event, void*) override
	{
		vtkRenderWindowInteractor* interactor = this->GetInteractor();

		int lastPos[2];
		
		interactor->GetLastEventPosition(lastPos);

		int currPos[2];
		interactor->GetEventPosition(currPos);

		if (event == vtkCommand::LeftButtonPressEvent)
		{
			this->Slicing = 1;
		}
		else if (event == vtkCommand::LeftButtonReleaseEvent)
		{
			this->Slicing = 0;
		}
		else if(event == vtkCommand::MouseWheelBackwardEvent)
		{
			sliceNr++;
			std::cout << "backward" << std::endl;
			vtkImageReslice* reslice = this->ImageReslice;
			double sliceSpacing = reslice->GetOutput()->GetSpacing()[2];
			vtkMatrix4x4* matrix = reslice->GetResliceAxes();

			double point[4];
			double center[4];
			point[0] = 0.0;
			point[1] = 0.0;
			point[2] = sliceSpacing * sliceNr;
			point[3] = 1.0;
			matrix->MultiplyPoint(point, center);
			matrix->SetElement(0, 3, center[0]);
			matrix->SetElement(1, 3, center[1]);
			matrix->SetElement(2, 3, center[2]);
			interactor->Render();
		}
		else if (event == vtkCommand::MouseWheelForwardEvent)
		{
			sliceNr--;
			std::cout << "backward" << std::endl;
			vtkImageReslice* reslice = this->ImageReslice;
			double sliceSpacing = reslice->GetOutput()->GetSpacing()[2];
			vtkMatrix4x4* matrix = reslice->GetResliceAxes();

			double point[4];
			double center[4];
			point[0] = 0.0;
			point[1] = 0.0;
			point[2] = sliceSpacing * sliceNr;
			point[3] = 1.0;
			matrix->MultiplyPoint(point, center);
			matrix->SetElement(0, 3, center[0]);
			matrix->SetElement(1, 3, center[1]);
			matrix->SetElement(2, 3, center[2]);
			interactor->Render();
		}
		else if (event == vtkCommand::MouseMoveEvent)
		{
			if (this->Slicing)
			{
				vtkImageReslice* reslice = this->ImageReslice;

				// Increment slice position by deltaY of mouse
				int deltaY = lastPos[1] - currPos[1];

				reslice->Update();
				double sliceSpacing = reslice->GetOutput()->GetSpacing()[2];
				vtkMatrix4x4* matrix = reslice->GetResliceAxes();
				// move the center point that we are slicing through
				double point[4];
				double center[4];
				point[0] = 0.0;
				point[1] = 0.0;
				point[2] = sliceSpacing * deltaY;
				point[3] = 1.0;
				matrix->MultiplyPoint(point, center);
				matrix->SetElement(0, 3, center[0]);
				matrix->SetElement(1, 3, center[1]);
				matrix->SetElement(2, 3, center[2]);
				interactor->Render();
			}
			else
			{
				vtkInteractorStyle* style =
					vtkInteractorStyle::SafeDownCast(interactor->GetInteractorStyle());
				if (style)
				{
					style->OnMouseMove();
				}
			}
		}
	}

private:
	// Actions (slicing only, for now)
	int Slicing;
	int sliceNr;
	// Pointer to vtkImageReslice
	vtkImageReslice* ImageReslice;

	// Pointer to the interactor
	vtkRenderWindowInteractor* Interactor;
};


class vtkResliceCursorCallback : public vtkCommand
{
public:
	static vtkResliceCursorCallback* New()
	{
		return new vtkResliceCursorCallback;
	}

	void Execute(vtkObject* caller, unsigned long ev,
		void* callData) override
	{

		if (ev == vtkResliceCursorWidget::WindowLevelEvent ||
			ev == vtkCommand::WindowLevelEvent ||
			ev == vtkResliceCursorWidget::ResliceThicknessChangedEvent)
		{
			// Render everything
			for (int i = 0; i < 2; i++)
			{
				this->RCW[i]->Render();
			}
			this->IPW[0]->GetInteractor()->GetRenderWindow()->Render();
			return;
		}

		vtkImagePlaneWidget* ipw =
			dynamic_cast<vtkImagePlaneWidget*>(caller);
		if (ipw)
		{
			double* wl = static_cast<double*>(callData);

			if (ipw == this->IPW[0])
			{
				this->IPW[1]->SetWindowLevel(wl[0], wl[1], 1);
			}
			else if (ipw == this->IPW[1])
			{
				this->IPW[0]->SetWindowLevel(wl[0], wl[1], 1);
			}

		}

		vtkResliceCursorWidget* rcw = dynamic_cast<
			vtkResliceCursorWidget*>(caller);
		if (rcw)
		{
			vtkResliceCursorLineRepresentation* rep = dynamic_cast<
				vtkResliceCursorLineRepresentation*>(rcw->GetRepresentation());
			// Although the return value is not used, we keep the get calls
			// in case they had side-effects
			rep->GetResliceCursorActor()->GetCursorAlgorithm()->GetResliceCursor();
			for (int i = 0; i < 2; i++)
			{
				vtkPlaneSource* ps = static_cast<vtkPlaneSource*>(
					this->IPW[i]->GetPolyDataAlgorithm());
				ps->SetOrigin(this->RCW[i]->GetResliceCursorRepresentation()->
					GetPlaneSource()->GetOrigin());
				ps->SetPoint1(this->RCW[i]->GetResliceCursorRepresentation()->
					GetPlaneSource()->GetPoint1());
				ps->SetPoint2(this->RCW[i]->GetResliceCursorRepresentation()->
					GetPlaneSource()->GetPoint2());

				// If the reslice plane has modified, update it on the 3D widget
				this->IPW[i]->UpdatePlacement();
			}
		}

		// Render everything
		for (int i = 0; i < 2; i++)
		{
			this->RCW[i]->Render();
		}
		this->IPW[0]->GetInteractor()->GetRenderWindow()->Render();
	}

	vtkResliceCursorCallback() {}
	vtkImagePlaneWidget* IPW[2];
	vtkResliceCursorWidget* RCW[2];
};



LiveDataHeatSimulation::LiveDataHeatSimulation(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::LiveDataHeatSimulation)
{
	
	ui->setupUi(this);
	
	qDebug() << "hello Ui";
	
	/*
	*	Set default value for file path.
	*/
	//https://github.com/martijnkoopman/Qt-VTK-viewer/issues/3
	ui->qt_lineEdit_filename->setText("C:/Users/maxro/Documents/MasterInformatik/Einzelprojekt/Phantome/Phantom_3");
	//ui->qt_lineEdit_filename->setText("C:/Users/Maximilian Roetzer/Downloads/2.5D_Thermometry_Study/Phantom_2/HeatMap");
	/*
	*	Initialize members.
	*/

	/*
	*	Connect signals.
	*/
	
	//root_direc = "C:/2.5D_Thermometry_Study/";
	root_direc = "C:/Users/Maximilian Roetzer/Downloads/2022-03-26_Leber_10min/Ablation/";//2022-03-26_Leber_8min

	timer = new QTimer(this);
	fileID = 0;
	
	isCropped = true;
	connect(timer, SIGNAL(timeout()), this, SLOT(updateImage()));
	QObject::connect(ui->qt_pushButton_browse, SIGNAL(clicked()), this, SLOT(slotButtonBrowse()));
	QObject::connect(ui->qt_button_StartSimulation, SIGNAL(clicked()), this, SLOT(slotButtonStartSimulation()));

	Test();

	std::cout << std::endl<< "endstartLiveDataSimulation" << std::endl;

}

LiveDataHeatSimulation::~LiveDataHeatSimulation()
{
	delete ui;
	delete pht;
}

/************************************************************************************
*								    slotButtonBrowse
*************************************************************************************/
void LiveDataHeatSimulation::slotButtonBrowse()
{
	/*
	*	Open file browser to specify file and set the text field.
	*/

	QString q_c_filename = m_q_dialog->getExistingDirectory();
	//QString q_c_filename = m_q_dialog->getOpenFileName();

	ui->qt_lineEdit_filename->setText(q_c_filename);
}
/************************************************************************************
*								    slotButtonStartSimulation
*************************************************************************************/
void LiveDataHeatSimulation::slotButtonStartSimulation()
{
	std::cout << "buttonStartSimulation" << std::endl;
	currentPhantomIndex = 0;
	//ui->qt_lineEdit_filename->setText(root_direc + phantoms[currentPhantomIndex]);
	ui->qt_lineEdit_filename->setText(root_direc);
	currentPhantomIndex = 0;

	startSimulation();
	

}
/************************************************************************************
*								    startSimulation
*************************************************************************************/
void LiveDataHeatSimulation::startSimulation()
{
	std::cout << "StartSimulation" << std::endl;
	int timeBetweenData = 1;

	q_c_dicomFilePath = ui->qt_lineEdit_filename->text() + "/HeatMap";
	vesselFilePath = ui->qt_lineEdit_filename->text().toStdString() + "/GroundTruth/Tubes.dcm";

	std::string la = ui->qt_lineEdit_filename->text().toStdString() + "/HeatMap/157_5/";
	/*DIR* dp;
	int i = 0;
	struct dirent* ep;
	dp = opendir(la.c_str());

	if (dp != NULL)
	{
		while (ep = readdir(dp))
			i++;

		(void)closedir(dp);
	}
	else
		perror("Couldn't open the directory");
	nrOfTimeSteps = i - 2;
	printf("There's %d files in the current directory.\n", nrOfTimeSteps);*/
	/*auto dirIter = std::experimental::filesystem::directory_iterator(la);
	int fileCount = 0;*/

	/*for (auto& entry : dirIter)
	{
			++fileCount;
	}
	std::cout << "filecount " << fileCount<<std::endl;
	nrOfTimeSteps = fileCount-1;*/
	nrOfTimeSteps = 22;//nrOfTimeSteps_list[currentPhantomIndex];

	//
	//Definieren des Volumens (Ausdehnung, Auflösung, Position, Orientierung, Nadelachse) durch Berücksichtigung der Metadaten der Baselineaufnahmen
	//
	volume = DataVolume(q_c_dicomFilePath);
	UIVolume = DataVolume(q_c_dicomFilePath);
	


	fileID = 0;
	timestep = 0;

#if VISUAL
	startSimulationUI(&UIVolume);
#endif

	unsigned int dim[3] = { volume.getImageData()->GetDimensions()[0], volume.getImageData()->GetDimensions()[1] ,volume.getImageData()->GetDimensions()[2] };
	//std::cout << "dim " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
	pht = new PennesHeatTransfer_Gpu(dim);
	//Anfangsparameter der Simulation
	pht->setDefaultParameters();
	double* spacing = volume.getImageData()->GetSpacing();
	
	for (int i = 0; i < 3; i++)
	{
		pht->par_h->dx_2[i] = spacing[i] * spacing[i] * 0.001 * 0.001;/* 0.001 * 0.001*/
	}
	 
	pht->par_h->dt = GlobalVariables::dt;

	//Baseline Temperatur als Temperaturinitialwert 
	pht->setAll2(GlobalVariables::baseLineT);
	pht->saveTimestep();
	
	int Nt = GlobalVariables::timeBetweenMeas / pht->par_h->dt;
	std::cout << "dt " << pht->par_h->dt << std::endl;
	
	
	//
	// Load VesselMap
	//
	VesselMap vesselMap = VesselMap(volume.getWorld2VoxelConverter(), volume.getImageData()->GetDimensions());
	if (q_c_dicomFilePath.contains(new QString("Perfusion")))
	{
		qDebug() << "isPerfusion";
		vesselMap.load(vesselFilePath);
	}
	pht->setPerfusion(*(vesselMap.getMap()));
	

	//
	//Baseline Karten werden eingeladen
	//
	for (int i = 0; i < slices_time.size(); i++)
	{
		volume.addSlice(slices_time[i], QString::fromStdString((std::to_string(timestep))));
	}


	//Initialer Wert für die Nadeltemperatur
	
	volume.defaultTempAlongNeedle();
	parOpt->par2opti.T_ind_max = 120;
	for (int i = 0; i < volume.posNeedleInVolume()->size(); i++)
	{
		if(i > volume.posNeedleInVolume()->size()/4 || i < volume.posNeedleInVolume()->size() - volume.posNeedleInVolume()->size() / 4)//weil Ränder werden nie erreicht
			pht->par_h->T_heat[i] = 100;
		pht->par_h->heat_x[i] = volume.posNeedleInVolume()->at(i)[0];
		pht->par_h->heat_y[i] = volume.posNeedleInVolume()->at(i)[1];
		pht->par_h->heat_z[i] = volume.posNeedleInVolume()->at(i)[2];
		/*std::cout << pht->par_h->T_heat[i] << " ";*/
	}
	pht->par_h->N_heat = volume.getQAlongNeedle().size();
	//std::cout<< pht->par_h->N_heat <<std::endl;
	

	timestep = 1;

	pht->saveTimestep();//TODO

	
	parOpt = new ParameterOptimization(*pht);//
	thermD = new ThermalDose(dim[0] * dim[1] * dim[2]);
	timer->start(timeBetweenData);//1000
	begin_simul = std::chrono::steady_clock::now();
	std::cout << "end Setup" << std::endl;
}

void LiveDataHeatSimulation::startSimulationUI(DataVolume* _volume)
{
	std::cout << "startSimulationUI" << std::endl;
	//setup UI


	//initSlicerVolumeViewer(_volume->getImageData());
	//TODO vielleicht noch ein Copy machen, ui und data unterschiedlich
	//_volume->getImageData()->DeepCopy();

	init3DVolumeViewer(_volume->getImageData());
	initSlicerVolumeViewer(_volume->getImageData());
	initHeatMapViewer(volume.getcurrentSlice());
	//for (int i = 0; i < 2; i++)
	//{
	//	riw[i]->SetResliceMode( 1 );
	//	riw[i]->GetRenderer()->ResetCamera();
	//	riw[i]->Render();
	//}
	//for (int i = 0; i < 2; i++)
	//{
	//	riw[i]->Render();
	//}
	//this->ui->qt_viewer_slicer->GetRenderWindow()->Render();
}


void LiveDataHeatSimulation::updateImage()
{
	std::cout << std::endl;
	qDebug() << "New Slice -> Update Volume";

	/*vtkSmartPointer<vtkRenderWindow> q_vtk_renderWindow_Right = vtkSmartPointer<vtkRenderWindow>::New();
	q_vtk_renderWindow_Right = ui->q_vtk_viewerRight->GetRenderWindow();

	vtkSmartPointer<vtkRenderWindow> q_vtk_renderWindow_Left = vtkSmartPointer<vtkRenderWindow>::New();
	q_vtk_renderWindow_Left = ui->q_vtk_viewerLeft->GetRenderWindow();*/

	//Trigger Visualiserung
#if VISUAL
	updateUI();
#endif
	

	for (int i = 0; i < pht->par_h->N_heat; i++)
	{
		parOpt->par2opti.Q_rel[i] = volume.getQAlongNeedle().at(i);
		pht->par_h->T_heat[i] = GlobalVariables::baseLineT + parOpt->par2opti.Q_rel[i] * (parOpt->par2opti.T_ind_max - GlobalVariables::baseLineT);
		
	}


	pht->saveTimestep();//TODO
	
	

	//Hier verstreicht die Zeit zwischen den Optimierugen
	//hier müsste man immer wieder Daten  von GPU fetchen um kontinuierliche Visualierung zu erhalten
	// also eher for-schleife mit pht->finiteStep(1);
	int Nt = GlobalVariables::timeBetweenMeas / GlobalVariables::dt;
	pht->finiteStep(Nt);

	//Update UIVolume
#if DEBUG
	Matrix3D<float>* Tnew_H = pht->getTnew_h();
	int* dim = UIVolume.getImageData()->GetDimensions();
	for (int i = 0; i < dim[0]; i++)
	{
		for (int j = 0; j < dim[1]; j++)
		{
			for (int k = 0; k < dim[2]; k++)
			{
				unsigned char* v = static_cast<unsigned char*>(UIVolume.getImageData()->GetScalarPointer(i, j, k));
				unsigned char l = (unsigned char)(int)std::round((*Tnew_H)[{i, j, k}]);
				v[0] = l;
			}
		}

	}
	UIVolume.getImageData()->Modified();

#endif



	//CEM43 modell
	thermD->update(pht->getTnew_h()->getData(), 0.1);


	GlobalVariables::slice_index = fileID;
	//Neue Schicht
	// Extrahieren der Isotherme
	//Bestimmen von Position und relative Stärke der Wärmequellen
	volume.addSlice(slices_time[fileID], QString::fromStdString((std::to_string(timestep))));


	//OPtimierung der Koeffizietnen
	if (true)
	{
		


		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		parOpt->setMeasurement(volume.currentMeasValues(), volume.getcurrentMeasPos(), volume.getcurrentSlicePos());
		parOpt->setParam(pht->par_h);
		parOpt->optimize();
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		float* optiParam = parOpt->getOptiParamter(GlobalVariables::slice_index);
		double squaredError = parOpt->getSquaredError(GlobalVariables::slice_index);

		ParameterOptimization::modifyParameter(pht->par_h, optiParam);

		//das ist nur für Debug
		std::array<double, GlobalVariables::numberOfOptiParam> par_arr;
		for (int i = 0; i < GlobalVariables::numberOfOptiParam; i++)
		{
			std::cout << " ------------------- " << optiParam[i];
			par_arr[i] = optiParam[i];
		}
		std::cout << "--Error-- " << squaredError;
		std::cout << std::endl;
		optiParamCache.push_back(par_arr);
		squaredErrorCache.push_back(squaredError);
		delete optiParam;

		optiTimeCache.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin));
		std::cout << "Time difference - OPTI = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
		//parOpt->deleteStuff();
	}
	

//im nächsten Schritt nächste Wärmekarte
	fileID++;
	qDebug() << "id: " << fileID << "step: " << timestep;
	if (fileID >= slices_time.size())
	{
		//nächster Zeitschritt
		timestep++;
		fileID = 0;


		if (timestep == 2 || timestep == 9)
		{
			std::cout << "------------------------------------------------------------------------------" << std::endl;
			qDebug() << q_c_dicomFilePath;
			int* dim = UIVolume.getImageData()->GetDimensions();
			//Wärmekarte
			Matrix3D<float>* Tnew_h = pht->getTnew_h();
			for (int i = 0; i < dim[0]; i++)
			{
				for (int j = 0; j < dim[1]; j++)
				{
					for (int k = 0; k < dim[2]; k++)
					{


						unsigned char* v = static_cast<unsigned char*>(UIVolume.getImageData()->GetScalarPointer(i, j, k));
						unsigned char l = (unsigned char)(int)std::round((*Tnew_h)[{i, j, k}]);
						v[0] = l;
					}
				}
			}
			UIVolume.getImageData()->Modified();
			if(timestep == 2)
				UIVolume.writeToFile(q_c_dicomFilePath + "/output1");
			if(timestep == 9)
				UIVolume.writeToFile(q_c_dicomFilePath + "/output8");
		}
	
	}
	if (timestep > nrOfTimeSteps-1)//14
	{
		//letzer Zeitschritt
		timer->stop();
#if VISUAL
		updateUI();
#endif
		qDebug() << q_c_dicomFilePath;
		int* dim = UIVolume.getImageData()->GetDimensions();
		//Wärmekarte
		Matrix3D<float>* Tnew_h = pht->getTnew_h();
		for (int i = 0; i < dim[0]; i++)
		{
			for (int j = 0; j < dim[1]; j++)
			{
				for (int k = 0; k < dim[2]; k++)
				{


					unsigned char* v = static_cast<unsigned char*>(UIVolume.getImageData()->GetScalarPointer(i, j, k));
					unsigned char l = (unsigned char)(int)std::round((*Tnew_h)[{i, j, k}]);
					v[0] = l;
				}
			}
		}
		UIVolume.getImageData()->Modified();
		UIVolume.writeToFile(q_c_dicomFilePath);

		//Nekrosekarte
		//thermD->computeNecrosisZones();
		//Matrix3D<float>* nec = new Matrix3D<float>(dim[0], dim[1], dim[2], thermD->getNecrosiszones());//
		////UIVolume.getImageData()->AllocateScalars(VTK_FLOAT, 1);//
		////UIVolume.getImageData()->Modified();
		////UIVolume.getImageData()->SetScalarType(VTK_FLOAT);
		//for (int i = 0; i < dim[0]; i++)
		//{
		//	for (int j = 0; j < dim[1]; j++)
		//	{
		//		for (int k = 0; k < dim[2]; k++)
		//		{


		//			unsigned char* v = static_cast<unsigned char*>(UIVolume.getImageData()->GetScalarPointer(i, j, k));//
		//			//std::cout << (*nec)[{i, j, k}] << " ";
		//			unsigned char l = (unsigned char)(int)std::round((*nec)[{i, j, k}]);//
		//			v[0] = l;

		//		//	UIVolume.getImageData()->SetScalarComponentFromFloat(i, j, k, 0, (*nec)[{i, j, k}]);
		//			

		//			/*if(i==GlobalVariables::pointOnNeedle_volume[0] && k == GlobalVariables::pointOnNeedle_volume[2])
		//				std::cout <<" nadle: " <<(*P_new)[{i, j, k}].T<<" "<< (*P_n)[{i, j, k}].T ;*/
		//				//sumHeat += (*P_new)[{i, j, k}].T;


		//		}
		//		/*if (i == (int)dim[1] / 2)
		//			std::cout << std::endl;*/
		//	}

		//}
		//UIVolume.getImageData()->Modified();
		//UIVolume.writeToFile(q_c_dicomFilePath);
		

		end_simul = std::chrono::steady_clock::now();
		std::cout << "Time difference - SIMUL = " << std::chrono::duration_cast<std::chrono::seconds> (end_simul - begin_simul).count() << "[s]" << std::endl;
		
		for (int i = 0; i < optiParamCache.size(); i++)
		{
			std::cout << optiParamCache[i][0] << " " << optiParamCache[i][1] << " error: " <<squaredErrorCache[i]<< std::endl;
		}
		for (int i = 0; i < optiTimeCache.size(); i++)
		{
			std::cout << optiTimeCache[i].count() << std::endl;
		}
		if (currentPhantomIndex < phantoms.size())
		{
			currentPhantomIndex++;
			ui->qt_lineEdit_filename->setText(root_direc + phantoms[currentPhantomIndex]);
			startSimulation();
		}

	}
}

void LiveDataHeatSimulation::updateUI()
{
	qDebug() << "updateUI";
	


	qt_vtk_renderWindow_3Dvolume->Render();
	qt_vtk_renderWindow_SlicerVolume->Render();


	qt_vtk_renderWindow_HeatMap->Render();
	initHeatMapViewer(volume.getcurrentSlice());


	qt_vtk_renderWindow_HeatMap->Render();

	initSlicerVolumeViewer(volume.getImageData());




}








void LiveDataHeatSimulation::init3DVolumeViewer(vtkSmartPointer<vtkImageData> _volume)
{
	int low = 30;// 10
	int middle =50;
	int high = 200;//220

	qt_vtk_renderWindow_3Dvolume = vtkSmartPointer<vtkRenderWindow>::New();
	qt_vtk_renderWindow_3Dvolume = ui->qt_viewer_3DVolume->GetRenderWindow();

	/*
	*	vtkNamedColors bietet eine einfache Farbdatenbank mit entsprechenden Namen, die man sich einfach ausgeben lassen kann.
	*	Beim Setzen der Hintergrundfarbe habe ich zum Beispiel die Farbe "Wheat" genommen.
	*/
	vtkSmartPointer<vtkNamedColors> namedColors = vtkSmartPointer<vtkNamedColors>::New();
	/*
	*	Hier wird das entsprechende RenderWindow und der dazu gehrige Renderer geholt.
	*/
	vtkSmartPointer<vtkRenderer> q_vtk_renderer = vtkSmartPointer<vtkRenderer>::New();
	qt_vtk_renderWindow_3Dvolume->AddRenderer(q_vtk_renderer);
	/*
	*	Hier gibts den Interactor und den 3D InteractorStyle.
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> q_vtk_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> q_vtk_interactorStyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	q_vtk_interactor->SetInteractorStyle(q_vtk_interactorStyle);
	q_vtk_interactor->SetRenderWindow(qt_vtk_renderWindow_3Dvolume);
	q_vtk_interactor->Modified();
	/*
	*	Als Mapper hab ich mich in der Flle von Mglichkeiten fr den OpenGL GPU Volume Cast Mapper entschieden, der in der Lage ist
	*	IsoSurfaces vernnftig darzustellen weil die entsprechenden Shader bereits integriert sind.
	*/
	//vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	//vtkFixedPointVolumeRayCastMapper
	volumeMapper->SetInputData(_volume);

	/*
	*	Die automatische Sample Distanz knnte je nach Auflsung des Datensatzes zu Problemen fhren. Ich hab sie jetzt auf 0.5 gesetzt aber
	*	da knnt ihr einfach mal ein bisschen rumspielen mit den Einstellungen. Wann sieht es gut aus und ist trotzdem performant?
	*	Anschlieend wird der BlendMode auf IsoSurface gesetzt. Da knnt ihr aber auch nochmal rumspielen, vielleicht findet ihr was cooleres.
	*/
	volumeMapper->AutoAdjustSampleDistancesOff();
	volumeMapper->SetSampleDistance(1.0);//0.5
	volumeMapper->SetBlendModeToIsoSurface();
	
	volumeMapper->Modified();
	/*
	*	Jetzt wird es ekelhaft. Um vernnftig zu visualisieren msst ihr eine richtige Transferfunktion definieren fr die Isowerte. Das geschieht mit
	*	zwei sogenannten LookupTables. In der vtkColorTransferFunction definiert ihr, welcher Voxelwert welche Farbe bekommt. In diesem Fall ist alles
	*	mit dem Wert 10 grn und alles mit dem Wert 220 rot.
	*/
	vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
	color->RemoveAllPoints();
	for (int i = 0; i < low; i++)
		color->AddRGBPoint(i, 0.0, 1.0, 0.0);

	for (int j = low; j < high; j++)
		color->AddRGBPoint(j, 1.0, 0.0, 0.0);
	color->AddRGBPoint(middle, 0.0, 0.0, 1.0);
	color->Modified();
	//color->AddRGBPoint(1, 0.0, 1.0, 0.0);
	//color->AddRGBPoint(0, 1.0, 0.0, 0.0);
	/*
	*	Mit der vtkPiecewiseFunction definiert ihr die Opacitt eines entsprechenden Wertes. Die Punkte, die ihr hier eintragt mssen identisch sein mit den
	*	Werten in der vtkColorTransferFunction. Hier haben Voxel mit dem Wert 10 eine Opazitt von 30% und Voxel mit einem Wert von 220 eine Opazitt von 60%.
	*/
	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->RemoveAllPoints();

	for (int i = 0; i < low; i++)
		compositeOpacity->AddPoint(i, 0.0);
	
	for (int j = low; j < high; j++)
		compositeOpacity->AddPoint(j, 1.0);
	compositeOpacity->AddPoint(middle, 1.0);
	compositeOpacity->Modified();
	//compositeOpacity->AddPoint(1, 0.3);
	//compositeOpacity->AddPoint(220, 0.6);
	/*
	*	Jetzt mssen wir die VolumeProperty definieren, dazu stellen wir das Shading an und eine einfache lineare interpolation. Bei der Interpolation gibt
	*	es bestimmt auch noch andere Mglichkeiten. SetColor bekommt die vtkColorTransferFunction und setScalarOpacity bekommt die
	*	vtkPiecewiseFunction.
	*/
	vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->ShadeOn();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetColor(color);
	volumeProperty->SetScalarOpacity(compositeOpacity);
	volumeProperty->Modified();
	/*
	*	Dann mssen wir das Volumen definieren mit dem entsprechenden Mapper und den vorab eingestellenten Properties.
	*/
	vtkSmartPointer<vtkVolume> vtk_volume = vtkSmartPointer<vtkVolume>::New();
	vtk_volume->SetMapper(volumeMapper);
	vtk_volume->SetProperty(volumeProperty);
	vtk_volume->Modified();
	/*
	*	Jetzt schmeien wir alles in den Renderer und setzen am Ende das DepthPeeling auf True, damit der Renderer wei
	*	in welcher Reihenfolge er die Strukturen visualisieren muss.
	*/
	q_vtk_renderer->AddVolume(vtk_volume);
	q_vtk_renderer->SetBackground(namedColors->GetColor3d("Wheat").GetData());//
	//q_vtk_renderer->GetActiveCamera()->UpdateViewport(q_vtk_renderer);
	//q_vtk_renderer->SetViewPoint(0, 1, 0);

	if (false)
	{
		q_vtk_renderer->GetActiveCamera()->SetFocalPoint(vtk_volume->GetCenter());
		q_vtk_renderer->GetActiveCamera()->SetPosition(0, 200, 0);
	}


	q_vtk_renderer->SetUseDepthPeeling(true);
	/*
	*	Ganz am Ende weisen wir unseren entsprechenden Isowerten auch noch die Isovalues zu fr die interne Verarbeitung.
	*	Hier reprsentiert unsere IsoSurface 0 den Wert 10 und unsere IsoSurface 1 den Wert 220. Auch diese Werte mssen
	*	identisch sein mit den Werten in der Transfer Funktion.
	*/
	volumeProperty->GetIsoSurfaceValues()->SetValue(0, low);
	volumeProperty->GetIsoSurfaceValues()->SetValue(2, high);

	qt_vtk_renderWindow_3Dvolume->Modified();
}

void LiveDataHeatSimulation::initHeatMapViewer(vtkSmartPointer<vtkImageData> _currentSlice)
{
	//qt_vtk_renderWindow_HeatMap = ui->qt_viewer_heatMap->GetRenderWindow();
	//vtkNew<vtkImageViewer> image_view;
	//// use our render window with image_view
	//image_view->SetRenderWindow(qt_vtk_renderWindow_HeatMap);
	//image_view->SetInputData(_currentSlice);
	//image_view->SetupInteractor(qt_vtk_renderWindow_HeatMap->GetInteractor());
	//image_view->SetColorLevel(138.5);
	//image_view->SetColorWindow(233);
	//image_view->Modified();

	//image_view->Render();

	//qt_vtk_renderWindow_HeatMap->Render();

	qt_vtk_renderWindow_HeatMap = ui->qt_viewer_heatMap->GetRenderWindow();
	/*
	*	Create container for image data properties.
	*/
	int extent[6];
	double origin[3];
	double spacing[3];
	/*
	*	Get the current interactor and set the interactor style to handle 2D images.
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> interactor;

	{
		interactor = qt_vtk_renderWindow_HeatMap->GetInteractor();
		vtkSmartPointer<vtkInteractorStyleImage> interactorStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
		interactor->SetInteractorStyle(interactorStyle);
	}


	/*
	*	Create 2D image viewer, connect it to the reader output port and specify the render window.
	*/
	imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	//imageViewer->SetupInteractor(interactor);
	 //  interactor->SetInteractorStyle(imageViewer->GetInteractorStyle());
	//imageViewer->SetInputData(reader->GetOutput());
	imageViewer->SetInputData(volume.getcurrentSlice());
	imageViewer->Modified();

	
	imageViewer->SetRenderWindow(qt_vtk_renderWindow_HeatMap);
	
	imageViewer->Render();
	/*
	*	Get the current active camera from the scene and turn on parallel projection.
	*/
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	camera = imageViewer->GetRenderer()->GetActiveCamera();
	camera->ParallelProjectionOn();
	/*
	*	Get the loaded image data properties.
	*/
	_currentSlice->GetExtent(extent);
	_currentSlice->GetOrigin(origin);
	_currentSlice->GetSpacing(spacing);
	/*
	*	Compute the image center in world coordinates and set the focal point of the camera.
	*/
	float f_centerX = origin[0] + 0.5 * (extent[0] + extent[1]) * spacing[0];
	float f_centerY = origin[1] + 0.5 * (extent[2] + extent[3]) * spacing[1];
	camera->SetFocalPoint(f_centerX, f_centerY, 0.0);
	/*
	*	Compute the center distance in world coordinates and set the parallel scale.
	*/
	float f_centerDistance = (extent[3] - extent[2] + 1) * spacing[1];
	camera->SetParallelScale(0.5f * f_centerDistance);
	/*
	*	Get the current camera distance and set camera position and distance.
	*/
	float f_cameraDistance = camera->GetDistance();
	camera->SetPosition(f_centerX, f_centerY, +f_cameraDistance);
	/*
	*	Rerender the render window to apply all the changes.
	*/

	{
		qt_vtk_renderWindow_HeatMap->Render();
	}


	//-----------------------------------------------------------------------------------------------------------


	//vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
	//actor->GetMapper()->SetInputData(volume.getImageData());
	////actor->GetMapper()->SetInputData(volume.getImageData());//my

	//vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	//renderer->AddActor(actor);

	//qt_vtk_renderWindow_HeatMap = ui->qt_viewer_slicer->renderWindow();
	//qt_vtk_renderWindow_HeatMap->AddRenderer(renderer);

	//// Set up the interaction
	//vtkSmartPointer<vtkInteractorStyleImage> imageStyle =
	//	vtkSmartPointer<vtkInteractorStyleImage>::New();
	//imageStyle->Modified();

	//vtkSmartPointer<vtkRenderWindowInteractor> interactor =
	//	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//interactor->SetInteractorStyle(imageStyle);
	//qt_vtk_renderWindow_HeatMap->SetInteractor(interactor);
	//qt_vtk_renderWindow_HeatMap->Render();


	//----------------------------------------------------------------------------------------------------------------------

	/*vtkSmartPointer<vtkRenderWindow> renderWindow = ui->qt_viewer_slicer->GetRenderWindow();
	renderWindow->SetSize(900, 900);
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(0, 1, 0);
	renderer->ResetCamera();
	renderWindow->AddRenderer(renderer);


	


	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = ui->qt_viewer_slicer->GetInteractor();
	vtkSmartPointer<vtkResliceImageViewer> resliceImageViewer = vtkSmartPointer<vtkResliceImageViewer>::New();
	resliceImageViewer->SetupInteractor(renderWindowInteractor);
	resliceImageViewer->SetRenderWindow(renderWindow);
	resliceImageViewer->SetSliceOrientationToXY();
	resliceImageViewer->SetRenderer(renderer);
	resliceImageViewer->SetInputData(volume.getImageData());
	resliceImageViewer->SetSliceOrientation(2);
	resliceImageViewer->SetResliceModeToAxisAligned();
	resliceImageViewer->SetColorLevel(40.0);
	resliceImageViewer->SetColorWindow(400.0);
	resliceImageViewer->SetResliceMode(0);
	resliceImageViewer->GetRenderer()->ResetCamera();
	resliceImageViewer->Render();



	resliceImageViewer->GetRenderWindow()->Render();
	renderWindowInteractor->Initialize();
	resliceImageViewer->GetRenderWindow()->Render();
	renderWindowInteractor->Start();*/
	


}
//https://stackoverflow.com/questions/31802089/vtkrenderwindowinteractor-doesnt-start-and-cause-program-freezes-java
//https://stackoverflow.com/questions/25318210/renderwindowinteractor-in-qvtkwidget
//http://ustk-doc.inria.fr/doxygen/ustk-daily/usMedicalImageViewer_8cpp_source.html

//https://programmer.group/modified-vtk-sample-code-of-vtk-dicom.html
//https://stackoverflow.com/questions/55120374/how-to-get-the-position-of-2d-dicom-image-slice-in-3d-surface-rendered-output-in

//https://programmerwiki.com/article/2716165429/
void LiveDataHeatSimulation::initSlicerVolumeViewer(vtkSmartPointer<vtkImageData> _volume)
{
	//https://stackoverflow.com/questions/38225343/vtkresliceimageviewer-display-incorrect-too-dark-dicom-image
//https://gitlab.kitware.com/vtk/vtk/-/blob/master/Examples/GUI/Qt/FourPaneViewer/QtVTKRenderWindows.cxx

	int extent[6];
	double spacing[3];
	double origin[3];

	_volume->GetExtent(extent);
	_volume->GetOrigin(origin);
	_volume->GetSpacing(spacing);

	double center[3];
	center[0] = origin[0] + spacing[0] * 0.5 * (extent[0] + extent[1]);
	center[1] = origin[1] + spacing[1] * 0.5 * (extent[2] + extent[3]);
	center[2] = origin[2] + spacing[2] * 0.5 * (extent[4] + extent[5]);
	//std::cout << "cent: "<<center[0] << center[1] << center[2] << std::endl;
	// Matrices for axial, coronal, sagittal, oblique view orientations
	static double axialElements[16] = {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 };

	static double coronalelements[16] = {
			1, 0, 0, 0,
			0, 0, 1, 0,
			0,-1, 0, 0,
			0, 0, 0, 1 };

	static double sagittalElements[16] = {
	0, 0, -1, 0, //
	1, 0, 0, 0,  //
	0, -1, 0, 0, //
	0, 0, 0, 1   //
	};

	static double obliqueElements[16] = {
			1, 0, 0, 0,
			0, 0.866025, -0.5, 0,
			0, 0.5, 0.866025, 0,
			0, 0, 0, 1 };


	

	// Set the slice orientation
	vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
	resliceAxes->DeepCopy(axialElements);//axial
	// Set the point through which to slice
	resliceAxes->SetElement(0, 3, center[0]);
	resliceAxes->SetElement(1, 3, center[1]);
	resliceAxes->SetElement(2, 3, center[2]);

	// Extract a slice in the desired orientation
	vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
	reslice->SetInputData(_volume);
	reslice->SetOutputDimensionality(2);
	reslice->SetResliceAxes(resliceAxes);
	reslice->SetInterpolationModeToLinear();
	

	// Create a greyscale lookup table
	vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
	table->SetRange(0, 200.0);            // image intensity range
	table->SetValueRange(0.0, 1.0);      // from black to white
	table->SetSaturationRange(0.0, 0.0); // no color saturation
	table->SetRampToLinear();
	table->Build();

	// Map the image through the lookup table
	vtkSmartPointer<vtkImageMapToColors> color = vtkSmartPointer<vtkImageMapToColors>::New();
	color->SetLookupTable(table);
	color->SetInputConnection(reslice->GetOutputPort());

	// Display the image
	vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
	actor->GetMapper()->SetInputConnection(color->GetOutputPort());
	//actor->GetMapper()->SetInputData(volume.getImageData());//my

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);

	////cameera
	//vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	//camera = renderer->GetActiveCamera();
	//camera->ParallelProjectionOn();
	///*
	//*	Compute the image center in world coordinates and set the focal point of the camera.
	//*/
	//float f_centerX = origin[0] + 0.5 * (extent[0] + extent[1]) * spacing[0];
	//float f_centerY = origin[1] + 0.5 * (extent[2] + extent[3]) * spacing[1];
	//float f_centerZ = origin[2] + 0.5 * (extent[4] + extent[5]) * spacing[1];
	//camera->SetFocalPoint(f_centerX, f_centerY, f_centerZ);
	///*
	//*	Compute the center distance in world coordinates and set the parallel scale.
	//*/
	//float f_centerDistance = (extent[3] - extent[2] + 1) * spacing[1];
	//camera->SetParallelScale(0.5f * f_centerDistance);
	///*
	//*	Get the current camera distance and set camera position and distance.
	//*/
	//float f_cameraDistance = camera->GetDistance();
	//camera->SetPosition(f_centerX, f_centerY, +f_cameraDistance);




	qt_vtk_renderWindow_SlicerVolume = ui->qt_viewer_slicer->GetRenderWindow();
	qt_vtk_renderWindow_SlicerVolume->AddRenderer(renderer);

	// Set up the interaction
	vtkSmartPointer<vtkInteractorStyleImage> imageStyle =
	vtkSmartPointer<vtkInteractorStyleImage>::New();
	imageStyle->SetInteractionModeToImage3D();
	imageStyle->Modified();

	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetInteractorStyle(imageStyle);
	qt_vtk_renderWindow_SlicerVolume->SetInteractor(interactor);
	qt_vtk_renderWindow_SlicerVolume->Render();
	
	vtkSmartPointer<vtkImageInteractionCallback> callback =
		vtkSmartPointer<vtkImageInteractionCallback>::New();
	callback->SetImageReslice(reslice);
	callback->SetInteractor(interactor);

	imageStyle->AddObserver(vtkCommand::MouseMoveEvent, callback);
	imageStyle->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);
	imageStyle->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
	imageStyle->AddObserver(vtkCommand::LeftButtonPressEvent, callback);
	imageStyle->AddObserver(vtkCommand::LeftButtonReleaseEvent, callback);

	// Start interaction
	// The Start() method doesn't return until the window is closed by the user
	//interactor->Start();



}



void LiveDataHeatSimulation::Test()
{
	//Test cropping

	DicomHandler* dicomHandler = new DicomHandler();
	ExtractROI* cpMap = new ExtractROI();
	NecrosisMapComputation* necComp = new NecrosisMapComputation();
	//QString fileDataName = "C:/Users/Maximilian Roetzer/Downloads/2022-03-26_Leber_8min/Ablation/Phase sortiert/1/IM00009.dcm";
	int numberOfSlices = GlobalVariables::numberOfSlices;
	//{ "135", "157_5", "0", "22_5", "45", "67_5", "90", "112_5"}
	QString fileName = "C:/Users/Maximilian Roetzer/Downloads/2022-03-26_Leber_10min/Baseline/Phase sortiert/alle/";
	for (int o = 1; o <= 8; o++)
	{
		std::vector< vtkSmartPointer<vtkImageData>> _phaseRefereceImages;
		for (int i = 0; i <= 2; i++)
		{
			int nr = o + numberOfSlices * i;
			QString fileDataName;
			if ((nr) < 10)
				fileDataName = fileName + "IM0000" + QString::number(nr) + ".dcm";
			else
				fileDataName = fileName + "IM000" + QString::number(nr) + ".dcm";
			//qDebug() << fileDataName;
			_phaseRefereceImages.push_back(dicomHandler->loadDicom(fileDataName));


			int* dims = _phaseRefereceImages[i]->GetDimensions();
			unsigned int sum = 0;
			for (int z = 0; z < dims[2]; z++)
			{
				for (int y = 0; y < dims[1]; y++)
				{
					for (int x = 0; x < dims[0]; x++)
					{
						unsigned short* pixel = static_cast<unsigned short*>(_phaseRefereceImages[i]->GetScalarPointer(x, y, z));
						sum += unsigned int(pixel[0]);
					}
				}
			}
			//std::cout << "sum " << sum/(dims[0]*dims[1]*dims[2]) << std::endl;

		}


	//std::cout << "geht" << std::endl;
	//QString fileDataName = "C:/Users/maxro/Documents/testCrop.dcm";
	
	
	//QString fileDataName = fileName + "/" + _angle  +".dcm";//.dcm _14.IMA
	
	
	necComp->computeReference(_phaseRefereceImages);

	//std::cout << "gehtRef" << std::endl;
	/*int* dims = _phaseRefereceImages[0]->GetDimensions();
	std::cout << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
	unsigned int sum = 0;
	for (int z = 0; z < dims[2]; z++)
	{
		for (int y = 0; y < dims[1]; y++)
		{
			for (int x = 0; x < dims[0]; x++)
			{
				unsigned short* pixel = static_cast<unsigned short*>(_phaseRefereceImages[0]->GetScalarPointer(x, y, z));
				sum += unsigned int(pixel[0]);
			}
		}
	}
	std::cout << "sum " << sum / (dims[0] * dims[1] * dims[2]) << std::endl;*/



	//vtkSmartPointer<vtkImageData> data = dicomHandler->loadDicom(fileDataName);



	cpMap->defineROI(_phaseRefereceImages[0], o-1);
	
	//int* extractedDims = _phaseRefereceImages[0]->GetDimensions();
	//std::cout << "Dims: "
	//	<< " x: " << extractedDims[0] << " y: " << extractedDims[1]
	//	<< " z: " << extractedDims[2] << std::endl;
	//cv::String str = "gallo";


	//--------------------------------------------------------------------------
	// 
	//vtkSmartPointer<vtkDICOMReader> reader = vtkSmartPointer<vtkDICOMReader>::New();
	//QString filename0 = "C:/Users/Maximilian Roetzer/Downloads/2022-03-26_Leber_8min/Ablation/Phase sortiert/1/IM00161.dcm";
	//reader->SetFileName(filename0.toStdString().c_str());
	//reader->Update();

	//DicomHandler dHandler;
	//DicomHandler::dicomDataProperties properties0;
	//dHandler.getImageDataProperties(reader, &properties0);
	//std::cout << "Halloooo" << std::endl;
	//vtkSmartPointer<vtkImageData> h = vtkSmartPointer<vtkImageData>::New();
	//necComp->computeHeatMap(_phaseRefereceImages[0], dicomHandler->loadDicom(filename0), _phaseRefereceImages[1], properties0.magneticFieldStrength, properties0.echoTime);
	////cpMap->crop(_phaseRefereceImages[1], 0);
	////int* dims = _phaseRefereceImages[0]->GetDimensions();
	////Matrix3D<double>* _b = new Matrix3D<double>(dims[0], dims[1], dims[2]);
	////for (int z = 0; z < dims[2]; z++)
	////{
	////	for (int y = 0; y < dims[1]; y++)
	////	{
	////		for (int x = 0; x < dims[0]; x++)
	////		{
	////			double b = _phaseRefereceImages[0]->GetScalarComponentAsDouble(x, y, z, 0);
	////			(*_b)[{x, y, z}] = b;
	////		}
	////	}
	////}

	////Debugger::debugViewer(_b,str);
	//vtkSmartPointer<vtkDICOMWriter> writer_d =
	//	vtkSmartPointer<vtkDICOMWriter>::New();
	//// Create a generator for MR images.
	//vtkNew<vtkDICOMMRGenerator> generator;

	//vtkSmartPointer<vtkMatrix4x4> p = vtkSmartPointer<vtkMatrix4x4>::New();
	//p->GlobalWarningDisplayOff();
	//p->Identity();
	//p->Modified();

	//generator->Modified();
	//writer_d->SetGenerator(generator.GetPointer());

	//// Set the output filename format as a printf-style string.
	//writer_d->SetFilePattern("%s/IM-0001-%04.3d.dcm");

	//std::string l = "C:/Users/Maximilian Roetzer/Documents/";
	//writer_d->SetFilePrefix(l.c_str());

	//writer_d->SetInputData(_phaseRefereceImages[1]);

	//writer_d->Modified();
	//writer_d->Write();

	//qDebug() << "write Dicom -Done";

	//--------------------------------------------------------------------------

	
	}

	delete cpMap;;
	delete dicomHandler;
	delete necComp;







	//Test GPU
	
	//lookup table für vessel position
	//

	//Test_Cuda cu = Test_Cuda();

	//cu.doCuda();
	//
	//std::cout << "--------------------------------------------------------------" << std::endl;

	//unsigned int dim[3] = { 60,60,60 };
	//int n = dim[0]*dim[1]*dim[2];
	////float* T_h = new float[n];
	////float* Tnew_h = new float[n];

	//double spacing[3] = { 0.005,0.005,0.005 };
	//Matrix3D<float>* T_h = new Matrix3D<float>(dim[0], dim[1], dim[2]);
	//Matrix3D<float>* Tnew_h = new Matrix3D<float>(dim[0], dim[1], dim[2]);
	//Matrix3D<PennesEquationParameter>* m_is = new Matrix3D<PennesEquationParameter>(dim[0], dim[1], dim[2]);
	//
	//m_is->setSpacing(spacing);
	//long float sum = 0;
	////for (int i = 0; i < n; i++)
	////{
	////	T_h[i] = 21;// i * 0.01;
	////	sum += T_h[i];
	////	Tnew_h[i] = 5;
	////	//m_is[i];
	////}

	//for (int i = 0; i < dim[0]; i++)
	//{
	//	for (int j = 0; j < dim[1]; j++)
	//	{
	//		for (int k = 0; k < dim[2]; k++)
	//		{
	//			(*m_is)[{i, j, k}].T =21;
	//			
	//			(*T_h)[{i, j, k}] = 21;
	//			sum += (*T_h)[{i, j, k}];
	//		}
	//	}
	//}

	//PennesHeatTransfer_Gpu pht_gpu = PennesHeatTransfer_Gpu(dim);
	//pht_gpu.setAll2(21);
	////pht_gpu.setT_h(T_h);
	//pht_gpu.setDefaultParameters();
	//pht_gpu.par_h->dx_2[0] = { 1.4 * 1e-6};
	//pht_gpu.par_h->dx_2[1] = { 1.4 * 1e-6 };
	//pht_gpu.par_h->dx_2[2] = { 1.4 * 1e-6 };
	//pht_gpu.par_h->N_heat = 60;
	//pht_gpu.par_h->dt = 1.0;
	//std::cout << "bound " << pht_gpu.par_h->boundCond << std::endl;
	//for (int i = 0; i < 60; i++)
	//{
	//	pht_gpu.par_h->heat_x[i] = 30;
	//	pht_gpu.par_h->heat_y[i] = i;
	//	pht_gpu.par_h->heat_z[i] = 30;
	//	pht_gpu.par_h->T_heat[i] = 25;

	//	(*m_is)[{30, i, 30}].Q_r_rel = 1;
	//	(*m_is)[{30, i, 30}].T_ind_max = 100;

	//}
	//PennesEquationParameter pn;
	//m_is->setAmbientValue(pn);
	//double spa[3] = { 1.4 * 1e-6 , 1.4 * 1e-6 , 1.4 * 1e-6 };
	//m_is->setSpacing(spa);
	//printf("sumSoll: %f \n", sum);
	//
	//
	//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	////for (int i = 0; i < 100; i++)
	////pht_gpu.finiteStep(100);
	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl; 
	////pht_gpu.getTnew_h();
	//pht_gpu.updateTnew_h();

	//
	//
	////PennesHeatTransfer pht = PennesHeatTransfer(false);
	////pht.setInitValues(m_is);
	////pht.setTimeSteps(0.1, 1);
	/////*std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();*/
	//////for(int i = 0; i<100;i++)
	////pht.finiteStep(0);
	/////*std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	////std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << "[ms]" << std::endl;
	////*/
	////sum = 0;
	////for (int i = 0; i < dim[0]; i++)
	////{
	////	for (int j = 0; j < dim[1]; j++)
	////	{
	////		for (int k = 0; k < dim[2]; k++)
	////		{
	////			sum += (*pht.getNewValues())[{i, j, k}].T;
	////		}
	////	}
	////}
	////printf(" \n sumIsRef: %f %f \n", sum, sum/(dim[0]* dim[1]* dim[2]));


	//delete m_is;





	////Test alglib
	//alglib::real_2d_array a, b, c;
	//int n = 2000;
	//int i, j;
	//double timeneeded, flops;

	//// Initialize arrays
	//a.setlength(n, n);
	//b.setlength(n, n);
	//c.setlength(n, n);
	//for (i = 0; i < n; i++)
	//	for (j = 0; j < n; j++)
	//	{
	//		a[i][j] = alglib::randomreal() - 0.5;
	//		b[i][j] = alglib::randomreal() - 0.5;
	//		c[i][j] = 0.0;
	//	}

	//// Set global threading settings (applied to all ALGLIB functions);
	//// default is to perform serial computations, unless parallel execution
	//// is activated. Parallel execution tries to utilize all cores; this
	//// behavior can be changed with alglib::setnworkers() call.
	//alglib::setglobalthreading(alglib::parallel);

	//// Perform matrix-matrix product.
	//flops = 2 * pow((double)n, (double)3);
	//timeneeded = counter();
	//alglib::rmatrixgemm(
	//	n, n, n,
	//	1.0,
	//	a, 0, 0, 0,
	//	b, 0, 0, 1,
	//	0.0,
	//	c, 0, 0);
	//timeneeded = counter() - timeneeded;

	//// Evaluate performance
	//printf("Performance is %.1f GFLOPS\n", (double)(1.0E-9 * flops / timeneeded));


	//TEST Cylinder Coordinates

	int a = 60;
	//Matrix3D<int> matrix = Matrix3D<int>(a, 1, a);
	//

	//for (int i = 0; i < a; i++)
	//{
	//	for (int j = 0; j < a; j++)
	//	{
	//		std::cout << matrix.getAngle({i,0,j}) << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout<<std::endl;
	//for (int i = 0; i < a; i++)
	//{
	//	std::cout << matrix.getAngle({ i,0, GlobalVariables::pointOnNeedle_volume[2] }) << " ";
	//}
	////Test Islinien
	//DicomHandler dH_i;
	//vtkSmartPointer<vtkImageData> img_iso = dH_i.loadDicom("C:/Users/maxro/Documents/MasterInformatik/Einzelprojekt/Phantome/Phantom_2/HeatMap/90/10.IMA");
	//std::cout << img_iso->GetScalarRange()[0] << " "<< DataVolume::encode(img_iso->GetScalarRange()[1]) << std::endl;
	//const size_t ROW= 60;
	//const size_t COL = 60;
	//vector<vector<int>> grid;
	//grid.resize(ROW);
	//for (int r = 0; r < ROW; r++)
	//	grid[r].resize(COL);
	////array<array<int, 10>, 9> grid{
	////{ { { 1, 0, 1, 1, 1, 1, 0, 1, 1, 1 } },
	////  { { 1, 1, 1, 0, 1, 1, 1, 0, 1, 1 } },
	////  { { 1, 1, 1, 0, 1, 1, 0, 1, 0, 1 } },
	////  { { 0, 0, 1, 0, 1, 0, 0, 0, 0, 1 } },
	////  { { 1, 1, 1, 0, 1, 1, 1, 0, 1, 0 } },
	////  { { 1, 0, 1, 0, 1, 1, 0, 1, 0, 0 } },
	////  { { 1, 0, 0, 0, 0, 1, 0, 0, 0, 1 } },
	////  { { 1, 0, 1, 1, 1, 1, 0, 1, 1, 1 } },
	////  { { 1, 1, 1, 0, 0, 0, 1, 0, 0, 1 } } }
	////};
	//Pathfinder::Pair start = { 0,30 };
	//Pathfinder::Pair end = { 59,30 };
	//for (size_t i = 0; i < 60; i++)
	//{
	//	for (size_t j = 0; j < 60; j++)
	//	{

	//		double zahl = (double)DataVolume::encode(static_cast<unsigned char*>(img_iso->GetScalarPointer(i, j, 0))[0]);

	//		if (zahl == 30)
	//			grid[i][j] = 1;
	//		else
	//			grid[i][j] = 4;
	//		std::cout << grid[i][j] << " ";
	//	}
	//	std::cout<<std::endl;
	//}
	//std::cout << std::endl;

	///*Pathfinder::Pair start = { 8,0 };
	//Pathfinder::Pair end = { 0,0 };*/
	//array <int, 3> t = { 1,2,3 };
	//Pathfinder pf;
	//pf.help(grid,start,end);

	////pf.help1();
	//vector<std::pair<int, int>> res;
	////pf.aStarSearch(grid, start, end,res);








	//Test explicit scheme;
//	DicomHandler dH_t;
//	vtkSmartPointer<vtkImageData> a_ = dH_t.loadDicom("C:/Users/maxro/Documents/MasterInformatik/Einzelprojekt/Phantome/Phantom_2/HeatMap/0/0.IMA");
//	PennesHeatTransfer pht_exp = PennesHeatTransfer(true);
//	int gl = 1;
//	Matrix3D<PennesEquationParameter>* initM = new Matrix3D<PennesEquationParameter>(60, 60, 60);
//	
//	//for (int i = 0; i < 60; i++)
//	//{
//	//	for (int j = 0; j < 60; j++)
//	//	{
//	//		for (int k = 0; k < gl; k++)
//	//		{
//	//			(*initM)[{i, j, k}].T = 100;
//	//		}
//	//	}
//	//}
//	(*initM)[{30, 30, 30}].T = 27;
////	(*initM)[{30, 30, 30}].keepTconstant = true;
//	initM->setSpacing(a_->GetSpacing());
//
//	pht_exp.setTimeSteps(0.01,600);
//	pht_exp.setInitValues(initM);
//	
////pht_exp.finiteStep(0);
//
//	ParameterOptimization::debugViewer(pht_exp.getNewValues());
//
//	pht_exp.setInitValues(initM);
//	pht_exp.setTimeSteps(3, 2);
//	pht_exp.finiteStep(0);
//	
//
//	ParameterOptimization::debugViewer(pht_exp.getNewValues());
//
//	delete initM;
//
//	ParameterOptimization* parOpt = new ParameterOptimization();
//	DicomHandler dH;
//	vtkSmartPointer<vtkImageData> a = dH.loadDicom("C:/Users/maxro/Documents/MasterInformatik/Einzelprojekt/Phantome/Phantom_2/HeatMap/0/0.IMA");
//	vtkSmartPointer<vtkImageData> b = dH.loadDicom("C:/Users/maxro/Documents/MasterInformatik/Einzelprojekt/Phantome/Phantom_2/HeatMap/0/1.IMA");
//
//
//	Matrix3D<double>* mTSoll = new Matrix3D<double>(60, 60, 1);
//	Matrix3D<PennesEquationParameter>* mIs = new Matrix3D<PennesEquationParameter>(60, 60, 1);
//
//	double heatsumA = 0;
//	double heatsumB = 0;
//	for (size_t i = 0; i < 60; i++)
//	{
//		for (size_t j = 0; j < 60; j++)
//		{
//
//			(*mIs)[{i, j, 0}].T = (double)DataVolume::encode(static_cast<unsigned char*>(a->GetScalarPointer(i, j, 0))[0]);
//			(*mTSoll)[{i, j, 0}] = (double)DataVolume::encode(static_cast<char*>(b->GetScalarPointer(i, j, 0))[0]);
//
//			heatsumA += (*mIs)[{i, j, 0}].T;
//			heatsumB += (*mTSoll)[{i, j, 0}];
//		}
//	}
//	std::cout << "heatsumA: " << heatsumA / 3600 << " hEatsumB: " << heatsumB / 3600 << std::endl;
//	mIs->setSpacing(a->GetSpacing());
//	mTSoll->setSpacing(b->GetSpacing());
//
//
//	PennesHeatTransfer pht_test = PennesHeatTransfer(true);
//	pht_test.setInitValues(mIs);
//	/*this->pht->finiteStep();
//	std::cout << heatDeviation() << std::endl;;*/
//	double dt = 0.1;
//	int time = 60;
//	int Nt = (int)time/ dt;
//	pht_test.setTimeSteps(dt, Nt);
//	double k= 10;
//	double Q_r[60];
//	std::fill_n(Q_r, 60, 50000000);
//	for (int i = 0; i < 60; i++)
//	{
//		std::cout<<Q_r[i];
//	}
//	std::cout << std::endl;
//	int Np = 1;
//	int j[Np];
//	pht_test.modifyInitValues(&k,j,Np);
//
//	int* dim = mIs->getDimensions();
//	/*dim[1] = 3;
//	dim[0] = 60;*/
//	int T[60];
//	for (int i = 0; i < 60;i++)
//	{
//		//if (i >= 30 && i <= 30)
//			//T[i] = 100;
//		T[i] = 50;//std::max({ (*mTSoll)[{i, 32, 0}], (*mTSoll)[{i, 31, 0}], (*mTSoll)[{i, 33, 0}] });
//	//	else
//		//	T[i] = 21;
//		//std::cout << T[i] << " ";
//	}
//	std::cout << std::endl;
//	int pN = 1;
//	int N = dim[0] * dim[1] * dim[2] * pN;
//	N = 1;
//	int* l = new int[2];
//	l[0] = 0;
//	l[1] = 32;
//	pht_test.finiteStep_opti();
//
//
//	//ParameterOptimization::debugViewer(pht_test.getNewValues());
//
//
//	for (int i = 0; i < dim[0]; i++)
//	{
//		for (int j = 0; j < dim[1]; j++)
//		{
//			(*mTSoll)[{i, j, 0}] = (*pht_test.getNewValues())[{i, j, 0}].T;
//		}
//	}
//
//
//	HasToBeOptimized _htbOpti;
//	//parOpt->optimizeTest();
//	////parOpt->optimize_LM(mIs,mTSoll, _htbOpti, 6);
//	parOpt->optimize_LM(mIs,mTSoll , _htbOpti, time);
//	delete mTSoll;
//	delete mIs;
//	delete parOpt;
}















