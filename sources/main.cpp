#include "LiveDataHeatSimulation.h"
#include <QSurfaceFormat>
#include "QVTKWidget.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication q_application(argc, argv);
	LiveDataHeatSimulation q_application_window;
//	vtkObject::SetGlobalWarningDisplay(0);
	/*
	*	Open and set stylesheet.
	*/
	//QFile File("..\\style\\stylesheet.qss");
	QFile File("C:/2.5DThermometryReconstruction/style/stylesheet.qss");
	File.open(QFile::ReadOnly);
	QString StyleSheet = QLatin1String(File.readAll());
	File.close();
	qApp->setStyleSheet(StyleSheet);
	//QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
	
	/*
	*	Show the window and return the execute callback from the QApplication.
	*/
	q_application_window.show();
    return q_application.exec();
}
