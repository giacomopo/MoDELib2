/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLoopActor_h_
#define model_NetworkLoopActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QSlider>
#include <QColorDialog>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QLineEdit>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>
#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <PeriodicGlidePlaneFactory.h>
#include <DefectiveCrystal.h>

namespace model
{
    struct NetworkLoopActor : public QWidget
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;

        QGridLayout* mainLayout;
        QCheckBox* showLoops;
        
        QGroupBox* slippedAreaBox;
        QSlider* sliderSlippedArea;

        QGroupBox* meshAreaBox;

        vtkSmartPointer<vtkPolyData> loopPolyData;
        vtkSmartPointer<vtkPolyDataMapper> loopMapper;
        vtkSmartPointer<vtkActor> loopActor;

        vtkSmartPointer<vtkPolyData> areaPolyData;
        vtkSmartPointer<vtkTriangleFilter> areaTriangleFilter; // needed to handle non-convex area polygons
        vtkSmartPointer<vtkPolyDataMapper> areaMapper;
        vtkSmartPointer<vtkActor> areaActor;
        
        vtkSmartPointer<vtkPolyData> meshPolydata;
        vtkSmartPointer<vtkPolyDataMapper> meshMapper;
        vtkSmartPointer<vtkActor> meshActor;
        QLineEdit* globalMeshSizeEdit;
        QLineEdit* localMeshSizeEdit;

        
        public:
                
        const DefectiveCrystal<3>& defectiveCrystal;
        const std::shared_ptr<DislocationNetwork<3,0>> dislocationNetwork;

        NetworkLoopActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,const DefectiveCrystal<3>& defectiveCrystal_in);
        void updateConfiguration();
        
    };
    
} // namespace model
#endif
