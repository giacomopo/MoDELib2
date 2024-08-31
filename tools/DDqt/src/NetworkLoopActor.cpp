/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLoopActor_cpp_
#define model_NetworkLoopActor_cpp_

#include <iostream>
#include <deque>
#include <string>

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
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkLine.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkLoopActor.h>
#include <DislocationLoopPatches.h>

namespace model
{


/**********************************************************************/
NetworkLoopActor::NetworkLoopActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,
                                   const DefectiveCrystal<3>& defectiveCrystal_in) :
/* init */ renderWindow(renWin)
/* init */,mainLayout(new QGridLayout(this))
/* init */,showLoops(new QCheckBox(this))
/* init */,slippedAreaBox(new QGroupBox(tr("&Slipped Area")))
/* init */,sliderSlippedArea(new QSlider(this))
/* init */,meshAreaBox(new QGroupBox(tr("&Slipped Mesh")))
/* init */,loopPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,loopMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,loopActor(vtkSmartPointer<vtkActor>::New())
/* init */,areaPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,areaTriangleFilter(vtkSmartPointer<vtkTriangleFilter>::New())
/* init */,areaMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,areaActor(vtkSmartPointer<vtkActor>::New())
/* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,globalMeshSizeEdit(new QLineEdit("100000000"))
/* init */,localMeshSizeEdit(new QLineEdit("1.0"))
/* init */,defectiveCrystal(defectiveCrystal_in)
/* init */,dislocationNetwork(defectiveCrystal.template getUniqueTypedMicrostructure<DislocationNetwork<3,0>>())
{
    showLoops->setChecked(false);
    showLoops->setText("Loops");
    
    slippedAreaBox->setCheckable(true);
    slippedAreaBox->setChecked(false);
    
    meshAreaBox->setCheckable(true);
    meshAreaBox->setChecked(false);

    
    //            showSlippedArea->setChecked(false);
    //            showSlippedArea->setText("SlippedArea");
    //            sliderSlippedArea->setEnabled(false);
    sliderSlippedArea->setMinimum(0);
    sliderSlippedArea->setMaximum(10);
    sliderSlippedArea->setValue(5);
    sliderSlippedArea->setOrientation(Qt::Horizontal);
    
    //            colorSlippedArea->setWindowFlags(Qt::Widget );
    
    QVBoxLayout *slippedAreaLayout = new QVBoxLayout();
    slippedAreaLayout->addWidget(sliderSlippedArea);
    //            slippedAreaLayout->addWidget(colorSlippedArea);
    slippedAreaBox->setLayout(slippedAreaLayout);
    
    /* a few options that we must set for it to work nicely */
    //            colorSlippedArea->setOptions(
    //                            /* do not use native dialog */
    //                            QColorDialog::DontUseNativeDialog
    //                            /* you don't need to set it, but if you don't set this
    //                                the "OK" and "Cancel" buttons will show up, I don't
    //                                think you'd want that. */
    ////                            | QColorDialog::NoButtons
    //                );
    
    //            QGroupBox *groupBox = new QGroupBox(tr("&Push Buttons"));
    //                groupBox->setCheckable(true);
    //                groupBox->setChecked(true);
    
    QGridLayout *meshAreaLayout = new QGridLayout();
    meshAreaLayout->addWidget(globalMeshSizeEdit,0,0,1,1);
    meshAreaLayout->addWidget(localMeshSizeEdit,1,0,1,1);
    meshAreaBox->setLayout(meshAreaLayout);

    mainLayout->addWidget(showLoops,0,0,1,1);
    mainLayout->addWidget(slippedAreaBox,1,0,1,1);
    mainLayout->addWidget(meshAreaBox,2,0,1,1);

    //            mainLayout->addWidget(sliderSlippedArea,1,1,1,1);
    //            mainLayout->addWidget(colorSlippedArea,1,2,1,1);
    
    
    this->setLayout(mainLayout);
    
    connect(showLoops,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    //            connect(slippedAreaBox,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    connect(slippedAreaBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
    
    connect(sliderSlippedArea,SIGNAL(valueChanged(int)), this, SLOT(modify()));
//    connect(globalMeshSizeEdit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));
//    connect(localMeshSizeEdit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));

    
    loopMapper->SetInputData(loopPolyData);
    loopActor->SetMapper(loopMapper);
    loopActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
    loopActor->SetVisibility(false);
    
    areaTriangleFilter->SetInputData(areaPolyData);
    //            filter->SetInputData(polygonPolyData);
    //            areaTriangleFilter->Update();
    
    //            areaMapper->SetInputData(areaPolyData);
    areaMapper->SetInputData(areaTriangleFilter->GetOutput());
    areaActor->SetMapper(areaMapper);
    areaActor->GetProperty()->SetColor(0.0, 1.0, 1.0); //(R,G,B)
    areaActor->SetVisibility(false);
    
    meshPolydata->Allocate();
    meshMapper->SetInputData(meshPolydata);
    meshActor->SetMapper ( meshMapper );
    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.

    
    renderer->AddActor(loopActor);
    renderer->AddActor(areaActor);
    renderer->AddActor(meshActor);

    //            std::vector<Eigen::Vector2d> clippedPolygonTemp(SutherlandHodgman::clip(polygon2d,box2d));
    
}



void NetworkLoopActor::updateConfiguration()
{
    std::cout<<"Updating loops..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    if(dislocationNetwork)
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//        std::map<size_t,size_t> loopNodesMap; // nodeID,nodePositionInDDauxIO
//        size_t k=0;
        for(const auto& weakNode : dislocationNetwork->loopNodes())
        {
            const auto& loopNode(weakNode.second.lock());
  //          loopNodesMap.emplace(loopNode->sID,k);
            const auto& P(loopNode->get_P());
            points->InsertNextPoint(P(0),P(1),P(2));
  //          k++;
        }
        
        vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
        for(const auto& loopLink : dislocationNetwork->loopLinks())
        {
            vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
            line->GetPointIds()->SetId(0, loopLink.second.source->networkID()); // the second 0 is the index of the Origin in linesPolyData's points
            line->GetPointIds()->SetId(1, loopLink.second.  sink->networkID());
            cells->InsertNextCell(line);

            
//            const auto sourceIter(loopNodesMap.find(loopLink.first.first));
//            if(sourceIter!=loopNodesMap.end())
//            {
//                const auto sinkIter(loopNodesMap.find(loopLink.first.second));
//                if(sinkIter!=loopNodesMap.end())
//                {
//                    vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
//                    line->GetPointIds()->SetId(0, sourceIter->second); // the second 0 is the index of the Origin in linesPolyData's points
//                    line->GetPointIds()->SetId(1, sinkIter->second);
//                    cells->InsertNextCell(line);
//                }
//                else
//                {
//                    throw std::runtime_error("Sink vertex not found in nodeMap");
//                }
//            }
//            else
//            {
//                throw std::runtime_error("Source vertex not found in nodeMap");
//            }
        }
        loopPolyData->SetPoints(points);
        loopPolyData->SetLines(cells);
        loopPolyData->Modified();
        
        // Slipped area
        vtkNew<vtkPoints> areaPoints;
        vtkNew<vtkCellArray> areaPolygons;
        size_t areaPointID(0);
        for(const auto& weakloop : dislocationNetwork->loops())
        {
            const auto& loop(weakloop.second.lock());
            for(const auto& pair : loop->patches().globalPatches())
            {
                vtkNew<vtkPolygon> polygon;
                for(const auto& globalPos : pair.second)
                {
                    areaPoints->InsertNextPoint(globalPos(0),
                                                globalPos(1),
                                                globalPos(2));
                    polygon->GetPointIds()->InsertNextId(areaPointID);
                    areaPointID++;
                }
                areaPolygons->InsertNextCell(polygon);
            }
        }
        areaPolyData->SetPoints(areaPoints);
        areaPolyData->SetPolys(areaPolygons);
        areaTriangleFilter->Update();
        
        
        if(meshAreaBox->isChecked())
        {
            vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
            vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
            meshColors->SetNumberOfComponents(3);

            const double meshSize(std::atof(globalMeshSizeEdit->text() .toStdString().c_str()));
            const double localMeshSize(std::atof(localMeshSizeEdit->text() .toStdString().c_str()));
            size_t ptsIncrement(0);
            for(const auto& weakloop : dislocationNetwork->loops())
            {
                const auto& loop(weakloop.second.lock());
                for(const auto& patchMesh : loop->meshed(meshSize,localMeshSize))
                {
                    for(const auto& point3d : patchMesh.points)
                    {
                        meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
                    }
                    for(const auto& tri : patchMesh.triangles)
                    {
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,tri(0)+ptsIncrement);
                        triangle->GetPointIds()->SetId (1,tri(1)+ptsIncrement);
                        triangle->GetPointIds()->SetId (2,tri(2)+ptsIncrement);
                        meshTriangles->InsertNextCell ( triangle );
                        const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
                        meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
//                        meshColors->InsertNextTuple3(153,153,255); // use this to assig color to each vertex

                    }
                    ptsIncrement+=patchMesh.points.size();
                }
            }
            meshPolydata->SetPoints ( meshPts );
            meshPolydata->SetPolys ( meshTriangles );
            meshPolydata->GetCellData()->SetScalars(meshColors);
            meshPolydata->Modified();
            meshMapper->SetScalarModeToUseCellData();
//            renWin->Render();

            
        }
//        renderWindow->Render();
    }
    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void NetworkLoopActor::modify()
{    
    loopActor->SetVisibility(showLoops->isChecked());
    areaActor->SetVisibility(slippedAreaBox->isChecked());
    areaActor->GetProperty()->SetOpacity(sliderSlippedArea->value()/10.0);
    
    renderWindow->Render();
}

} // namespace model
#endif
