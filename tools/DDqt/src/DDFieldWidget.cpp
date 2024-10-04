/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDFieldWidget_cpp_
#define model_DDFieldWidget_cpp_

#include <sstream>      // std::ostringstream
#include <string>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <QVBoxLayout>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkAxis.h>
#include <DDFieldWidget.h>

namespace model
{

template <typename FloatType>
std::string to_string_exact(const FloatType& x, const int& n)
{
  std::ostringstream os;
  os << std::scientific<<std::setprecision(n)<< x;
  return os.str();
}

FieldDataPnt::FieldDataPnt(const DefectiveCrystal<3>& defectiveCrystal,const Eigen::Matrix<double,3,1>& Pin,const ElementType* const ele_in):
/* init */ voigtTraits(defectiveCrystal.ddBase.voigtTraits)
/* init */,P(Pin)
/* init */,ele(ele_in)
/* init */,displacement(defectiveCrystal.microstructures().size(),Eigen::Matrix<double,3,1>::Zero())
/* init */,stress(defectiveCrystal.microstructures().size(),Eigen::Matrix<double,3,3>::Zero())
/* init */,mobileConcentration(defectiveCrystal.microstructures().size(),Eigen::Matrix<double,ClusterDynamicsParameters<3>::mSize,1>::Zero())
/* init */,immobileClusters(defectiveCrystal.microstructures().size(),Eigen::Matrix<double,ClusterDynamicsParameters<3>::iSize,1>::Zero())
{
    
}

void FieldDataPnt::compute(const DefectiveCrystal<3>& defectiveCrystal)
{
    int mID(0);
    for(const auto& mStruct : defectiveCrystal.microstructures())
    {
        displacement[mID]=mStruct->displacement(P,nullptr,ele,nullptr);
        stress[mID]=mStruct->stress(P,nullptr,ele,nullptr);
        mobileConcentration[mID]=mStruct->mobileConcentration(P,nullptr,ele,nullptr);
        mID++;
    }
}


double FieldDataPnt::value(const int& valID,const std::vector<QCheckBox*>& microstructuresCheck) const
{
    //    std::cout<<valID<<std::endl;
    if(valID<3)
    {
        double temp(0.0);
        for(size_t k=0;k<displacement.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=displacement[k](valID);
            }
        }
        return temp;
    }
    else if(valID<6+3)
    {
        const int i(voigtTraits.tensorIndex(valID-3,0));
        const int j(voigtTraits.tensorIndex(valID-3,1));
        double temp(0.0);
        for(size_t k=0;k<stress.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=stress[k](i,j);
            }
        }
        return temp;
    }
    else if(valID==6+3)
    {
        Eigen::Matrix<double,3,3> temp(Eigen::Matrix<double,3,3>::Zero());
        for(size_t k=0;k<stress.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=stress[k];
            }
        }
        return temp.trace();
    }
    else if(valID==7+3)
    {
        Eigen::Matrix<double,3,3> temp(Eigen::Matrix<double,3,3>::Zero());
        for(size_t k=0;k<stress.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=stress[k];
            }
        }
        const Eigen::Matrix<double,3,3> tempDev(temp-temp.trace()/3.0*Eigen::Matrix<double,3,3>::Identity());
        return std::sqrt((tempDev*tempDev).trace()*1.5);
    }
//    else if(valID==8+3)
//    {
//        return solidAngle*useDD;
//    }
    else if(valID>7+3 && valID<=7+3+mSize)
    {
        double temp(0.0);
        for(size_t k=0;k<stress.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=mobileConcentration[k](valID-8-3);
            }
        }
        return temp;
    }
    else if(valID>7+3+mSize && valID<=7+3+mSize+iSize)
    {
        double temp(0.0);
        for(size_t k=0;k<stress.size();++k)
        {
            if(microstructuresCheck[k]->isChecked())
            {
                temp+=immobileClusters[k](valID-8-3-mSize);
            }
        }
        return temp;
    }
    else
    {
        return 0.0;
    }
}

DDFieldWidget::DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,
                             vtkRenderer* const renderer_in,
                             const  DefectiveCrystal<3>& defectiveCrystal_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,boxLabel(new QLabel(tr("# of planes")))
/* init */,spinBox(new QSpinBox(this))
/* init */,groupBox(new QGroupBox(tr("&Planes")))
/* init */,linesBoxLabel(new QLabel(tr("# of lines")))
/* init */,linesSpinBox(new QSpinBox(this))
/* init */,linesGroupBox(new QGroupBox(tr("&Lines")))
/* init */,computeButton(new QPushButton(tr("Compute")))
/* init */,fieldComboBox(new QComboBox(this))
///* init */,dislocationsCheck(new QCheckBox("dislocations",this))
///* init */,inclusionsCheck(new QCheckBox("inclusions",this))
///* init */,cdCheck(new QCheckBox("cluster dynamics",this))
///* init */,edCheck(new QCheckBox("elastic deformation",this))
/* init */,customScaleBox(new QGroupBox(tr("&Custom scale")))
/* init */,minScale(new QLineEdit(tr("0")))
/* init */,maxScale(new QLineEdit(tr("1")))
/* init */,scaleBarBox(new QGroupBox(tr("&ScaleBar")))
/* init */,lut(vtkSmartPointer<vtkLookupTable>::New())
/* init */,scalarBar(vtkSmartPointer<vtkScalarBarActor>::New())
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,defectiveCrystal(defectiveCrystal_in)
{
    
    QVBoxLayout* groupBoxLayout = new QVBoxLayout();
    groupBox->setLayout(groupBoxLayout);

    QVBoxLayout* linesGroupBoxLayout = new QVBoxLayout();
    linesGroupBox->setLayout(linesGroupBoxLayout);

    
    QGridLayout* autoscaleLayout = new QGridLayout();
    customScaleBox->setCheckable(true);
    customScaleBox->setChecked(false);
    
    autoscaleLayout->addWidget(minScale,0,0,1,1);
    autoscaleLayout->addWidget(maxScale,0,1,1,1);
    customScaleBox->setLayout(autoscaleLayout);
    
    scaleBarBox->setCheckable(true);
    scaleBarBox->setChecked(false);

    

    
    for(int k=0;k<3;++k)
    {
        fieldComboBox->insertItem(k,QString::fromStdString("displacement_"+std::to_string(k+1)));
    }
    for(int k=0;k<defectiveCrystal.ddBase.voigtTraits.tensorIndex.rows();++k)
    {
        const int& i(defectiveCrystal.ddBase.voigtTraits.tensorIndex(k,0));
        const int& j(defectiveCrystal.ddBase.voigtTraits.tensorIndex(k,1));
        fieldComboBox->insertItem(k+3,QString::fromStdString("stress_"+std::to_string(i+1)+std::to_string(j+1)));
    }
    fieldComboBox->insertItem(6+3,"tr(stress)");
    fieldComboBox->insertItem(7+3,"stress_VM");
    for(int k=0;k<ClusterDynamicsParameters<3>::mSize;++k)
    {
        fieldComboBox->insertItem(8+3+k,"v");
    }
    
    for(const auto& mstruct : defectiveCrystal.microstructures())
    {
        microstructuresCheck.emplace_back(new QCheckBox(QString::fromStdString(mstruct->tag),this));
        microstructuresCheck.back()->setChecked(true);
        connect(microstructuresCheck.back(),SIGNAL(stateChanged(int)), this, SLOT(plotField()));
    }
    
    connect(spinBox,SIGNAL(valueChanged(int)), this, SLOT(resetPlanes()));
    connect(linesSpinBox,SIGNAL(valueChanged(int)), this, SLOT(resetLines()));
    connect(computeButton,SIGNAL(released()), this, SLOT(compute()));
    connect(fieldComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(plotField()));
    connect(minScale,SIGNAL(returnPressed()), this, SLOT(plotField()));
    connect(maxScale,SIGNAL(returnPressed()), this, SLOT(plotField()));
//    connect(dislocationsCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
//    connect(inclusionsCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
//    connect(cdCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
//    connect(edCheck,SIGNAL(stateChanged(int)), this, SLOT(plotField()));
    connect(customScaleBox,SIGNAL(toggled(bool)), this, SLOT(plotField()));
    connect(scaleBarBox,SIGNAL(toggled(bool)), this, SLOT(plotField()));

    
    
    mainLayout->addWidget(boxLabel,0,0,1,1);
    mainLayout->addWidget(spinBox,0,1,1,1);
    mainLayout->addWidget(groupBox,1,0,1,2);
    
    mainLayout->addWidget(linesBoxLabel,2,0,1,1);
    mainLayout->addWidget(linesSpinBox,2,1,1,1);
    mainLayout->addWidget(linesGroupBox,3,0,1,2);
    
    mainLayout->addWidget(computeButton,4,0,1,1);
    mainLayout->addWidget(fieldComboBox,5,0,1,1);
    for(size_t k=0;k<microstructuresCheck.size();++k)
    {
        mainLayout->addWidget(microstructuresCheck[k],4+k,1,1,1);
    }
//    mainLayout->addWidget(dislocationsCheck,2,1,1,1);
//    mainLayout->addWidget(inclusionsCheck,3,1,1,1);
//    mainLayout->addWidget(cdCheck,4,1,1,1);
//    mainLayout->addWidget(edCheck,5,1,1,1);
    mainLayout->addWidget(customScaleBox,4+microstructuresCheck.size(),0,1,2);
    mainLayout->addWidget(scaleBarBox,5+microstructuresCheck.size(),0,1,2);

    this->setLayout(mainLayout);
    
    
    lut->SetHueRange(0.66667, 0.0);
    lut->Build();

    scalarBar->VisibilityOff();
    scalarBar->SetNumberOfLabels(4);
    scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
    scalarBar->SetLookupTable( lut );
    
    renderer->AddActor2D(scalarBar);

}

void DDFieldWidget::compute()
{
    if(groupBox)
    {
        if(groupBox->layout())
        {
            for(int k=0;k<groupBox->layout()->count();++k)
            {
                QLayoutItem *item = groupBox->layout()->itemAt(k);
                QWidget* widget = item->widget();
                if(widget)
                {
                    auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                    if (ddPlaneField)
                    {
                        ddPlaneField->compute(defectiveCrystal);
                    }
                }
            }
        }
    }
    
    if(linesGroupBox)
    {
        if(linesGroupBox->layout())
        {
            for(int k=0;k<linesGroupBox->layout()->count();++k)
            {
                QLayoutItem *item = linesGroupBox->layout()->itemAt(k);
                QWidget* widget = item->widget();
                if(widget)
                {
                    auto* ddLineField = dynamic_cast<DDLineField*>(widget);
                    if (ddLineField)
                    {
                        ddLineField->compute(defectiveCrystal);
                    }
                }
            }
        }
    }
    
    plotField();

}

void DDFieldWidget::plotField()
{
    
//        std::vector<bool> checkedMicrostructs;
//        for(const auto& mStruct : microstructuresCheck)
//        {
//            checkedMicrostructs.push_back(mStruct->isChecked());
//        }
        const int valID(fieldComboBox->currentIndex());
        
        // Compute lut range if customScaleBox is off
        if(!customScaleBox->isChecked())
        {// compute range
            double minValue=std::numeric_limits<double>::max();
            double maxValue=-std::numeric_limits<double>::max();
            if(groupBox->layout())
            {
                for(int k=0;k<groupBox->layout()->count();++k)
                {
                    QLayoutItem *item = groupBox->layout()->itemAt(k);
                    QWidget* widget = item->widget();
                    if(widget)
                    {
                        auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                        if (ddPlaneField)
                        {
                            if(ddPlaneField->groupBox->isChecked())
                            {
                                for(const auto& vtx : ddPlaneField->dataPnts())
                                {
                                    //                                const double value(vtx.value(valID,dislocationsCheck->isChecked(),inclusionsCheck->isChecked(),cdCheck->isChecked(),edCheck->isChecked()));
                                    const double value(vtx.value(valID,microstructuresCheck));
                                    minValue=std::min(minValue,value);
                                    maxValue=std::max(maxValue,value);
                                }
                            }
                        }
                    }
                }
        }
            
            if(linesGroupBox->layout())
            {
                for(int k=0;k<linesGroupBox->layout()->count();++k)
                {
                    QLayoutItem *item = linesGroupBox->layout()->itemAt(k);
                    QWidget* widget = item->widget();
                    if(widget)
                    {
                        auto* ddLineField = dynamic_cast<DDLineField*>(widget);
                        if (ddLineField)
                        {
                            if(ddLineField->groupBox->isChecked())
                            {
                                for(const auto& vtx : ddLineField->dataPnts())
                                {
                                    const double value(vtx.value(valID,microstructuresCheck));
                                    minValue=std::min(minValue,value);
                                    maxValue=std::max(maxValue,value);
                                }
                            }
                        }
                    }
                }
        }
            
            minScale->setText(QString::fromStdString(to_string_exact(minValue,7)));
            maxScale->setText(QString::fromStdString(to_string_exact(maxValue,7)));
        }
        
        // Set lut range
        try
        {
            double lutMin(std::stod(minScale->text().toStdString().c_str()));
            //            minScale->setStyleSheet("color: black");
            try
            {
                double lutMax(std::stod(maxScale->text().toStdString().c_str()));
                //                maxScale->setStyleSheet("color: black");
                lut->SetTableRange(lutMin, lutMax);
                
                for(int k=0;k<groupBox->layout()->count();++k)
                {
                    QLayoutItem *item = groupBox->layout()->itemAt(k);
                    QWidget* widget = item->widget();
                    if(widget)
                    {
                        auto* ddPlaneField = dynamic_cast<DDPlaneField*>(widget);
                        if (ddPlaneField)
                        {
                            if(ddPlaneField->groupBox->isChecked())
                            {
                                ddPlaneField->plotField(valID,microstructuresCheck,lut);
                            }
                        }
                    }
                }
                
                for(int k=0;k<linesGroupBox->layout()->count();++k)
                {
                    QLayoutItem *item = linesGroupBox->layout()->itemAt(k);
                    QWidget* widget = item->widget();
                    if(widget)
                    {
                        auto* ddLineField = dynamic_cast<DDLineField*>(widget);
                        if (ddLineField)
                        {
                            if(ddLineField->groupBox->isChecked())
                            {
                                ddLineField->plotField(valID,microstructuresCheck,fieldComboBox);
                            }
                        }
                    }
                }
                
                scalarBar->SetVisibility(scaleBarBox->isChecked());
                renWin->Render();
            }
            catch (const std::exception& e)
            {
                //                maxScale->setStyleSheet("color: red");
            }
        }
        catch (const std::exception& e)
        {
            //            minScale->setStyleSheet("color: red");
        }
    
}

void DDPlaneField::compute(const DefectiveCrystal<3>& defectiveCrystal)
{
    std::cout<<"DDPlaneField computing "<<dataPnts().size()<<" points..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
//    const Simplex<3,3>* gessSimplex(&defectiveCrystal.ddBase.mesh.simplices().begin()->second);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t vtkID=0;vtkID<dataPnts().size();++vtkID)
    {
        auto& vtx(dataPnts()[vtkID]);
        vtx.compute(defectiveCrystal);
    }
    std::cout<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void DDFieldWidget::clearLayout(QLayout *layout)
{
    if(layout)
    {
        while(layout->count() > 0)
        {
            QLayoutItem *item = layout->takeAt(0);
            QWidget* widget = item->widget();
            if(widget)
            {
                delete widget;
            }
            delete item;
        }
    }
}

void DDFieldWidget::resetPlanes()
{
    clearLayout(groupBox->layout());
    std::cout<<"Resetting "<<spinBox->value()<<" planes"<<std::endl;
    for(int k=0;k<spinBox->value();++k)
    {
        DDPlaneField* planeField(new DDPlaneField(renWin,renderer,defectiveCrystal));
        planeField->resetPlane();
        groupBox->layout()->addWidget(planeField);
    }
    renWin->Render();
}

void DDFieldWidget::resetLines()
{
    clearLayout(linesGroupBox->layout());
    std::cout<<"Resetting "<<linesSpinBox->value()<<" lines"<<std::endl;
    for(int k=0;k<linesSpinBox->value();++k)
    {
        DDLineField* lineField(new DDLineField(renWin,renderer,defectiveCrystal));
        lineField->resetLine();
        linesGroupBox->layout()->addWidget(lineField);
    }
    renWin->Render();
}

void DDPlaneField::plotField(const int& valID,const std::vector<QCheckBox*>& microstructuresCheck,const vtkSmartPointer<vtkLookupTable>& lut)
{
    vtkSmartPointer<vtkUnsignedCharArray> fieldColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
    fieldColors->SetNumberOfComponents(3);
    for(auto& vtx : dataPnts())
    {
        const double value(vtx.value(valID,microstructuresCheck));
        double dclr[3];
        lut->GetColor(value, dclr);
        unsigned char cclr[3];
        for(unsigned int j = 0; j < 3; j++)
        {
            cclr[j] = static_cast<unsigned char>(255.0 * dclr[j]);
        }
        fieldColors->InsertNextTypedTuple(cclr);
    }
    meshPolydata->GetPointData()->SetScalars(fieldColors);
    meshPolydata->Modified();
    meshMapper->SetScalarModeToUsePointData();
}

DDPlaneField::DDPlaneField(vtkGenericOpenGLRenderWindow* const renWin_in,
                           vtkRenderer* const renderer_in,
                           const DefectiveCrystal<3>& defectiveCrystal_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,groupBox(new QGroupBox(tr("&Plane")))
/* init */,posEdit(new QLineEdit("0 0 0"))
/* init */,normalEdit(new QLineEdit("1 1 1"))
/* init */,meshSizeEdit(new QLineEdit("1000"))
/* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,defectiveCrystal(defectiveCrystal_in)
///* init */,poly(poly_in)
{
    const VectorDim c(0.5*(defectiveCrystal.ddBase.mesh.xMax()+defectiveCrystal.ddBase.mesh.xMin()));
    posEdit->setText(QString::fromStdString(std::to_string(c(0))+" "+std::to_string(c(1))+" "+std::to_string(c(2))));
    
    const VectorDim n(Eigen::Matrix<double,3,1>::Random());
    normalEdit->setText(QString::fromStdString(std::to_string(n(0))+" "+std::to_string(n(1))+" "+std::to_string(n(2))));
    
    QGridLayout* groupBoxLayout = new QGridLayout();
    groupBoxLayout->addWidget(posEdit,0,0,1,1);
    groupBoxLayout->addWidget(normalEdit,1,0,1,1);
    groupBoxLayout->addWidget(meshSizeEdit,2,0,1,1);
    groupBox->setLayout(groupBoxLayout);
    
    groupBox->setCheckable(true);
    mainLayout->addWidget(groupBox,0,0,1,1);
    this->setLayout(mainLayout);
    
    connect(posEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(normalEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(meshSizeEdit,SIGNAL(returnPressed()), this, SLOT(resetPlane()));
    connect(groupBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
    
    meshPolydata->Allocate();
    meshMapper->SetInputData(meshPolydata);
    meshActor->SetMapper ( meshMapper );
    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.
    renderer->AddActor(meshActor);
}

DDPlaneField::~DDPlaneField()
{
    renderer->RemoveActor(meshActor);
}

const std::deque<FieldDataPnt>& DDPlaneField::dataPnts() const
{
    return *this;
}

std::deque<FieldDataPnt>& DDPlaneField::dataPnts()
{
    return *this;
}

void DDPlaneField::modify()
{
    meshActor->SetVisibility(groupBox->isChecked());
    renWin->Render();
}

void DDPlaneField::resetPlane()
{
    
    double meshSize;
    std::stringstream ssM(meshSizeEdit->text().toStdString());
    if(ssM >> meshSize)
    {
        meshSizeEdit->setStyleSheet("background-color: white");
        
        VectorDim P;
        std::stringstream ssP(posEdit->text().toStdString());
        if(ssP >> P(0) && ssP >> P(1) && ssP >> P(2))
        {
            posEdit->setStyleSheet("background-color: white");
//            std::cout<<"P="<<P.transpose()<<std::endl;
            
            VectorDim N;
            std::stringstream ssN(normalEdit->text().toStdString());
            if(ssN >> N(0) && ssN >> N(1) && ssN >> N(2))
            {
                normalEdit->setStyleSheet("background-color: white");
//                std::cout<<"N="<<N.transpose()<<std::endl;
                const double nNorm(N.norm());
                if(nNorm>FLT_EPSILON)
                {
                    plane.reset(new MeshPlane<3>(defectiveCrystal.ddBase.mesh,P,N));
                    std::deque<Eigen::Matrix<double,2,1>> boundaryPts;
                    std::deque<Eigen::Matrix<double,2,1>> internalPts;
                    for(const auto& bndLine : plane->meshIntersections)
                    {
                        boundaryPts.push_back(plane->localPosition(bndLine->P0));
                        
                    }
                    this->reMesh(boundaryPts,internalPts,meshSize,"pazq");
                    
                    dataPnts().clear();
                    vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
                    for(const auto& point2d : this->vertices())
                    {
                        const auto point3d(plane->globalPosition(point2d));
                        meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
                        const ElementType* ele(nullptr);
                        if(defectiveCrystal.ddBase.fe)
                        {
                            const auto searchPair(defectiveCrystal.ddBase.mesh.search(point3d));
                            if(searchPair.first)
                            {
                                ele=&defectiveCrystal.ddBase.fe->elements().at(searchPair.second->xID);
                            }
                        }
                        dataPnts().emplace_back(defectiveCrystal,point3d,ele);
                    }
                    
                    vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
                    vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
                    meshColors->SetNumberOfComponents(3);
                    for(const auto& tri : this->triangles())
                    {
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,tri(0));
                        triangle->GetPointIds()->SetId (1,tri(1));
                        triangle->GetPointIds()->SetId (2,tri(2));
                        meshTriangles->InsertNextCell ( triangle );
                        const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
                        meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
                    }
                    
                    meshPolydata->SetPoints ( meshPts );
                    meshPolydata->SetPolys ( meshTriangles );
                    meshPolydata->GetCellData()->SetScalars(meshColors);
                    meshPolydata->Modified();
                    meshMapper->SetScalarModeToUseCellData();
                    renWin->Render();
                }
                else
                {
                    normalEdit->setStyleSheet("background-color: red");
                }
            }
            else
            {
                normalEdit->setStyleSheet("background-color: red");
            }
        }
        else
        {
            posEdit->setStyleSheet("background-color: red");
        }
    }
    else
    {
        meshSizeEdit->setStyleSheet("background-color: red");
    }
}


DDLineField::DDLineField(vtkGenericOpenGLRenderWindow* const renWin_in,
                           vtkRenderer* const renderer_in,
                           const DefectiveCrystal<3>& defectiveCrystal_in):
/* init */ mainLayout(new QGridLayout(this))
/* init */,groupBox(new QGroupBox(tr("&Line")))
/* init */,posEdit(new QLineEdit("0 0 0"))
/* init */,directionEdit(new QLineEdit("1 0 0"))
/* init */,numPoinEdit(new QLineEdit("1000"))
///* init */,clrDialog(new QColorDialog(this))
///* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
///* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
///* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,renWin(renWin_in)
/* init */,renderer(renderer_in)
/* init */,chart(vtkSmartPointer<vtkChartXY>::New())
/* init */,chartScene(vtkSmartPointer<vtkContextScene>::New())
/* init */,chartActor(vtkSmartPointer<vtkContextActor>::New())
/* init */,table(vtkSmartPointer<vtkTable>::New())
/* init */,colors(vtkSmartPointer<vtkNamedColors>::New())
/* init */,points(chart->AddPlot(vtkChart::LINE))
/* init */,lineSource(vtkSmartPointer<vtkLineSource>::New())
/* init */,lineMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,lineActor(vtkSmartPointer<vtkActor>::New())
/* init */,defectiveCrystal(defectiveCrystal_in)
{
    const VectorDim c(0.5*(defectiveCrystal.ddBase.mesh.xMax()+defectiveCrystal.ddBase.mesh.xMin()));
    posEdit->setText(QString::fromStdString(std::to_string(c(0))+" "+std::to_string(c(1))+" "+std::to_string(c(2))));
    
    const VectorDim n(Eigen::Matrix<double,3,1>::Random());
    directionEdit->setText(QString::fromStdString(std::to_string(n(0))+" "+std::to_string(n(1))+" "+std::to_string(n(2))));
    
    QGridLayout* groupBoxLayout = new QGridLayout();
    groupBoxLayout->addWidget(posEdit,0,0,1,2);
    groupBoxLayout->addWidget(directionEdit,1,0,1,2);
    groupBoxLayout->addWidget(numPoinEdit,2,0,1,1);
//    groupBoxLayout->addWidget(clrDialog,2,1,1,1);
    groupBox->setLayout(groupBoxLayout);
    
    groupBox->setCheckable(true);
    mainLayout->addWidget(groupBox,0,0,1,1);
    this->setLayout(mainLayout);
    
    connect(posEdit,SIGNAL(returnPressed()), this, SLOT(resetLine()));
    connect(directionEdit,SIGNAL(returnPressed()), this, SLOT(resetLine()));
    connect(numPoinEdit,SIGNAL(returnPressed()), this, SLOT(resetLine()));
    connect(groupBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
    
    vtkColor3d color3d1 = colors->GetColor3d("Peacock");
    points->SetColorF(color3d1.GetRed(), color3d1.GetGreen(), color3d1.GetBlue());
    points->SetWidth(5.0);
    
    chart->GetAxis(0)->SetNumberOfTicks(10);
    chart->GetAxis(1)->SetNumberOfTicks(10);
    chart->GetAxis(0)->SetGridVisible(false);
    chart->GetAxis(1)->SetGridVisible(false);

    
    chartScene->AddItem(chart);
    chartActor->SetScene(chartScene);
    renderer->AddActor(chartActor);
    chartScene->SetRenderer(renderer);

    lineMapper->SetInputConnection(lineSource->GetOutputPort());
    lineActor->SetMapper(lineMapper);
    lineActor->GetProperty()->SetLineWidth(4);
    lineActor->GetProperty()->SetColor(colors->GetColor3d("Peacock").GetData());
    renderer->AddActor(lineActor);

//    meshPolydata->Allocate();
//    meshMapper->SetInputData(meshPolydata);
//    meshActor->SetMapper ( meshMapper );
//    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.
//    renderer->AddActor(meshActor);
}

DDLineField::~DDLineField()
{
    renderer->RemoveActor(chartActor);
    renderer->RemoveActor(lineActor);

}

const std::deque<FieldDataPnt>& DDLineField::dataPnts() const
{
    return *this;
}

std::deque<FieldDataPnt>& DDLineField::dataPnts()
{
    return *this;
}

const std::deque<double>& DDLineField::abscissa() const
{
    return *this;
}

std::deque<double>& DDLineField::abscissa()
{
    return *this;
}

void DDLineField::resetLine()
{
    
    int numPoints;
    std::stringstream ssM(numPoinEdit->text().toStdString());
    if(ssM >> numPoints)
    {
        if(numPoints>2)
        {
        numPoinEdit->setStyleSheet("background-color: white");
        
        VectorDim P;
        std::stringstream ssP(posEdit->text().toStdString());
        if(ssP >> P(0) && ssP >> P(1) && ssP >> P(2))
        {
            posEdit->setStyleSheet("background-color: white");
//            std::cout<<"P="<<P.transpose()<<std::endl;
            
            VectorDim N;
            std::stringstream ssN(directionEdit->text().toStdString());
            if(ssN >> N(0) && ssN >> N(1) && ssN >> N(2))
            {
                directionEdit->setStyleSheet("background-color: white");
//                std::cout<<"d="<<N.transpose()<<std::endl;
                const double nNorm(N.norm());
                if(nNorm>FLT_EPSILON)
                {
                    const VectorDim unitDir(N/nNorm);
                    line.reset(new MeshLine<3>(defectiveCrystal.ddBase.mesh,P,unitDir));
                    lineSource->SetPoint1(line->P0.data());
                    lineSource->SetPoint2(line->P1.data());
                    dataPnts().clear();
                    abscissa().clear();
                    const VectorDim chord(line->P1-line->P0);
                    const VectorDim dP(chord/(numPoints-1));
                    for(int k=0;k<numPoints;++k)
                    {
                        const VectorDim point3d(line->P0+k*dP);
                        const ElementType* ele(nullptr);
                        if(defectiveCrystal.ddBase.fe)
                        {
                            const auto searchPair(defectiveCrystal.ddBase.mesh.search(point3d));
                            if(searchPair.first)
                            {
                                ele=&defectiveCrystal.ddBase.fe->elements().at(searchPair.second->xID);
                            }
                        }
                        dataPnts().emplace_back(defectiveCrystal,point3d,ele);
                        abscissa().emplace_back((point3d-P).dot(unitDir));

                    }
                    
//                    vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
//                    for(const auto& point2d : this->vertices())
//                    {
//                        const auto point3d(plane->globalPosition(point2d));
//                        meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
//                        const ElementType* ele(nullptr);
//                        if(defectiveCrystal.ddBase.fe)
//                        {
//                            const auto searchPair(defectiveCrystal.ddBase.mesh.search(point3d));
//                            if(searchPair.first)
//                            {
//                                ele=&defectiveCrystal.ddBase.fe->elements().at(searchPair.second->xID);
//                            }
//                        }
//                        dataPnts().emplace_back(defectiveCrystal,point3d,ele);
//                    }
                    
//                    vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
//                    vtkSmartPointer<vtkUnsignedCharArray> meshColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
//                    meshColors->SetNumberOfComponents(3);
//                    for(const auto& tri : this->triangles())
//                    {
//                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
//                        triangle->GetPointIds()->SetId (0,tri(0));
//                        triangle->GetPointIds()->SetId (1,tri(1));
//                        triangle->GetPointIds()->SetId (2,tri(2));
//                        meshTriangles->InsertNextCell ( triangle );
//                        const auto triColor(Eigen::Matrix<int,1,3>::Random()*255);
//                        meshColors->InsertNextTuple3(triColor(0),triColor(1),triColor(2)); // use this to assig color to each vertex
//                    }
//                    
//                    meshPolydata->SetPoints ( meshPts );
//                    meshPolydata->SetPolys ( meshTriangles );
//                    meshPolydata->GetCellData()->SetScalars(meshColors);
//                    meshPolydata->Modified();
//                    meshMapper->SetScalarModeToUseCellData();
//                    renWin->Render();
                }
                else
                {
                    directionEdit->setStyleSheet("background-color: red");
                }
            }
            else
            {
                directionEdit->setStyleSheet("background-color: red");
            }
        }
        else
        {
            posEdit->setStyleSheet("background-color: red");
        }
    }
    else
    {
        numPoinEdit->setStyleSheet("background-color: red");
    }
    }
    else
    {
        numPoinEdit->setStyleSheet("background-color: red");
    }
    modify();
}



void DDLineField::compute(const DefectiveCrystal<3>&)
{
    std::cout<<"DDLineField computing "<<dataPnts().size()<<" points..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t vtkID=0;vtkID<dataPnts().size();++vtkID)
    {
        auto& vtx(dataPnts()[vtkID]);
        vtx.compute(defectiveCrystal);
    }
    std::cout<<magentaColor<<"["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void DDLineField::plotField(const int& valID,const std::vector<QCheckBox*>& microstructuresCheck,const QComboBox* const fieldComboBox)
{
    
    table->Initialize();
    vtkNew<vtkFloatArray> xArray;
    xArray->SetName("x");
    table->AddColumn(xArray);
    vtkNew<vtkFloatArray> yArray;
    yArray->SetName(fieldComboBox->itemText(valID).toStdString().c_str());
    table->AddColumn(yArray);

    for(size_t k=0;k<abscissa().size();++k)
    {                    
        const auto row=table->InsertNextBlankRow(2);
        table->SetValue(row, 0, abscissa()[k]);
        table->SetValue(row, 1, dataPnts()[k].value(valID,microstructuresCheck));
    }
    
    

    
    points->SetInputData(table, 0, 1);
    table->Modified();
    chart->RecalculateBounds();
    
    //    chart->GetAxis(0)->SetTitle(yComboBox->currentText().toStdString());
        chart->GetAxis(0)->SetTitle(fieldComboBox->itemText(valID).toStdString());

    //    chart->GetAxis(0)->SetLogScale(yAxisLog->isChecked());
    //    chart->GetAxis(1)->SetLogScale(xAxisLog->isChecked());


//        chart->GetAxis(0)->SetNumberOfTicks(10);
//        chart->GetAxis(1)->SetNumberOfTicks(10);

        chart->GetAxis(0)->GetLabelProperties()->SetFontSize(16);
        chart->GetAxis(1)->GetLabelProperties()->SetFontSize(16);

        chart->GetAxis(0)->GetTitleProperties()->SetFontSize(16);
        chart->GetAxis(1)->GetTitleProperties()->SetFontSize(16);
    
    chart->Modified();
}

void DDLineField::modify()
{
    lineActor->SetVisibility(groupBox->isChecked());
    chartActor->SetVisibility(groupBox->isChecked());
    renWin->Render();

}

} // namespace model
#endif
