/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalActor_cpp_
#define model_DefectiveCrystalActor_cpp_

#include <algorithm>
#include <QString>
#include <QApplication>
#include <QFuture>
#include <QtCore>
#include <QtConcurrent/QtConcurrentRun>
#include <QImage>

#include <DefectiveCrystalActor.h>
#include <SimplicialMesh.h>

namespace model
{
        
        /**********************************************************************/
        DefectiveCrystalActor::DefectiveCrystalActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,
                                 QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                                 DislocationDynamicsBase<3>& ddBase_in) :
        /* init */ DDconfigIO<3>(ddBase_in.simulationParameters.traitsIO.evlFolder)
        /* init */,DDauxIO<3>(ddBase_in.simulationParameters.traitsIO.auxFolder)
        /* init */,renderWindow(renWin)
        /* init */,qvtkGLwidget(qvtkGLwidget_in)
        /* init */,ddBase(ddBase_in)
//        /* init */,configFields(ddBase,*this)
        /* init */,defectiveCrystal(ddBase)
        /* init */,nodes(new NetworkNodeActor(renWin,ren))
        /* init */,segments(new NetworkLinkActor(renWin,ren,defectiveCrystal))
        /* init */,loops(new NetworkLoopActor(renWin,ren,defectiveCrystal))
        /* init */,inclusions(new InclusionActor(renWin,ren,defectiveCrystal))
//        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,defectiveCrystal))
        /* init */,quadrature(new QuadratureActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,chartActor(new ChartActor(ddBase.simulationParameters.traitsIO,renderWindow,ren))
        /* init */,ddField(new DDFieldWidget(renderWindow,ren,defectiveCrystal))
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,frameIDedit(new QLineEdit("0"))
        /* init */,plusFrameButton(new QPushButton(">"))
        /* init */,playFrameButton(new QPushButton(">>"))
        /* init */,minusFrameButton(new QPushButton("<"))
        /* init */,frameIncrementEdit(new QLineEdit("1"))
        /* init */,saveImage(new QCheckBox(this))
        /* init */,tabWidget(new QTabWidget())
        {
            
            tabWidget->addTab(nodes, tr(std::string("Nodes").c_str()));
            tabWidget->addTab(segments, tr(std::string("Segments").c_str()));
            tabWidget->addTab(loops, tr(std::string("Loops").c_str()));
            tabWidget->addTab(inclusions, tr(std::string("Inclusions").c_str()));
            tabWidget->addTab(glidePlanes, tr(std::string("GlidePlanes").c_str()));
            tabWidget->addTab(quadrature, tr(std::string("Quadrature").c_str()));
            tabWidget->addTab(chartActor, tr(std::string("Chart").c_str()));
            tabWidget->addTab(ddField, tr(std::string("Fields").c_str()));

            saveImage->setText("save PNG");


            mainLayout->addWidget(frameIDedit,0,0,1,1);
            mainLayout->addWidget(minusFrameButton,0,1,1,1);
            mainLayout->addWidget(plusFrameButton,0,2,1,1);
            mainLayout->addWidget(playFrameButton,0,3,1,1);
            mainLayout->addWidget(frameIncrementEdit,0,4,1,1);
            mainLayout->addWidget(saveImage,1,0,1,1);


            mainLayout->addWidget(tabWidget,2,0,1,5);
            this->setLayout(mainLayout);

            connect(frameIDedit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));
            connect(plusFrameButton,SIGNAL(pressed()), this, SLOT(nextConfiguration()));
            connect(minusFrameButton,SIGNAL(pressed()), this, SLOT(prevConfiguration()));
            connect(playFrameButton,SIGNAL(pressed()), this, SLOT(playConfigurations()));

            connect(frameIDedit,SIGNAL(returnPressed()), ddField, SLOT(compute()));
            connect(plusFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));
            connect(minusFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));
            connect(playFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));

            
            QApplication::processEvents();
            updateConfiguration(0);
        }

        const DDconfigIO<3>& DefectiveCrystalActor::configIO() const
        {
            return *this;
        }

        DDconfigIO<3>& DefectiveCrystalActor::configIO()
        {
            return *this;
        }

        const DDauxIO<3>& DefectiveCrystalActor::auxIO() const
        {
            return *this;
        }

        DDauxIO<3>& DefectiveCrystalActor::auxIO()
        {
            return *this;
        }

        void DefectiveCrystalActor::playConfigurations()
        {
//            bool updated(true);
//            while (updated)
//            {
//                updated=nextConfiguration();
////                QFuture<bool> future = QtConcurrent::run(&DefectiveCrystalActor::nextConfiguration,this);
////                QFuture<bool> future = QtConcurrent::run(nextConfiguration);
////                future.waitForFinished();
////                updated=future.result();
////                QApplication::processEvents();
////                std::chrono::seconds dura( 1);
////                std::this_thread::sleep_for( dura );
////                pause(1);
//            }
        }

        bool DefectiveCrystalActor::nextConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID+currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            const bool updated(updateConfiguration());
            if(!updated)
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
            return updated;
        }

        bool DefectiveCrystalActor::prevConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID-currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            const bool updated(updateConfiguration());
            if(!updated)
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
            return updated;
        }

        bool DefectiveCrystalActor::updateConfiguration()
        {
            const size_t frameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            if(updateConfiguration(frameID))
            {
                renderWindow->Render();
                return true;
            }
            else
            {
                return false;
            }
        }

        bool DefectiveCrystalActor::updateConfiguration(const size_t& frameID)
        {
            try
            {
                DislocationNode<3,0>::force_count(0);
                DislocationLoopNode<3,0>::force_count(0);
                DislocationLoop<3,0>::force_count(0);
                EshelbyInclusionBase<3>::force_count(0);

                configIO().read(frameID);
                auxIO().read(frameID);
                //configFields.updateConfiguration();
                defectiveCrystal.initializeConfiguration(*this);
                nodes->updateConfiguration(*this);
                segments->updateConfiguration(nodes->nodePolyData);
                loops->updateConfiguration();
                inclusions->updateConfiguration();
//                glidePlanes->updateConfiguration(*this);
                glidePlanes->updateConfiguration();
                quadrature->updateConfiguration(*this);
                chartActor->updateConfiguration(frameID);
//                ddField->compute(); // this causes segFault at DDFieldWidget::plotField::renWin->Render();
                
                if(saveImage->isChecked())
                {
                    QImage img=qvtkGLwidget->grabFramebuffer();
                    img.save(QString::fromStdString(ddBase.simulationParameters.traitsIO.evlFolder+"/img_"+std::to_string(frameID)+".png"), "PNG", -1);
                }
                return true;
            }
            catch(const std::exception& e)
            {
                std::cout<<e.what()<<std::endl;
                return false;
            }
        }
        
        void DefectiveCrystalActor::modify()
        {
//            nodes.modify();
//            segments.modify();
//            loops.modify();
        }
        
}
#endif
