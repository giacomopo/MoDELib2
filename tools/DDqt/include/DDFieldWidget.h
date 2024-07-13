/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDFieldWidget_H_
#define model_DDFieldWidget_H_

#include <map>

#include <Eigen/Dense>

#include <memory>
#include <QGroupBox>
#include <QGridLayout>
#include <QWidget>
#include <QSpinBox>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkPointData.h>


#include <DDtraitsIO.h>
#include <SimplicialMesh.h>
#include <MeshPlane.h>
#include <TriangularMesh.h>
#include <DDconfigIO.h>
#include <Polycrystal.h>
#include <DefectiveCrystal.h>

namespace model
{

    struct FieldDataPnt
    {
        
        typedef  typename MicrostructureBase<3>::ElementType ElementType;

        static constexpr int mSize=ClusterDynamicsParameters<3>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<3>::iSize;
        const SymmetricVoigtTraits<3>& voigtTraits;
        const Eigen::Matrix<double,3,1> P;
        const ElementType* const ele;
        std::vector<Eigen::Matrix<double,3,1>> displacement;
        std::vector<Eigen::Matrix<double,3,3>> stress;
        std::vector<Eigen::Matrix<double,ClusterDynamicsParameters<3>::mSize,1>> mobileConcentration;
        std::vector<Eigen::Matrix<double,ClusterDynamicsParameters<3>::iSize,1>> immobileClusters;
        
        FieldDataPnt(const DefectiveCrystal<3>& defectiveCrystal,const Eigen::Matrix<double,3,1>& Pin,const ElementType* const ele_in);
        double value(const int& valID,const std::vector<QCheckBox*>& microstructuresCheck) const;
    };


    struct DDPlaneField : public QWidget
    /*                */, public TriangularMesh
    /*                */, public std::deque<FieldDataPnt>
    {
        
        Q_OBJECT
        
    public:
        typedef  typename MicrostructureBase<3>::ElementType ElementType;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        QGridLayout* mainLayout;
        QGroupBox* groupBox;
        
        QLineEdit* posEdit;
        QLineEdit* normalEdit;
        QLineEdit* meshSizeEdit;
        
        vtkSmartPointer<vtkPolyData> meshPolydata;
        vtkSmartPointer<vtkPolyDataMapper> meshMapper;
        vtkSmartPointer<vtkActor> meshActor;
        
        vtkGenericOpenGLRenderWindow* const renWin;
        vtkRenderer* const renderer;
        const DefectiveCrystal<3>& defectiveCrystal;
        std::shared_ptr<MeshPlane<3>> plane;
        
    public slots:
        void resetPlane();
        void modify();
        
    public:
        DDPlaneField(vtkGenericOpenGLRenderWindow* const renWin_in,vtkRenderer* const renderer_in,const DefectiveCrystal<3>& defectiveCrystal_in);
        ~DDPlaneField();
        const std::deque<FieldDataPnt>& dataPnts() const;
        std::deque<FieldDataPnt>& dataPnts();
        void compute(const DefectiveCrystal<3>&);
        void plotField(const int& valID,const std::vector<QCheckBox*>& microstructuresCheck,const vtkSmartPointer<vtkLookupTable>& lut);
    };

    struct DDFieldWidget : public QWidget
    {
        
        Q_OBJECT
        
        typedef  typename MicrostructureBase<3>::ElementType ElementType;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        QGridLayout* mainLayout;
        QLabel* boxLabel;
        QSpinBox* spinBox;
        QGroupBox* groupBox;
        QPushButton* computeButton;
        QComboBox* fieldComboBox;
        
        std::vector<QCheckBox*> microstructuresCheck;

        QGroupBox* customScaleBox;
        QLineEdit* minScale;
        QLineEdit* maxScale;
        QGroupBox* scaleBarBox;
        vtkSmartPointer<vtkLookupTable> lut;
        vtkSmartPointer<vtkScalarBarActor> scalarBar;

        vtkGenericOpenGLRenderWindow* const renWin;
        vtkRenderer* const renderer;
        const DefectiveCrystal<3>& defectiveCrystal;
        
    private slots:
        void clearLayout(QLayout *layout);
        void resetPlanes();
        void plotField();

   public slots:
        void compute();

    public:
        
        DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,
                      vtkRenderer* const renderer_in,
                      const  DefectiveCrystal<3>& defectiveCrystal_in);
    };

} // namespace model
#endif







