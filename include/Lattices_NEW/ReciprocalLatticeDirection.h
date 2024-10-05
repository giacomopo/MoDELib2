/* This file is part of model.
 *
 * MoDELib is distributed without any warranty under the MIT License.
 */


#ifndef model_ReciprocalLatticeDirection_h_
#define model_ReciprocalLatticeDirection_h_

#include <LatticeModule.h>

namespace model
{
    template <int dim>
    struct ReciprocalLatticeDirection :
    /* inherits */ public ReciprocalLatticeVector<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;


        ReciprocalLatticeDirection(const ReciprocalLatticeVector<dim>& v) ;
                
    };
    
} // end namespace
#endif
