/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODEL_ParticleSystem_h_
#define _MODEL_ParticleSystem_h_

#ifdef _MODEL_MPI_
#include <ParticleSystemParallel.h>
#else
#include <ParticleSystemSerial.h>
#endif

//#include <SystemProperties.h>

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief Class template that stores _ParticleType objects organized
     * in SpatialCell(s) and computes nearest-neighbor and far-field interactions
     * among them. 
     * 
     * ParticleSystem can be compiled in serial and parallel version, depending
     * whether _MODEL_MPI_ is defined or not.
     */
    template <typename _ParticleType >
    struct ParticleSystem :
#ifdef _MODEL_MPI_
    /* inheritance */  public ParticleSystemParallel<_ParticleType>
#else
    /* inheritance */  public ParticleSystemSerial<_ParticleType>
#endif
    {
        
        /*****************************************/
        ParticleSystem()
        {/*!
          */
        }
        
#ifdef _MODEL_MPI_
        /*****************************************/
        ParticleSystem(int argc, char* argv[]) :
        /* init list */  ParticleSystemParallel<_ParticleType>(argc,argv)
        {/*!
          */
        }
#else
        /*****************************************/
        ParticleSystem(int argc, char* argv[])
        {/*! Dummy constructor, does nothing
          */
            argc+=0; // avoid unused warning
            argv=argv; // avoid unused warning

        }
        
        static void initMPI(int argc, char* argv[])
        {/*! Dummy static function, does nothing 
          */
            argc+=0; // avoid unused warning
            argv=argv; // avoid unused warning
        }
#endif
        
    };
    
} // end namespace
#endif

