/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TriangularMesh_H_
#define model_TriangularMesh_H_

#include <deque>
#include <Eigen/Dense>
#include <iostream>
#include <triangle.h>

namespace model
{
	
    class TriangularMesh : public std::deque<Eigen::Vector2d>
    /*                  */,public std::deque<Eigen::Vector3i>
    
    {
        

        
        
    public:
        
        const Eigen::Vector2d& vertex(const size_t& vID);
        std::deque<Eigen::Vector2d>& vertices();
        const std::deque<Eigen::Vector2d>& vertices() const;
        std::deque<Eigen::Vector3i>& triangles();
        const std::deque<Eigen::Vector3i>& triangles() const;
        void reMesh(const std::deque<Eigen::Matrix<double,2,1>>& boundaryPts,
                    const std::deque<Eigen::Matrix<double,2,1>>& internalPts,
                    const double& meshSize,
                    const std::string& flags);
        void reMesh(const std::vector<Eigen::Matrix<double,2,1>>& points,
                    const std::vector<Eigen::Matrix<int,2,1>>& segments,
                    const double& meshSize,
                    const std::string& flags);
        
        
    };
	
} // namespace model
#endif
