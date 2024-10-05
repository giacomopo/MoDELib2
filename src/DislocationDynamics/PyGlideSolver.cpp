/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * Written in 2024 by Matthew Maron <mlm335@miami.edu>
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PyGlideSolver_cpp_
#define model_PyGlideSolver_cpp_
#include <random>
#include <filesystem>
#include <PyGlideSolver.h>
#include <TerminalColors.h>
#include <Eigen/StdVector>
#include <vector>
#include <list>


// Functions for Adjusting Lists and Outputting Data
void printNestedList(std::list<std::list<int>> connections) // Print the connections
{
    std::cout << "[\n";
    std::list<std::list<int> >::iterator connections_itr;
    for (connections_itr = connections.begin();
        connections_itr != connections.end();
        ++connections_itr)
    {
        std::cout << " [";
        std::list<int>::iterator single_list_itr;
        std::list<int>& single_list_pointer = *connections_itr;
        for (single_list_itr = single_list_pointer.begin();
            single_list_itr != single_list_pointer.end();
            single_list_itr++)
        {
            std::cout << " " << *single_list_itr << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]";
}

void printNestedListToFile(std::list<std::list<int>>& connections, std::ofstream& myFile)
{
    std::list<std::list<int>>::iterator connections_itr;
    myFile << "[\n";
    for (connections_itr = connections.begin();
        connections_itr != connections.end();
        ++connections_itr)
    {
        myFile << " [";
        std::list<int>::iterator single_list_itr;
        std::list<int>& single_list_pointer = *connections_itr;
        for (single_list_itr = single_list_pointer.begin();
            single_list_itr != single_list_pointer.end();
            single_list_itr++)
        {
            myFile << " " << *single_list_itr << " ";
        }
        myFile << "]\n";
    }
    myFile << "]";
}

void sequenceOuterConnections(std::list<std::list<int>>& connections)
{
    std::set<int> unique_numbers;
    for (const auto& sublist : connections)
    {
        if (!sublist.empty())
        {
            unique_numbers.insert(sublist.front());
            unique_numbers.insert(sublist.back());
        }
    }
    std::vector<int> sorted_numbers(unique_numbers.begin(), unique_numbers.end());
    for (size_t i = 1; i < sorted_numbers.size(); ++i)
    {
        if (sorted_numbers[i] != sorted_numbers[i - 1] + 1)
        {
            int diff = sorted_numbers[i] - (sorted_numbers[i - 1] + 1);
            for (size_t j = i; j < sorted_numbers.size(); ++j)
            {
                sorted_numbers[j] -= diff;
            }
        }
    }
    std::map<int, int> adjust_map;
    for (size_t i = 0; i < sorted_numbers.size(); ++i)
    {
        adjust_map[*next(unique_numbers.begin(), i)] = sorted_numbers[i];
    }
    for (auto& sublist : connections)
    {
        for (auto& num : sublist)
        {
            if (adjust_map.find(num) != adjust_map.end())
            {
                num = adjust_map[num];
            }
        }
    }
}

void adjustQPConnections(std::list<std::list<int>>& connections)
{
    // First pass:
    int largest_outer_num = std::numeric_limits<int>::min();
    for (const auto& sublist : connections)
    {
        largest_outer_num = std::max(largest_outer_num, sublist.front());
        largest_outer_num = std::max(largest_outer_num, sublist.back());
    }
    int current_value = largest_outer_num + 1;
    // Second pass:
    for (auto& sublist : connections)
    {
        if (sublist.size() > 2)
        {
            for (auto it = std::next(sublist.begin()); it != std::prev(sublist.end()); ++it)
            {
                *it = current_value++;
            }
        }
    }
}


int findMin(const std::list<std::list<int>>& connections)
{
    int minVal = std::numeric_limits<int>::max();
    for (const auto& sublist : connections) 
    {
        for (const auto& val : sublist) 
        {
            if (val < minVal) 
            {
                minVal = val;
            }
        }
    }
    return minVal;
}
void adjustMinValueToZero(std::list<std::list<int>>& connections)
{
    int minVal = findMin(connections);
    if (minVal == 0) 
    {
        return;
    }
    for (auto& sublist : connections) {
        for (auto& val : sublist) {
            val -= minVal;
        }
    }
}


namespace model
{

    #ifdef _MODEL_PYBIND11_ // COMPLIED WITH PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    /* init */,pyModuleName(this->DN.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/"+TextFileParser(this->DN.ddBase.simulationParameters.traitsIO.ddFile).readString("pyModuleName"))
    {
        std::filesystem::path pyModulePath(pyModuleName);
        std::string pyModuleDir(pyModulePath.parent_path());
        std::string pyModuleFile(pyModulePath.filename());
        std::cout<<greenBoldColor<<"Creating PyGlideSolver: dir= "<<pyModuleDir<<", file="<<pyModuleFile<<defaultColor<<std::endl;
        //pybind11::scoped_interpreter guard{};
        pybind11::module sys = pybind11::module::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(pyModuleDir.c_str());
        pyModule=pybind11::module::import(pyModuleFile.c_str());
    }

    template <typename DislocationNetworkType>
    Eigen::MatrixXd PyGlideSolver<DislocationNetworkType>::updateBoundaryNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& boundaryNodeIDmap) const
    {
        Eigen::MatrixXd updatedNodePythonVelocity = nodePythonVelocity;
        for (const auto& pair : boundaryNodeIDmap)
        {
            int rowIndex = pair.second;
            if (rowIndex < updatedNodePythonVelocity.rows())
            {
                updatedNodePythonVelocity.row(rowIndex) = Eigen::RowVectorXd::Zero(this->NdofXnode);
            }
        }
        return updatedNodePythonVelocity;
    }

    template <typename DislocationNetworkType>
    Eigen::MatrixXd PyGlideSolver<DislocationNetworkType>::updateUnconnectedNodePythonVelocity(const Eigen::MatrixXd& nodePythonVelocity, const std::map<size_t, size_t>& removedNodeIDmap) const
    {
        std::set<size_t> removedIndices;
        for (const auto& pair : removedNodeIDmap)
        {
            removedIndices.insert(pair.second);
        }
        Eigen::MatrixXd updatedNodePythonVelocity = Eigen::MatrixXd::Zero(nodePythonVelocity.rows() + removedIndices.size(), this->NdofXnode);
        
        int oldRow = 0;
        int newRow = 0;
        for (newRow = 0; oldRow < nodePythonVelocity.rows(); ++newRow)
        {
            if (removedIndices.find(newRow) != removedIndices.end())
            {
                updatedNodePythonVelocity.row(newRow) = Eigen::RowVectorXd::Zero(this->NdofXnode); // Insert zero row
            }
            else
            {
                updatedNodePythonVelocity.row(newRow) = nodePythonVelocity.row(oldRow); // Copy existing row
                ++oldRow;
            }
        }
        return updatedNodePythonVelocity;
    }


    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        const auto t0= std::chrono::system_clock::now();
        std::cout<<"PyGlideSolver solving " << "\n" << std::flush;
        
//        auto configIO(this->DN.io().configIO()); //evl files
//        const auto auxIO(this->DN.io().auxIO()); //DDaux files
//        const auto Temperature(this->DN.ddBase.poly.T);
//        configIO.finalize();
//        const auto segmentMap(configIO.segments());
//        
        DDconfigIO<dim> configIO(".");
        DDauxIO<dim> auxIO(".");
        std::ofstream f_file;
        std::ofstream F_labels;
        this->DN.output(configIO,auxIO,f_file,F_labels);
        const auto Temperature(this->DN.ddBase.poly.T);
        configIO.finalize();
        const auto segmentMap(configIO.segments());

        std::map<size_t, DislocationNodeIO<3>> nodeMap;
        std::map<size_t, size_t> boundaryNodeMap;
        int k(0);
        for (const auto &node : configIO.nodes())
        {
            nodeMap.emplace(node.sID, node);
            if (node.meshLocation==1)
            {
                boundaryNodeMap.emplace(node.sID, k);
            }
            k++;
        }

        std::map<std::pair<size_t, size_t>, std::map<size_t, DislocationQuadraturePointIO<3>>> qPointsMap;
        std::set<std::pair<size_t, size_t>> segmentConnetions;
        int qpSize(0);
        for (const auto &qp : auxIO.quadraturePoints())
        {
            const std::pair<size_t, size_t> outerKey(std::make_pair(qp.sourceID, qp.sinkID));
            const size_t innerKey(qp.qID);
            qPointsMap[outerKey].emplace(innerKey, qp);
            segmentConnetions.insert(outerKey);
            qpSize+=1;
        }
        int segmentsSize(segmentConnetions.size());
        
        std::map<size_t, size_t> nodesToRemove; // Unconnected nodes
        size_t j(0);
        for (const auto& node : nodeMap)
        {
            bool found = false;
            for (const auto& qp : qPointsMap) 
            {
                if (qp.first.first == node.first || qp.first.second == node.first)
                {
                    found = true;
                    break;
                }
            }
            if (!found) 
            {
                nodesToRemove.emplace(node.first, j);
            }
            j++;
        }
        for (const auto& nodeKey : nodesToRemove)
        {
            nodeMap.erase(nodeKey.first);
        }
        int nodeSize(nodeMap.size());
        
        
        Eigen::MatrixXd nodalPositions(Eigen::MatrixXd::Zero(nodeSize+qpSize,this->DN.NdofXnode));
        int positionsCount(0);
        for (const auto &node : nodeMap)
        {
            nodalPositions.row(positionsCount) << node.second.P.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
            positionsCount+=1;
        }
        for (const auto &outerPair : qPointsMap)
        {
            for (const auto &innerPair : outerPair.second)
            {
                const auto &qp(innerPair.second);
                nodalPositions.row(positionsCount) << qp.r.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
                positionsCount+=1;
            }
        }
        

        Eigen::MatrixXd nodalStresses(Eigen::MatrixXd::Zero(nodeSize+qpSize,this->DN.NdofXnode*this->DN.NdofXnode));
        int stressesCount(nodeMap.size());
        for (const auto &outerPair : qPointsMap)
        {
            for (const auto &innerPair : outerPair.second)
            {
                const auto &qp(innerPair.second);
                Eigen::Matrix<double, 3, 3> tempStress(Eigen::Matrix<double,3,3>::Zero());
                tempStress <<  qp.stress*this->DN.ddBase.poly.mu_SI*1e-9;
                nodalStresses.row(stressesCount) << tempStress.reshaped<Eigen::RowMajor>().transpose();
                stressesCount+=1;
            }
        }
        
        
        Eigen::MatrixXd nodalBurgersVector(Eigen::MatrixXd::Zero(segmentsSize,this->DN.NdofXnode));
        int burgersVectorCount(0);
        for (const auto &outerPair : qPointsMap)
        {
            nodalBurgersVector.row(burgersVectorCount) << segmentMap.at(outerPair.first).b.transpose()*this->DN.ddBase.poly.b_SI*1.0e10;
            burgersVectorCount+=1;
        }
        
        
        std::list<std::list<int>> segmentConnections;
        std::list<int> singleConnection;
        size_t qPointsCounter(1);
        for (const auto &outerPair : qPointsMap)
        {
            const auto &sourceNodeID(outerPair.first.first);
            const auto &sinkNodeID(outerPair.first.second);
            singleConnection.push_back(sourceNodeID);
            for (size_t k = 0; k < outerPair.second.size(); ++k)
            {
                singleConnection.push_back(qPointsCounter + k + nodeMap.rbegin()->first);
            }
            qPointsCounter += outerPair.second.size();
            singleConnection.push_back(sinkNodeID);
            segmentConnections.push_back(singleConnection);
            singleConnection.erase(singleConnection.begin(),singleConnection.end());
        }
        sequenceOuterConnections(segmentConnections);
        adjustQPConnections(segmentConnections);
        adjustMinValueToZero(segmentConnections);
        //printNestedList(segmentConnections);
        
        
        Eigen::MatrixXd nodePythonVelocity;
        if (nodeSize>0)
        {
            std::cout<< "Computing Nodal Velocities from ML..." << std::endl;
            nodePythonVelocity = pyModule.attr("PIGNN")().attr("MLmobility")(nodalPositions,nodalStresses,nodalBurgersVector,segmentConnections,Temperature).template cast<Eigen::MatrixXd>();
        }
        else
        {
            std::cout<< "There are No Remaining Dislocations..." << std::endl;
            nodePythonVelocity = Eigen::MatrixXd::Zero(nodeMap.size(), 3); // Returns an empty matrix if nodeSize is 0

        }

        Eigen::MatrixXd updatedVelo = updateUnconnectedNodePythonVelocity(nodePythonVelocity, nodesToRemove);
        Eigen::MatrixXd updatedVelocity = updateBoundaryNodePythonVelocity(updatedVelo, boundaryNodeMap);
        updatedVelocity /= 3.228757859764734e01; //Convert to Velocities
        
//        pybind11::print(pyModule.attr("PIGNN")().attr("MLmobility")(positionNodes,stressNodes,burgersNodes,connections,Temperature).template cast<Eigen::MatrixXd>());
//        const std::string directory = "/Users/matthewmaron/Documents/MoDELib2-ML/tutorials/DislocationDynamics/periodicDomains/uniformLoadControllerML_500K/evl/";
//        const std::string base = "MLData_";
//        const std::string extension = ".txt";
//        const double runid(this->DN.ddBase.simulationParameters.runID);
//        std::ostringstream oss;
//        oss << directory << base << runid << extension;
//        std::cout << "Writing " << oss.str() << std::endl;
//        std::ofstream myfile;
//        myfile.open(oss.str());
//          myfile << "Temperature:" << " \n" << Temperature << " \n"<< " \n"<<std::endl;
//          myfile << "Positions Nodes + QP:" << " \n" << nodalPositions << " \n"<< " \n"<<std::endl;
//          myfile << "Stress Nodes + QP:" << " \n" << std::setprecision(15) << nodalStresses<< " \n"<< " \n"<<std::endl;
//          myfile << "Burgers Vectors Nodes:" << " \n" << nodalBurgersVector<< " \n"<< " \n"<<std::endl;
//          myfile << "Cell Connections:" << " \n" << std::endl;
//          printNestedListToFile(segmentConnections,myfile);
//          myfile << " \n" << " \n" << std::endl;
//          myfile << "Output Velocity:" << " \n" << nodePythonVelocity << " \n"<< " \n"<<std::endl;
//          myfile.close();
//        
        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
       
        return updatedVelocity.reshaped<Eigen::RowMajor>();
        
    }

    #else // COMPLIED WITHOUT PYBIND11

    template <typename DislocationNetworkType>
    PyGlideSolver<DislocationNetworkType>::PyGlideSolver(const DislocationNetworkType& DN_in) :
    /* init */ DislocationGlideSolverBase<DislocationNetworkType>(DN_in)
    {
        std::cout<<greenBoldColor<<"Creating PyGlideSolver"<<defaultColor<<std::endl;
        throw std::runtime_error("PyGlideSolver must be compiled with pybind11");
    }

    template <typename DislocationNetworkType>
    Eigen::VectorXd PyGlideSolver<DislocationNetworkType>::getNodeVelocities() const
    {
        return Eigen::VectorXd();
    }
    #endif


    template class PyGlideSolver<DislocationNetwork<3,0>>;
    
}
#endif
