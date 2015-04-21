//
// Created by Jelmer Mulder on 13/04/15.
//

#ifndef LNA_H_
#define LNA_H_

#include <string>
#include <lemon/time_measure.h>
#include "input/matchinggraph.h"
#include "output/output.h"

namespace nina {
namespace gna {

template<typename GR, typename BGR>
class PairwiseLocalNetworkAlignment {

public:
    /// The graph type of the input graphs
    typedef GR Graph;
    /// The graph type of the bipartite matching graph
    typedef BGR BpGraph;

    typedef typename BpGraph::Node BpNode;
    typedef typename BpGraph::Edge BpEdge;
    typedef typename BpGraph::NodeIt BpNodeIt;
    typedef typename BpGraph::EdgeIt BpEdgeIt;
    typedef typename BpGraph::IncEdgeIt BpIncEdgeIt;
    typedef typename BpGraph::RedNode BpRedNode;
    typedef typename BpGraph::BlueNode BpBlueNode;
    typedef typename BpGraph::RedNodeIt BpRedNodeIt;
    typedef typename BpGraph::BlueNodeIt BpBlueNodeIt;

    /// Type of a map assigning a boolean to every matching edge
    typedef typename BpGraph::template EdgeMap<bool> BpBoolMap;
    /// Type of a matching map: maps a node to its matching edge
    typedef typename BpGraph::template NodeMap<BpEdge> BpMatchingMap;
    /// Type of the matching graph
    typedef MatchingGraph<Graph, BpGraph> MatchingGraphType;
    /// Type of the output
    typedef Output<Graph, BpGraph> OutputType;
    /// Type of a vector of output methods
    typedef std::vector<OutputType*> OutputVector;
    /// Iterator type of output vector
    typedef typename OutputVector::const_iterator OutputVectorIt;

protected:
    MatchingGraphType* _pMatchingGraph;
    OutputVector _outputs;

public:
    PairwiseLocalNetworkAlignment();

    virtual int solve() = 0;

    virtual void getSolution(BpBoolMap& m, int i = 0) const = 0;

    virtual void getSolution(BpMatchingMap& m, int i = 0) const = 0;

    virtual int getNumberOfSolutions() const = 0;

    virtual void generateOutput(const typename OutputType::OutputType outType = OutputType::MINIMAL,
                                const std::string& filename = std::string());

    void addOutput(OutputType* pOutput)
    {
        std::cout << "Adding an output type..." << std::endl;
        _outputs.push_back(pOutput);
    }
};

template<typename GR, typename BGR>
inline PairwiseLocalNetworkAlignment<GR, BGR>::
PairwiseLocalNetworkAlignment()
        : _pMatchingGraph(new MatchingGraphType())
        , _outputs()
{
}



template<typename GR, typename BGR>
inline void PairwiseLocalNetworkAlignment<GR, BGR>::
generateOutput(const typename OutputType::OutputType outType,
               const std::string& filename)
{

    int n = getNumberOfSolutions();
    for (int i = 0; i < n; i++)
    {
        BpMatchingMap matchingMap(_pMatchingGraph->getGm(), lemon::INVALID);
        getSolution(matchingMap, i);

        std::cout << "Writing solution " << i << "." << std::endl;
        for (OutputVectorIt it = _outputs.begin(); it != _outputs.end(); it++)
        {
            std::cout << "    - output type: " << (*it)->getExtension() << std::endl;
            if (filename.empty() || filename == "-")
            {
                // output to std::out
                (*it)->write(matchingMap, outType, std::cout);
            }
            else
            {
                // this is safe, log_10(MAXINT) < 128
                char buf[128];
                sprintf(buf, "-%d", i+1);
                std::string number(buf);

                std::string newFilename = filename + number;
                (*it)->write(matchingMap, outType, newFilename);
            }
        }
    }
}

} // namespace gna
} // namespace nina

#endif //LANA_H_
