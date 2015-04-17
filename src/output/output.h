//
// Created by Jelmer Mulder on 16/04/15.
//

#ifndef LANA_OUTPUT_H
#define LANA_OUTPUT_H

namespace nina {
namespace gna {

template<typename GR, typename BGR>
class Output {

public:
    /// The graph type of the input graphs
    typedef GR Graph;
    /// The graph type of the bipartite matching graph
    typedef BGR BpGraph;

    typedef typename BpGraph::Node BpNode;
    typedef typename BpGraph::Edge BpEdge;
    /// What to output
    typedef enum {
            MINIMAL,    // nodes and matching edges
            ORIG_EDGES, // nodes, matching and original edges
            FULL,       // all nodes, matching edges and all original edges
    } OutputType;
protected:
    typedef MatchingGraph<Graph, BpGraph> MatchingGraphType;
    typedef typename BpGraph::template NodeMap<BpEdge> BpMatchingMapType;

    const MatchingGraphType&  _matchingGraph;


public:
    Output(const MatchingGraphType& matchingGraph)
            : _matchingGraph(matchingGraph)
    {
    }

    virtual void write(const BpMatchingMapType& matchingMap,
                       OutputType outputType,
                       std::ostream& outFile) const = 0;

    virtual void write(const BpMatchingMapType& matchingMap,
                       OutputType outputType,
                       const std::string& filename) const;

    virtual std::string getExtension() const = 0;
};

template<typename GR, typename BGR>
inline void Output<GR, BGR>::write(const BpMatchingMapType& matchingMap,
                                   OutputType outputType,
                                   const std::string& filename) const
{
    std::ofstream out((filename + getExtension()).c_str());
    write(matchingMap, outputType, out);
}
}
}
#endif //LANA_OUTPUT_H
