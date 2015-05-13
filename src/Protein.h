//
// Created by Jelmer Mulder on 16/04/15.
//

#ifndef LANA_PROTEIN_H
#define LANA_PROTEIN_H

#include <lemon/core.h>

// TODO: Try to remove this class entirely.
template <typename GR>
class Protein
{
public:
    typedef GR Graph;


protected:
    TEMPLATE_GRAPH_TYPEDEFS(Graph);



public:
    /// Labels of the nodes
    typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;

    Protein(const Graph& g, const OrigLabelNodeMap& labelMap)
        : _g(g)
        , _labelMap(labelMap)
    {
    }


    const std::string& getLabel(Node n) const
    {
        assert(n != lemon::INVALID && _g.valid(n));
        return _labelMap[n];
    }



    const Graph& getGraph() const
    {
        return _g;
    }


protected:
    const Graph& _g;

    const OrigLabelNodeMap& _labelMap;

};

#endif //LANA_PROTEIN_H
