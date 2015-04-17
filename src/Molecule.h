//
// Created by Jelmer Mulder on 16/04/15.
//

#ifndef LANA_MOLECULE_H
#define LANA_MOLECULE_H

#include <lemon/core.h>

template <typename GR>
class Molecule {

public:
    typedef GR Graph;


protected:
    TEMPLATE_GRAPH_TYPEDEFS(Graph);



public:
    /// Labels of the nodes
    typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;

    Molecule(const Graph& g, const OrigLabelNodeMap& labelMap)
        : _g(g)
        , _labelMap(labelMap)
    {
    }


    const std::string& getLabel(Node n) const
    {
        assert(n != lemon::INVALID && _g.valid(n));
        return _labelMap[n];
    }


    const std::string& getLabel2(Node n) const
    {
        assert(n != lemon::INVALID && _g.valid(n));
        return _labelMap[n];
    }

    const Graph& getGraph() const
    {
        return _g;
    }

    int getAtomType(Node n) const
    {
        return 1;
    }


protected:
    const Graph& _g;

    const OrigLabelNodeMap& _labelMap;

};

#endif //LANA_MOLECULE_H
