//
// Created by Jelmer Mulder on 13/04/15.
//

#include "input/parser.h"
#include "input/csvparser.h"
#include "input/gmlparser.h"
#include "input/graphmlparser.h"
#include "input/ledaparser.h"
#include "input/lgfparser.h"
#include "input/stringparser.h"
#include "input/edgelistparser.h"
#include "input/bpparser.h"
#include "input/bpblastparser.h"
#include "input/bpcandlistparser.h"
#include "input/bplgfparser.h"
#include "lna.h"
#include "molecule.h"
#include "product.h"
#include "bronkerboschconnected.h"

#ifndef LANA_H_
#define LANA_H_

namespace nina {
namespace gna {

template<typename GR, typename BGR>
class Lana : public PairwiseLocalNetworkAlignment<GR, BGR>
{

public:

    /// The graph type of the input graphs
    typedef GR Graph;
    /// The graph type of the bipartite matching graph
    typedef BGR BpGraph;
    /// Base class type
    typedef typename nina::gna::Lana<GR, BGR>::Graph::Node GraphNodeType;
    typedef typename BpGraph::Node BpNode;
    typedef typename BpGraph::Edge BpEdge;
    typedef typename BpGraph::NodeIt BpNodeIt;
    typedef typename BpGraph::EdgeIt BpEdgeIt;

    typedef PairwiseLocalNetworkAlignment<GR, BGR> Parent;
    /// Type of a map assigning a boolean to every matching edge
    typedef typename BpGraph::template EdgeMap<bool> BpBoolMap;
    /// Type of a matching map: maps a node to its matching edge
    typedef typename BpGraph::template NodeMap<BpEdge> BpMatchingMap;

    typedef Product<Graph, BpGraph> ProductType;

    /// Type of input graph parser
    typedef Parser<Graph> ParserType;
    typedef CSVParser<Graph> ParserCsvType;
    typedef GMLParser<Graph> ParserGmlType;
    typedef GraphMLParser<Graph> ParserGraphMLType;
    typedef LedaParser<Graph> ParserLedaType;
    typedef LgfParser<Graph> ParserLgfType;
    typedef StringParser<Graph> ParserStringType;
    typedef EdgeListParser<Graph> ParserEdgeListType;

    /// Type of input matching graph parser
    typedef BpParser<Graph, BpGraph> BpParserType;
    typedef BpLgfParser<Graph, BpGraph> BpParserLgfType;
    typedef BpCandListParser<Graph, BpGraph> BpParserCandListType;
    typedef BpBlastParser<Graph, BpGraph> BpParserBlastType;

    using Parent::_pMatchingGraph;

    /// Output format type
    typedef enum {
        BP_OUT_DOT,
        BP_OUT_GML,
        BP_OUT_LGF,
        BP_OUT_SIF,
        BP_OUT_JSON,
        BP_OUT_NEATO,
        BP_OUT_CSV_MATCHED,
        BP_OUT_CSV_UNMATCHED_IN_G1,
        BP_OUT_CSV_UNMATCHED_IN_G2,
        BP_OUT_CSV_ALIGNMENT,
        BP_OUT_SIF_EDGE_ATTR,
        BP_OUT_SIF_NODE_ATTR,
    } OutputFormatEnum;

    /// Matching graph input format type
    typedef enum {
        BP_IN_CAND_LIST,
        BP_IN_BLAST,
        BP_IN_LGF,
    } BpInputFormatEnum;

    /// Input format type
    typedef enum {
        IN_GML,
        IN_GRAPHML,
        IN_STRING,
        IN_LGF,
        IN_CSV,
        IN_LEDA,
        IN_EDGE_LIST,
    } InputFormatEnum;


    /// Initializes the graphs using the specified parsers
    virtual bool init(ParserType* pParserG1,
                      ParserType* pParserG2,
                      BpParserType* pParserGm);


    static ParserType* createParser(const std::string& filename,
                                    InputFormatEnum fmt);


    static BpParserType* createBpParser(const std::string& filename,
                                        BpInputFormatEnum fmt,
                                        ParserType* pParserG1,
                                        ParserType* pParserG2);

    virtual int solve();

    virtual void getSolution(BpBoolMap &m, int i) const;

    virtual void getSolution(BpMatchingMap &m, int i) const;

    virtual int getNumberOfSolutions() const;

private:
    TEMPLATE_GRAPH_TYPEDEFS(Graph);
};




template<typename GR, typename BGR>
int Lana<GR, BGR>::solve() {

    // TODO: Properly initialize and pass these variables elsewhere.
    bool noAuto = false;

    Molecule<Graph> m1(_pMatchingGraph->getG1(), _pMatchingGraph->getMapLabelG1());
    Molecule<Graph> m2(_pMatchingGraph->getG2(), _pMatchingGraph->getMapLabelG2());

    // TODO: Use proper shell or remove it.
    ProductType prod(m1, m2, *_pMatchingGraph, 0);
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
        std::cerr << "Product graph has " << prod.getNumNodes()
        << " nodes and " << prod.getNumEdges()
        << " edges" << std::endl;
    }

    BronKerboschConnected<Graph, BpGraph> bk(prod);
    lemon::Timer t;
    bk.run(BronKerbosch<Graph>::BK_CLASSIC);
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
        std::cerr << "Time: " << t.realTime() << "s" << std::endl;
        std::cerr << "#max-cliques: " << bk.getNumberOfMaxCliques() << std::endl;
    }

    const std::vector< std::vector<GraphNodeType> > cliques = noAuto ? bk.getMaxCliques() : prod.removeAutomorphisms(bk.getMaxCliques());

    prod.printDOT(std::cerr);
    for (size_t i = 0; i < cliques.size(); ++i)
    {
        prod.printProductNodeVector(cliques[i], std::cerr);
    }


    return 0;
}

template<typename GR, typename BGR>
void Lana<GR, BGR>::getSolution(BpBoolMap &m, int i) const {

}

template<typename GR, typename BGR>
void Lana<GR, BGR>::getSolution(BpMatchingMap &m, int i) const {

}

template<typename GR, typename BGR>
int Lana<GR, BGR>::getNumberOfSolutions() const {
    return 0;
}


template<typename GR, typename BGR>
inline bool Lana<GR, BGR>::init(ParserType* pParserG1,
                                   ParserType* pParserG2,
                                   BpParserType* pParserGm)
{
    bool res = false;
    if (pParserG1 && pParserG2 && pParserGm)
        res = _pMatchingGraph->init(pParserG1, pParserG2, pParserGm);
    return res;
}


template<typename GR, typename BGR>
inline typename Lana<GR, BGR>::ParserType*
Lana<GR, BGR>::createParser(const std::string& filename,
                               InputFormatEnum fmt)
{
    ParserType* pParser = NULL;
    switch (fmt)
    {
        case IN_GML:
            pParser = new ParserGmlType(filename);
            break;
        case IN_GRAPHML:
            pParser = new ParserGraphMLType(filename);
            break;
        case IN_STRING:
            // TODO: Correct threshold?
            pParser = new ParserStringType(filename, 0);
            break;
        case IN_LGF:
            pParser = new ParserLgfType(filename);
            break;
        case IN_CSV:
            // TODO: Correct threshold?
            pParser = new ParserCsvType(filename, 0);
            break;
        case IN_LEDA:
            pParser = new ParserLedaType(filename);
            break;
        case IN_EDGE_LIST:
            pParser = new ParserEdgeListType(filename);
            break;
    }

    return pParser;
}

template<typename GR, typename BGR>
inline typename Lana<GR, BGR>::BpParserType*
Lana<GR, BGR>::createBpParser(const std::string& filename,
                                 BpInputFormatEnum fmt,
                                 ParserType* pParserG1,
                                 ParserType* pParserG2)
{
    BpParserType* pBpParser = NULL;
    switch (fmt)
    {
        case BP_IN_CAND_LIST:
            pBpParser = new BpParserCandListType(filename, pParserG1, pParserG2);
            break;
        case BP_IN_BLAST:
            // TODO: Correct threshold?
            pBpParser = new BpParserBlastType(filename, pParserG1, pParserG2, 0);
            break;
        case BP_IN_LGF:
            pBpParser = new BpParserLgfType(filename, pParserG1, pParserG2);
            break;
    }

    return pBpParser;
}


} // namespace gna
} // namespace nina


#endif //LANA_H_
