//
// Created by Jelmer Mulder on 13/04/15.
//

#include <output/outputsif.h>
#include <output/outputlgf.h>
#include <output/outputdot.h>
#include <output/outputgml.h>
#include <output/outputneato.h>
#include <output/outputcsv.h>
#include <output/outputeda.h>
#include <output/outputnoa.h>
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
    typedef std::vector<typename Graph::Node> NodeVector;
    typedef typename NodeVector::const_iterator NodeVectorIt;
    typedef std::vector<NodeVector> NodeVectorVector;
    typedef typename std::vector<NodeVector>::const_iterator NodeVectorVectorIt;
    typedef typename BpGraph::RedNodeIt BpRedNodeIt;

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

    /// Type of the output
    typedef typename Parent::OutputType OutputType;
    typedef OutputSif<Graph, BpGraph> OutputSifType;
    typedef OutputLgf<Graph, BpGraph> OutputLgfType;
    typedef OutputDot<Graph, BpGraph> OutputDotType;
    typedef OutputGml<Graph, BpGraph> OutputGmlType;
    typedef OutputNeato<Graph, BpGraph> OutputNeatoType;
    typedef OutputCsv<Graph, BpGraph> OutputCsvType;
    typedef OutputEda<Graph, BpGraph> OutputEdaType;
    typedef OutputNoa<Graph, BpGraph> OutputNoaType;

    using Parent::_pMatchingGraph;
    using Parent::_outputs;
    using Parent::addOutput;

    /// Output format type
    typedef enum {
        BP_OUT_DOT,
        BP_OUT_GML,
        BP_OUT_LGF,
        BP_OUT_SIF,
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

protected:
    NodeVectorVector _solutions;
    ProductType* _prod;


public:


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

    void parseOutputString(std::string const &str);

    void addOutput(OutputFormatEnum fmt);

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
    // TODO: Destroy?
    _prod = new ProductType(m1, m2, *_pMatchingGraph, 0);


    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
        std::cerr << "Product graph has " << _prod->getNumNodes()
        << " nodes and " << _prod->getNumEdges()
        << " edges" << std::endl;
    }

    BronKerboschConnected<Graph, BpGraph> bk(*_prod);
    lemon::Timer t;
    bk.run(BronKerbosch<Graph>::BK_CLASSIC);
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
        std::cerr << "Time: " << t.realTime() << "s" << std::endl;
        std::cerr << "#max-cliques: " << bk.getNumberOfMaxCliques() << std::endl;
    }

    // TODO: This makes a copy, a better way?
    _solutions = noAuto ? bk.getMaxCliques() : _prod->removeAutomorphisms(bk.getMaxCliques());



    // TODO: Remove or add argument for this.
//    _prod->printDOT(std::cout);

    for (size_t i = 0; i < _solutions.size(); ++i) {
        _prod->printProductNodeVector(_solutions.at(i), std::cout);
    }



    return 0;
}

template<typename GR, typename BGR>
void Lana<GR, BGR>::getSolution(BpBoolMap &m, int i) const {
    // TODO: Implement
    assert (false && "Not yet implemented");
}

template<typename GR, typename BGR>
void Lana<GR, BGR>::getSolution(BpMatchingMap &m, int i) const {
    NodeVector solution = _solutions.at(i);
    lemon::ArcLookUp<BpGraph> arcLookUpGm(_pMatchingGraph->getGm());

    for (NodeVectorIt it = solution.begin(); it != solution.end(); ++it)
    {
        // A node from the product graph.
        Node n = *it;

        // Nodes from the original G1 and G2
        Node n_g1 = _prod->getNodeG1(n);
        Node n_g2 = _prod->getNodeG2(n);

        // Nodes in the bipartite graph Gm
        BpNode n1_gm = _pMatchingGraph->mapG1ToGm(n_g1);
        BpNode n2_gm = _pMatchingGraph->mapG2ToGm(n_g2);

        // Add the edge from the bipartite graph Gm
        m[n1_gm] = m[n2_gm] = arcLookUpGm(n1_gm, n2_gm);
    }
}

template<typename GR, typename BGR>
int Lana<GR, BGR>::getNumberOfSolutions() const {
    return _solutions.size();
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

template<typename GR, typename BGR>
inline void Lana<GR, BGR>::parseOutputString(const std::string& str)
{
    size_t idx = 0, new_idx = 0;
    while ((new_idx = str.find(",", idx)) != std::string::npos)
    {
        std::string substr = str.substr(idx, new_idx - idx);
        if (substr == "0")
            addOutput(static_cast<OutputFormatEnum>(0));
        else if (substr == "1")
            addOutput(static_cast<OutputFormatEnum>(1));
        else if (substr == "2")
            addOutput(static_cast<OutputFormatEnum>(2));
        else if (substr == "3")
            addOutput(static_cast<OutputFormatEnum>(3));
        else if (substr == "4")
            addOutput(static_cast<OutputFormatEnum>(4));
        else if (substr == "5")
            addOutput(static_cast<OutputFormatEnum>(5));
        else if (substr == "6")
            addOutput(static_cast<OutputFormatEnum>(6));
        else if (substr == "7")
            addOutput(static_cast<OutputFormatEnum>(7));
        else if (substr == "8")
            addOutput(static_cast<OutputFormatEnum>(8));
        else if (substr == "9")
            addOutput(static_cast<OutputFormatEnum>(9));
        else if (substr == "10")
            addOutput(static_cast<OutputFormatEnum>(10));
        else if (substr == "11")
            addOutput(static_cast<OutputFormatEnum>(11));

        idx = new_idx + 1;
    }
}
template<typename GR, typename BGR>
inline void Lana<GR, BGR>::addOutput(OutputFormatEnum fmt)
{
    OutputType* pOutput = NULL;

    switch (fmt)
    {
        case BP_OUT_SIF:
            pOutput = new OutputSifType(*_pMatchingGraph);
            break;
        case BP_OUT_DOT:
            pOutput = new OutputDotType(*_pMatchingGraph);
            break;
        case BP_OUT_GML:
            pOutput = new OutputGmlType(*_pMatchingGraph);
            break;
        case BP_OUT_LGF:
            pOutput = new OutputLgfType(*_pMatchingGraph);
            break;
        case BP_OUT_NEATO:
            pOutput = new OutputNeatoType(*_pMatchingGraph);
            break;
        case BP_OUT_CSV_MATCHED:
            pOutput = new OutputCsvType(*_pMatchingGraph, OutputCsvType::CSV_MATCHED);
            break;
        case BP_OUT_CSV_UNMATCHED_IN_G1:
            pOutput = new OutputCsvType(*_pMatchingGraph, OutputCsvType::CSV_UNMATCHED_IN_G1);
            break;
        case BP_OUT_CSV_UNMATCHED_IN_G2:
            pOutput = new OutputCsvType(*_pMatchingGraph, OutputCsvType::CSV_UNMATCHED_IN_G2);
            break;
        case BP_OUT_CSV_ALIGNMENT:
            pOutput = new OutputCsvType(*_pMatchingGraph, OutputCsvType::CSV_ALIGNMENT);
            break;
        case BP_OUT_SIF_EDGE_ATTR:
            pOutput = new OutputEdaType(*_pMatchingGraph);
            break;
        case BP_OUT_SIF_NODE_ATTR:
            pOutput = new OutputNoaType(*_pMatchingGraph);
            break;
            // TODO: Implement and test (?) all outputs.
            assert (false && "Not yet implemented");
    }

    addOutput(pOutput);
}

} // namespace gna
} // namespace nina


#endif //LANA_H_
