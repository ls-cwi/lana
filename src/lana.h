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
#include "proteinnetwork.h"
#include "product.h"
#include "bronkerboschconnected.h"
#include "bronkerboschconnectedrelaxed.h"
#include "options.h"
#include <fstream>

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

    typedef ProductGraph<Graph, BpGraph> ProductType;
    typedef std::vector<typename Graph::Node> NodeVector;
    typedef typename NodeVector::const_iterator NodeVectorIt;
    typedef std::vector<NodeVector> NodeVectorVector;
    typedef typename std::vector<NodeVector>::const_iterator NodeVectorVectorIt;
    typedef typename BpGraph::RedNodeIt BpRedNodeIt;
    typedef typename Graph::template NodeMap<bool> NodeFilterMap;
    typedef typename lemon::FilterNodes<Graph, typename Graph::template NodeMap<bool> > NodeFilterGraph;
    typedef BronKerbosch<NodeFilterGraph, BpGraph, Graph> BronKerboschType;
    typedef BronKerboschConnected<NodeFilterGraph, BpGraph, Graph> BronKerboschConnectedType;
    typedef BronKerboschConnectedRelaxed<NodeFilterGraph, BpGraph, Graph> BronKerboschConnectedRelaxedType;

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


    Lana(const Options& options)
        : Parent()
        , _options(options)
    {
    }

    ~Lana()
    {
        delete _prod;
    }
protected:
    NodeVectorVector _solutions;
    ProductType* _prod;
    Options _options;


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
                                        ParserType* pParserG2,
                                        double evalCutOff);

    virtual int solve();

    virtual void getSolution(BpBoolMap &m, int i) const;

    virtual void getSolution(BpMatchingMap &m, int i) const;

    virtual int getNumberOfSolutions() const;

    void parseOutputString(std::string const &str);

    void printCliqueSizeFrequencies(std::ofstream &freq_file);

    void addOutput(OutputFormatEnum fmt);

    void checkDuplicates(NodeVectorVector vector);
private:
    TEMPLATE_GRAPH_TYPEDEFS(Graph);

};



template<typename GR, typename BGR>
int Lana<GR, BGR>::solve() {

    // Create two protein networks.
    ProteinNetwork<Graph> n1(_pMatchingGraph->getG1(), _pMatchingGraph->getMapLabelG1());
    ProteinNetwork<Graph> n2(_pMatchingGraph->getG2(), _pMatchingGraph->getMapLabelG2());

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
        std::cout << "Generating ProductGraph." << std::endl;
    }

    // Create the product graph
    _prod = new ProductType(n1, n2, *_pMatchingGraph, _options);

    // Prints the product graph in DOT format.
    if (!_options._prod_file_name->empty())
    {
        std::ofstream prod_file;
        prod_file.open(_options._prod_file_name->c_str());
        _prod->printDOT(prod_file);
        prod_file.close();
    }



    int numComponents = _prod->getNumComponents();

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
        std::cout << "ProductGraph graph has " << _prod->getNumNodes()
        << " nodes, " << _prod->getNumEdges()
        << " edges, and " << numComponents << " components." << std::endl;
    }

    // The product graph is split up into c-connected components.
    // Here we iterate over these components.
    for (int i=0; i<numComponents; i++)
    {
        int size = _prod->getComponentSize(i);

        if (i % 1000 == 0)
        {
            if (g_verbosity >= VERBOSE_ESSENTIAL)
            {
                std::cout << "Checking component " << i << "/" << numComponents << "." << std::endl;
            }
        }

        // Skip small componenents.
        // The vast majority of the componenets tend to be size 1, and are thus skipped.
        if (size <= std::max<int>(2, _options._minCliqueSize))
            continue;

        if (i % 1000 != 0)
        {
            if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
            {
                std::cout << "Checking component " << i << "/" << numComponents << " (size: " << size << ")" <<
                std::endl;
            }

        }

        _prod->enableComponent(i);

        if (g_verbosity >= VERBOSE_DEBUG)
        {
            std::cout << "Generating BronKerbosch object." << std::endl;
        }

        // If the mse parameter (max s-edges) is set to 0, we use BronKerboschConnected
        // instead of BronKerboschConnectedRelaxed because it's more efficient.
        BronKerboschConnectedType* bk;
        if (_options._nMaxSEdges > 0) {
            bk = new BronKerboschConnectedRelaxedType(*_prod, _options);
        }
        else
        {
            bk = new BronKerboschConnectedType(*_prod, _options);
        }


        if (g_verbosity >= VERBOSE_DEBUG)
        {
            std::cout << "Done generating BronKerbosch object." << std::endl;
            std::cout << "Running BronKerbosch algorithm." << std::endl;
        }

        lemon::Timer t;

        bk->run(BronKerboschType::BK_CLASSIC);

        if (g_verbosity >= VERBOSE_DEBUG)
        {
            std::cout << "Done running BronKerbosch algorithm." << std::endl;
        }

        if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
        {
            std::cout << "Time: " << t.realTime() << "s" << std::endl;
            std::cout << "Number of alignments found: " << bk->getNumberOfMaxCliques() << std::endl << std::endl;
        }

        NodeVectorVector x = _options._removeAutomorphisms ? _prod->removeAutomorphisms(bk->getMaxCliques()) : bk->getMaxCliques();
        _solutions.insert(_solutions.end(), x.begin(), x.end());
        _prod->disableComponent(i);

        delete bk;

    }
    checkDuplicates(_solutions);
    _prod->enableAllComponents();

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
        std::cout << "Algorithm has succesfully terminated. " << _solutions.size() << " alignments were found." << std::endl;
    }


    // Prints found alignements to STDOUT in human readable format if the -sol flag is given.
    if (_options._printProductVector)
    {
        for (size_t i=0; i<_solutions.size(); ++i) {
            _prod->printProductNodeVector(_solutions.at(i), std::cout);
        }
    }



    return 0;
}

template<typename GR, typename BGR>
void Lana<GR, BGR>::checkDuplicates(NodeVectorVector vec)
{
    if (g_verbosity >= VERBOSE_DEBUG)
    {
        std::cout << "Checking duplicates..." << std::endl;
    }
    unsigned long size = vec.size();
    std::sort( vec.begin(), vec.end() );
    vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
    unsigned long new_size = vec.size();

    if (new_size < size)
    {
        std::cerr << "Duplicates were detected (and removed). This should not occur!" << std::endl;
    }

    if (g_verbosity >= VERBOSE_DEBUG)
    {
        std::cout << "Old size: size: " << size << std::endl;
        std::cout << "New size: size: " << new_size << std::endl;
    }

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
        // TODO: Do these have to be destroyed?
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
                                 ParserType* pParserG2,
                                 double evalCutOff)
{
    BpParserType* pBpParser = NULL;
    switch (fmt)
    {
        case BP_IN_CAND_LIST:
            pBpParser = new BpParserCandListType(filename, pParserG1, pParserG2);
            break;
        case BP_IN_BLAST:
            pBpParser = new BpParserBlastType(filename, pParserG1, pParserG2, evalCutOff);
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


template<typename GR, typename BGR>
inline void Lana<GR, BGR>::printCliqueSizeFrequencies(std::ofstream &freq_file)
{
    int n_solutions = _solutions.size();
    int solution_sizes[n_solutions];
    int max_size = 0;

    for (int i = 0; i < n_solutions; i++)
    {
        solution_sizes[i] = _solutions[i].size();
        if (solution_sizes[i] > max_size)
        {
            max_size = solution_sizes[i];
        }
    }

    int size_frequencies[max_size + 1];
    memset(size_frequencies, 0, sizeof size_frequencies);

    for (int j = 0; j < n_solutions; j++)
    {
        size_frequencies[solution_sizes[j]] = size_frequencies[solution_sizes[j]] + 1;
    }

    for (int k = 0; k <= max_size; k++)
    {
        freq_file << k << ", " << size_frequencies[k] << std::endl;
    }
    freq_file.close();
}

} // namespace gna
} // namespace nina


#endif //LANA_H_
