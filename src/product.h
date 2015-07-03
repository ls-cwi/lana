/*
 * product.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 */

// TODO:
// - Remove deg-1 nodes
// - Remove automorphisms

#ifndef PRODUCT_H
#define PRODUCT_H


#include <set>
#include <vector>
#include <lemon/core.h>
#include <assert.h>
#include <bits/stl_map.h>
#include <lemon/maps.h>
#include <lemon/smart_graph.h>
#include <input/matchinggraph.h>
#include "Protein.h"

namespace nina {
namespace gna {

template<typename GR, typename BGR>
class Product
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// The graph type of the bipartite matching graph
  typedef BGR BpGraph;
  /// Type of the matching graph
  typedef MatchingGraph<Graph, BpGraph> MatchingGraphType;

  typedef typename BpGraph::Node BpNode;
  typedef typename BpGraph::Edge BpEdge;
  typedef Protein<Graph> ProteinType;
  typedef std::vector<typename Graph::Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  typedef std::set<typename Graph::Node> NodeSet;
  typedef std::pair<NodeSet, NodeSet> NodeSetPair;
  typedef typename Graph::template NodeMap<bool> NodeFilterMap;
  typedef typename lemon::FilterNodes<Graph, typename Graph::template NodeMap<bool> > NodeFilterGraph;

  typedef enum {
    PRODUCT_NO_EDGE,
    PRODUCT_BLACK_EDGE,
    PRODUCT_RED_EDGE,
    PRODUCT_BLUE_EDGE
  } ProductEdgeType;

  typedef typename Graph::template EdgeMap<ProductEdgeType> EdgeTypeEdgeMap;

private:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Graph::template NodeMap<Node> NodeNodeMap;
  typedef std::multiset<int> IntSet;
  typedef typename Graph::template NodeMap<IntSet> IntSetNodeMap;
  typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;

public:
  Product(const ProteinType &prot1, const ProteinType &prot2, const MatchingGraphType& matchingGraph)
    : _prot1(prot1)
    , _prot2(prot2)
    , _matchingGraph(matchingGraph)
    , _g()
    , _nodeFilter(_g, false)
    , _fg(_g, _nodeFilter)
    , _prot1ToG(prot1.getGraph(), lemon::INVALID)
    , _prot2ToG(prot2.getGraph(), lemon::INVALID)
    , _gToProt1(_g)
    , _gToProt2(_g)
    , _edgeType(_g, PRODUCT_NO_EDGE)
    , _g1ToDeg1Neighbors(prot1.getGraph())
    , _g2ToDeg1Neighbors(prot2.getGraph())
    , _numNodes(0)
    , _numEdges(0)
  {
    generate();
  }

  ~Product()
  {
    delete _components;
  }







  NodeMatrix removeAutomorphisms(const NodeMatrix& fragments) const
  {
    NodeMatrix result;

    typedef std::map<NodeSetPair, NodeVector> NodeSetPairMap;
    typedef typename NodeSetPairMap::const_iterator NodeSetPairMapIt;

    NodeSetPairMap uniqueCommonFragments;
    for (NodeMatrixIt fragmentIt = fragments.begin(); fragmentIt != fragments.end(); ++fragmentIt)
    {
      NodeSet fragment1, fragment2;

      for (NodeVectorIt it = fragmentIt->begin(); it != fragmentIt->end(); ++it)
      {
        Node uv = *it;
        Node u = _gToProt1[uv];
        Node v = _gToProt2[uv];

        fragment1.insert(u);
        fragment2.insert(v);
      }

      NodeSetPair fragmentPair = std::make_pair(fragment1, fragment2);

      if (uniqueCommonFragments.find(fragmentPair) == uniqueCommonFragments.end())
      {
        uniqueCommonFragments[fragmentPair] = *fragmentIt;
        result.push_back(*fragmentIt);
        //std::cerr << "Added fragment with size: " << fragmentIt->size() << std::endl;
        //printProductNodeVector(*fragmentIt, std::cerr);
      }
      else
      {
        //std::cerr << "Skipped fragment with size: " << fragmentIt->size() << std::endl;
        //printProductNodeVector(*fragmentIt, std::cerr);
      }
    }

    return result;
  }

  int getNumNodes() const { return _numNodes; }
  int getNumEdges() const { return _numEdges; }

  const NodeFilterGraph & getGraph() const { return _fg; }

  ProductEdgeType connectivityEdge(Edge e) const { return _edgeType[e]; }

  Node getNodeG1(Node uv) { return _gToProt1[uv]; }

  Node getNodeG2(Node uv) { return _gToProt2[uv]; }

  int getNumComponents() { return _components->size(); }

  int getComponentSize(int index) { return (*_components)[index].size(); }

  void enableComponent(int index) { setComponentState(index, true); }

  void disableComponent(int index) { setComponentState(index, false); }

  void enableAllComponents() { lemon::mapFill(_g, _nodeFilter, true); }

  void disableAllComponents() { lemon::mapFill(_g, _nodeFilter, false); }

private:
  const ProteinType &_prot1;
  const ProteinType &_prot2;
  const MatchingGraphType& _matchingGraph;
  Graph _g;
  NodeMatrix* _components;
  NodeFilterMap _nodeFilter;
  NodeFilterGraph _fg;
  NodeNodeMap _prot1ToG;
  NodeNodeMap _prot2ToG;
  NodeNodeMap _gToProt1;
  NodeNodeMap _gToProt2;
  EdgeTypeEdgeMap _edgeType;
  NodeVectorMap _g1ToDeg1Neighbors;
  NodeVectorMap _g2ToDeg1Neighbors;

  int _numNodes;
  int _numEdges;

  void generate();

  void generate(const ProteinType &prot,
                const IntNodeMap& deg,
                IntSetNodeMap& degSet);

  void generateDeg1NeighborSet(const Graph& g,
                               const IntNodeMap& deg,
                               NodeVectorMap& deg1NeighborMap);

  void determineDegrees(const Graph& g, IntNodeMap& deg);

  void setComponentState(int index, bool state);

public:


  void printDOT(std::ostream& out) const;

  void printProductNodeJSON(Node uv, std::ostream& out) const
  {
    Node u = _gToProt1[uv];
    Node v = _gToProt2[uv];

    out << "            {" << std::endl
    << "              \"id1\": " << _prot1.getLabel(u) << "," << std::endl
    << "              \"id2\": " << _prot2.getLabel(v) << "," << std::endl
    << "            }";

    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

    assert(uNeighbors.size() == vNeighbors.size());
    for (size_t i = 0; i < uNeighbors.size(); ++i)
    {
      out << "," << std::endl;
      out << "            {" << std::endl
      << "              \"id1\": " << _prot1.getLabel(uNeighbors[i]) << "," << std::endl
      << "              \"id2\": " << _prot2.getLabel(vNeighbors[i]) << "," << std::endl
      << "            }";
    }
  }

  void printProductNode(Node uv, std::ostream& out) const
  {
    Node u = _gToProt1[uv];
    Node v = _gToProt2[uv];

    out << "[" << _g.id(uv) << ": " << _prot1.getLabel(u) << " -> "
    << _prot2.getLabel(v) << "]";

//    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
//    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

//    assert(uNeighbors.size() == vNeighbors.size());
//    for (size_t i = 0; i < uNeighbors.size(); ++i)
//    {
//      out << ", ";
//      out << "[" << _prot1.getLabel(uNeighbors[i])
//      << " (" << _prot1.getLabel(uNeighbors[i]) << ") , "
//      << _prot2.getLabel(vNeighbors[i])
//      << " (" << _prot2.getLabel(vNeighbors[i]) << ")]";
//    }
  }



  void printProductNodeVector(const NodeVector& nodes,
                              std::ostream& out) const
  {
    out << nodes.size() << ": ";
    bool first = true;
    for (NodeVectorIt it = nodes.begin(); it != nodes.end(); ++it)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << ", ";
      }
      printProductNode(*it, out);
    }
    out << std::endl;
  }

  void printProductNodeVectorJSON(const NodeVector& nodes,
                                  std::ostream& out,
                                  const bool fst) const
  {
    if (fst)
    {
      out << "        {" << std::endl
      << "          \"pairs\": [" << std::endl;
    }
    else
    {
      out << "," << std::endl
      << "        {" << std::endl
      << "          \"pairs\": [" << std::endl;
    }

    bool first = true;
    for (NodeVectorIt it = nodes.begin(); it != nodes.end(); ++it)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << "," << std::endl;
      }
      printProductNodeJSON(*it, out);
    }
    /* TODO: implement score function */
    out << std::endl
    << "          ]," << std::endl
    << "          \"score\": " << nodes.size() << std::endl
    << "        }";
  }

};

    template<typename GR, typename BGR>
inline void Product<GR,BGR>::determineDegrees(const Graph& g, IntNodeMap& deg)
{
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      ++deg[v];
    }
  }
}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::generateDeg1NeighborSet(const Graph& g,
                                                 const IntNodeMap& deg,
                                                 NodeVectorMap& deg1NeighborMap)
{
  for (NodeIt u(g); u != lemon::INVALID; ++u)
  {
    if (deg[u] > 1)
    {
      for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
      {
        Node v = g.oppositeNode(u, e);
        if (deg[v] == 1) {
          deg1NeighborMap[u].push_back(v);
        }
      }
    }
  }
}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::generate()
{

  const Graph& g1 = _prot1.getGraph();
  const Graph& g2 = _prot2.getGraph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);
  lemon::ArcLookUp<BpGraph> arcLookUpGm(_matchingGraph.getGm());

  IntSetNodeMap degSet1(g1);
  IntSetNodeMap degSet2(g2);

  IntNodeMap deg1(g1, 0);
  determineDegrees(g1, deg1);
  IntNodeMap deg2(g2, 0);
  determineDegrees(g2, deg2);


  // determine degrees
  generate(_prot1, deg1, degSet1);
  generate(_prot2, deg2, degSet2);


  // generate nodes
  for (NodeIt u(g1); u != lemon::INVALID; ++u)
  {
    for (NodeIt v(g2); v != lemon::INVALID; ++v)
    {

      BpNode bp_u = _matchingGraph.mapG1ToGm(u);
      BpNode bp_v = _matchingGraph.mapG2ToGm(v);
      BpEdge uv = arcLookUpGm(bp_u, bp_v);
      if (uv != lemon::INVALID )//   && degSet1[u] == degSet2[v])
      {
//        assert(deg1[u] == deg2[v]);
        Node uv = _g.addNode();
        _prot1ToG[u] = uv;
        _prot2ToG[v] = uv;
        _gToProt1[uv] = u;
        _gToProt2[uv] = v;
        ++_numNodes;
      }
    }
  }

  // generate deg1 sets
  generateDeg1NeighborSet(g1, deg1, _g1ToDeg1Neighbors);
  generateDeg1NeighborSet(g2, deg2, _g2ToDeg1Neighbors);

  // generate c-edges ('red')
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1)
  {
    Node u1 = _gToProt1[u1v1];
    Node v1 = _gToProt2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2)
    {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _gToProt1[u2v2];
      Node v2 = _gToProt2[u2v2];


      if (u1 != u2 && v1 != v2)
      {

        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 && v1v2)
        {
          _edgeType[_g.addEdge(u1v1, u2v2)] = PRODUCT_RED_EDGE;
          ++_numEdges;
        }
      }
    }
  }


  IntNodeMap component_labels(_g);
  int n_components = lemon::connectedComponents(_g, component_labels);
  _components = new NodeMatrix (n_components, std::vector<Node>(0));



  for (NodeIt k(_g); k != lemon::INVALID; ++k)
  {
    (*_components)[component_labels[k]].push_back(k);
  }

  for (int i=0; i< n_components; i++)
  {
    unsigned long size = (*_components)[i].size();
    for (int j=0; j< size; j++)
    {
      Node u1v1 = (*_components)[i][j];
      Node u1 = _gToProt1[u1v1];
      Node v1 = _gToProt2[u1v1];
      for (int k=j+1; k< size;k++)
      {
        Node u2v2 = (*_components)[i][k];

        if (u1v1 == u2v2)
          continue;

        Node u2 = _gToProt1[u2v2];
        Node v2 = _gToProt2[u2v2];


        if (u1 != u2 && v1 != v2)
        {
          bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
          bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

          if (!u1u2 && !v1v2)
          {
            _edgeType[_g.addEdge(u1v1, u2v2)] = PRODUCT_BLACK_EDGE;
            ++_numEdges;
          } else if (!u1u2 || !v1v2)
          {

//            _edgeType[_g.addEdge(u1v1, u2v2)] = PRODUCT_BLUE_EDGE;
//            ++_numEdges;
          }
        }
      }
    }
  }

  int component_sizes [n_components];
  memset(component_sizes, 0, sizeof component_sizes);


  int max_size = 0;
  for (NodeIt k(_g); k != lemon::INVALID; ++k)
  {
    component_sizes[component_labels[k]] = component_sizes[component_labels[k]] + 1;
    if (component_sizes[component_labels[k]] > max_size)
    {
      max_size = component_sizes[component_labels[k]];
    }
  }


  int size_frequencies [max_size+1];
  memset(size_frequencies, 0, sizeof size_frequencies);
  for (int m=0; m<n_components; m++)
  {
    ++size_frequencies[component_sizes[m]];
  }


  std::ofstream freq_file;
  freq_file.open("/Users/jelmer/Desktop/comp_freq.dot");

  for (int n=0; n<=max_size; n++)
  {
    if (size_frequencies[n] > 0)
      freq_file << n << ": " << size_frequencies[n] << std::endl;
  }
  
  freq_file.close();




}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::generate(const ProteinType &prot,
                                           const IntNodeMap& deg,
                                           IntSetNodeMap& degSet)
{
  const Graph& g = prot.getGraph();
  BoolNodeMap visited(g, false);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    lemon::mapFill(g, visited, false);
    IntSet& ds = degSet[v];
    visited[v] = true;
    ds.insert(deg[v]);
  }
}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::setComponentState(int index, bool state)
{
  NodeVector &component = (*_components)[index];
  unsigned long size = component.size();

  for (int i = 0; i < size; i++)
  {
    _nodeFilter[component[i]] = state;
  }
}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::printDOT(std::ostream& out) const
{
  // header
  out << "graph G {" << std::endl
  << "\toverlap=scale" << std::endl
  << "\tlayout=neato" << std::endl;

  // nodes
  for (NodeIt uv(_g); uv != lemon::INVALID; ++uv)
  {
    Node u = _gToProt1[uv];
    Node v = _gToProt2[uv];
    out << "\t" << _g.id(uv) << " [label=\"[" << _prot1.getLabel(u) << " -> "
    << _prot2.getLabel(v) << "]" << "\\n" << _g.id(uv) << "\"]" << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
    if (connectivityEdge(e) == PRODUCT_RED_EDGE)
    {
      out << " [color=red]";
    }
    else if (connectivityEdge(e) == PRODUCT_BLUE_EDGE)
    {
      out << " [color=blue]";
    }
    out << std::endl;
  }

  out << "}" << std::endl;
}

} // namespace cgp
} // namespace nina

#endif // PRODUCT_H
