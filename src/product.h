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
#include "molecule.h"

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
  typedef Molecule<Graph> MoleculeType;
  typedef std::vector<typename Graph::Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  typedef std::set<typename Graph::Node> NodeSet;
  typedef std::pair<NodeSet, NodeSet> NodeSetPair;
  typedef std::set<NodeSetPair> NodeSetPairSet;

private:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Graph::template NodeMap<Node> NodeNodeMap;
  typedef std::multiset<int> IntSet;
  typedef typename Graph::template NodeMap<IntSet> IntSetNodeMap;
  typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;

public:
  Product(const MoleculeType& mol1, const MoleculeType& mol2, const MatchingGraphType& matchingGraph, int shell)
    : _mol1(mol1)
    , _mol2(mol2)
    , _matchingGraph(matchingGraph)
    , _shell(shell)
    , _g()
    , _mol1ToG(mol1.getGraph(), lemon::INVALID)
    , _mol2ToG(mol2.getGraph(), lemon::INVALID)
    , _gToMol1(_g)
    , _gToMol2(_g)
    , _connectivityEdge(_g)
    , _g1ToDeg1Neighbors(mol1.getGraph())
    , _g2ToDeg1Neighbors(mol2.getGraph())
    , _numNodes(0)
    , _numEdges(0)
  {
    generate();
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
        Node u = _gToMol1[uv];
        Node v = _gToMol2[uv];

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

  const Graph & getGraph() const { return _g; }

  bool connectivityEdge(Edge e) const { return _connectivityEdge[e]; }

private:
  const MoleculeType& _mol1;
  const MoleculeType& _mol2;
  const MatchingGraphType& _matchingGraph;
  const int _shell;
  Graph _g;
  NodeNodeMap _mol1ToG;
  NodeNodeMap _mol2ToG;
  NodeNodeMap _gToMol1;
  NodeNodeMap _gToMol2;
  BoolEdgeMap _connectivityEdge;
  NodeVectorMap _g1ToDeg1Neighbors;
  NodeVectorMap _g2ToDeg1Neighbors;

  int _numNodes;
  int _numEdges;

  void generate();

  void generate(const MoleculeType& mol,
                const IntNodeMap& deg,
                IntSetNodeMap& intSet,
                IntSetNodeMap& degSet);

  void generateDeg1NeighborSet(const Graph& g,
                               const IntNodeMap& deg,
                               NodeVectorMap& deg1NeighborMap);

  void dfs(const IntNodeMap& deg,
           const Node v, const int depth,
           const MoleculeType& mol,
           BoolNodeMap& visited,
           IntSet& s, IntSet& ds);

  void determineDegrees(const Graph& g, IntNodeMap& deg);

public:

  Node getNodeG1(Node uv) {
    return _gToMol1[uv];
  }

  Node getNodeG2(Node uv) {
    return _gToMol2[uv];
  }

  void printDOT(std::ostream& out) const;

  void printProductNodeJSON(Node uv, std::ostream& out) const
  {
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];

    out << "            {" << std::endl
    << "              \"id1\": " << _mol1.getLabel(u) << "," << std::endl
    << "              \"id2\": " << _mol2.getLabel(v) << "," << std::endl
    << "            }";

    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

    assert(uNeighbors.size() == vNeighbors.size());
    for (size_t i = 0; i < uNeighbors.size(); ++i)
    {
      out << "," << std::endl;
      out << "            {" << std::endl
      << "              \"id1\": " << _mol1.getLabel(uNeighbors[i]) << "," << std::endl
      << "              \"id2\": " << _mol2.getLabel(vNeighbors[i]) << "," << std::endl
      << "            }";
    }
  }

  void printProductNode(Node uv, std::ostream& out) const
  {
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];

    out << "[" << _g.id(uv) << ": " << _mol1.getLabel2(u) << " -> "
    << _mol2.getLabel2(v) << "]";

    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

    assert(uNeighbors.size() == vNeighbors.size());
    for (size_t i = 0; i < uNeighbors.size(); ++i)
    {
      out << ", ";
      out << "[" << _mol1.getLabel2(uNeighbors[i])
      << " (" << _mol1.getLabel(uNeighbors[i]) << ") , "
      << _mol2.getLabel2(vNeighbors[i])
      << " (" << _mol2.getLabel(vNeighbors[i]) << ")]";
    }
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
inline void Product<GR,BGR>::dfs(const IntNodeMap& deg,
                             const Node v, const int depth,
                             const MoleculeType& mol,
                             BoolNodeMap& visited,
                             IntSet& s, IntSet& ds)
{
  const Graph& g = mol.getGraph();
  visited[v] = true;
  s.insert(mol.getAtomType(v));
  ds.insert(deg[v]);

  if (depth < _shell)
  {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      Node w = g.oppositeNode(v, e);
      if (!visited[w])
      {
        dfs(deg, w, depth + 1, mol, visited, s, ds);
      }
    }
  }
}

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

  const Graph& g1 = _mol1.getGraph();
  const Graph& g2 = _mol2.getGraph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);
  lemon::ArcLookUp<BpGraph> arcLookUpGm(_matchingGraph.getGm());

  IntSetNodeMap set1(g1);
  IntSetNodeMap degSet1(g1);
  IntSetNodeMap set2(g2);
  IntSetNodeMap degSet2(g2);

  IntNodeMap deg1(g1, 0);
  determineDegrees(g1, deg1);
  IntNodeMap deg2(g2, 0);
  determineDegrees(g2, deg2);


  // determine degrees
  generate(_mol1, deg1, set1, degSet1);
  generate(_mol2, deg2, set2, degSet2);


  // generate nodes
  int64_t node_created = 0;
  for (NodeIt u(g1); u != lemon::INVALID; ++u)
  {
    for (NodeIt v(g2); v != lemon::INVALID; ++v)
    {

      BpNode bp_u = _matchingGraph.mapG1ToGm(u);
      BpNode bp_v = _matchingGraph.mapG2ToGm(v);
      BpEdge uv = arcLookUpGm(bp_u, bp_v);
      if ( set1[u] == set2[v] && uv != lemon::INVALID )//   && degSet1[u] == degSet2[v])
      {
        node_created++;
//        assert(deg1[u] == deg2[v]);
        // don't add product nodes for a pair of deg-1 nodes unless _shell == 0
        if (_shell == 0 || deg1[u] > 1)
        {
          Node uv = _g.addNode();
          _mol1ToG[u] = uv;
          _mol2ToG[v] = uv;
          _gToMol1[uv] = u;
          _gToMol2[uv] = v;
          ++_numNodes;
        }
      }
    }
  }

  // generate deg1 sets
  generateDeg1NeighborSet(g1, deg1, _g1ToDeg1Neighbors);
  generateDeg1NeighborSet(g2, deg2, _g2ToDeg1Neighbors);

  // generate c-edges ('red')
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1)
  {
    Node u1 = _gToMol1[u1v1];
    Node v1 = _gToMol2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2)
    {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _gToMol1[u2v2];
      Node v2 = _gToMol2[u2v2];

      assert(_mol1.getAtomType(u1) == _mol2.getAtomType(v1));

      if (u1 != u2 && v1 != v2)
      {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 && v1v2)
        {
          _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
          ++_numEdges;
        }
      }
    }
  }


  IntNodeMap components(_g);
  int n_components = lemon::connectedComponents(_g, components);
  std::vector<std::vector<Node> > component_lists(n_components, std::vector<Node>(0));



  for (NodeIt k(_g); k != lemon::INVALID; ++k)
  {
    component_lists[components[k]].push_back(k);
  }

  for (int i=0; i< n_components; i++)
  {
    unsigned long size = component_lists[i].size();
    // TODO: REMOVE THIS.
    if (size > 500) {continue;}
    for (int j=0; j< size; j++)
    {
      Node u1v1 = component_lists[i][j];
      Node u1 = _gToMol1[u1v1];
      Node v1 = _gToMol2[u1v1];
      for (int k=j+1; k< size;k++)
      {
        Node u2v2 = component_lists[i][k];

        if (u1v1 == u2v2)
          continue;

        Node u2 = _gToMol1[u2v2];
        Node v2 = _gToMol2[u2v2];

        assert(_mol1.getAtomType(u1) == _mol2.getAtomType(v1));

        if (u1 != u2 && v1 != v2)
        {
          bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
          bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

          if (!u1u2 && !v1v2)
          {
            _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
            ++_numEdges;
          }
        }
      }
    }
  }
  

}

template<typename GR, typename BGR>
inline void Product<GR,BGR>::generate(const MoleculeType& mol,
                                           const IntNodeMap& deg,
                                           IntSetNodeMap& intSet,
                                           IntSetNodeMap& degSet)
{
  const Graph& g = mol.getGraph();
  BoolNodeMap visited(g, false);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    lemon::mapFill(g, visited, false);
    IntSet& s = intSet[v];
    IntSet& ds = degSet[v];
    dfs(deg, v, 0, mol, visited, s, ds);
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
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];
    out << "\t" << _g.id(uv) << " [label=\"[" << _mol1.getLabel2(u) << " -> "
    << _mol2.getLabel2(v) << "]" << "\\n" << _g.id(uv) << "\"]" << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
    if (connectivityEdge(e))
    {
      out << " [color=red]";
    }
    out << std::endl;
  }

  out << "}" << std::endl;
}

} // namespace cgp
} // namespace nina

#endif // PRODUCT_H
