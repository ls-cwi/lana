/*
 * bronkerboschconnected.h
 *
 *  Created on: 23-jan-2014
 *      Author: M. El-Kebir
 */

#ifndef BRONKERBOSCHCONNECTED_H
#define BRONKERBOSCHCONNECTED_H

#include "bronkerbosch.h"
#include "product.h"
#include "options.h"
#include <boost/dynamic_bitset.hpp>
#include <lemon/core.h>

namespace nina {
namespace gna {

template<typename GR, typename BGR, typename PGR>
class BronKerboschConnected : public nina::BronKerbosch<GR,BGR,PGR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// The graph type of the bipartite matching graph
  typedef BGR BpGraph;
  /// The graph type of the first parameter type of the ProductGraph.
  typedef PGR ProductGraphType;
  /// Base class type
  typedef nina::BronKerbosch<GR,BGR,PGR> Parent;
  /// Product graph type
  typedef Product<ProductGraphType , BpGraph> ProductType;
  /// Solver type
  typedef typename Parent::SolverType SolverType;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Parent::BitSetNodeMap BitSetNodeMap;
  typedef typename Parent::BitSet BitSet;
  typedef typename Parent::NodeList NodeList;
  using Parent::_g;
  using Parent::_n;
  using Parent::_bitToNode;
  using Parent::_nodeToBit;
  using Parent::_bitNeighborhood;
  using Parent::computeDegeneracy;
  using Parent::report;
  using Parent::printBitSet;
  using Parent::_options;

public:
  BronKerboschConnected(const ProductType& product, const Options& options)
    : Parent(product, options)
    , _restrictedBitNeighborhood(_g, BitSet(_n))
  {
    // initialize restricted neighborhood mapping
    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      if (product.connectivityEdge(e) == ProductType::PRODUCT_RED_EDGE)
      {
        Node u = _g.u(e);
        Node v = _g.v(e);
        _restrictedBitNeighborhood[u][_nodeToBit[v]] = 1;
        _restrictedBitNeighborhood[v][_nodeToBit[u]] = 1;
      }
    }
  }

  virtual void run(SolverType type);

protected:
  BitSetNodeMap _restrictedBitNeighborhood;

private:
  void bkPivot(BitSet P, BitSet D, BitSet R, BitSet X, BitSet S);
};

template<typename GR, typename BGR, typename PGR>
void BronKerboschConnected<GR,BGR,PGR>::run(SolverType)
{
  NodeList order;
  std::cout << "Starting generation of Degeneracy" << std::endl;
  computeDegeneracy(order);
  std::cout << "Done with generation of Degeneracy" << std::endl;

  BitSet mask(_n);

  int i = 0;
  unsigned long size = order.size();
  for (typename NodeList::const_iterator it = order.begin(); it != order.end(); ++it)
  {
    std::cout << "Checked start node " << i << "/" << size << "(" << 100.*i/size <<  "%)" << std::endl;
    i++;
    Node v = *it;
    const BitSet& N_v = _bitNeighborhood[v];
    const BitSet& Nc_v = _restrictedBitNeighborhood[v];
    const BitSet Nd_v = N_v - Nc_v;

    // ~mask includes v but we're fine as N_v, Nc_v and Nd_v exclude v
    BitSet P = Nc_v & ~mask;
    BitSet D = Nd_v & ~mask;
    BitSet X = Nc_v & mask;
    BitSet S = Nd_v & mask;

    BitSet R(_n);
    R.set(_nodeToBit[v]);

    if (g_verbosity >= VERBOSE_DEBUG)
    {
      std::cerr << "P = ";
      printBitSet(P, std::cerr);
      std::cerr << ", D = ";
      printBitSet(D, std::cerr);
      std::cerr << ", R = ";
      printBitSet(R, std::cerr);
      std::cerr << ", X = ";
      printBitSet(X, std::cerr);
      std::cerr << ", S = ";
      printBitSet(S, std::cerr);
      std::cerr << std::endl;
    }

    bkPivot(P, D, R, X, S);
    mask.set(_nodeToBit[v]);
  }
}

template<typename GR, typename BGR, typename PGR>
void BronKerboschConnected<GR,BGR,PGR>::bkPivot(BitSet P, BitSet D,
                                        BitSet R,
                                        BitSet X, BitSet S)
{
  // all sets are pairwise disjoint
//  assert((P & D).none());
//  assert((P & R).none());
//  assert((P & X).none());
//  assert((P & S).none());
//  assert((D & R).none());
//  assert((D & X).none());
//  assert((D & S).none());
//  assert((R & X).none());
//  assert((R & S).none());
//  assert((X & S).none());

  if (g_verbosity >= VERBOSE_DEBUG)
  {
    std::cerr << "P = ";
    printBitSet(P, std::cerr);
    std::cerr << ", D = ";
    printBitSet(D, std::cerr);
    std::cerr << ", R = ";
    printBitSet(R, std::cerr);
    std::cerr << ", X = ";
    printBitSet(X, std::cerr);
    std::cerr << ", S = ";
    printBitSet(S, std::cerr);
    std::cerr << std::endl;
  }

  // Reports maximal c-cliques in P \cup R (but not in X and S)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none())
  {
    report(R);
  }
  //else if (P.none())
  //{
  //  report(R);
  //}
  else
  {
    // choose a pivot u from (P | X) s.t |P & N(u)| is maximum, Tomita et al. (2006)
    size_t maxBitCount = 0;
    Node max_u = lemon::INVALID;
    for (size_t i = 0; i < P.size(); ++i)
    {
      if (P_cup_X[i])
      {
        Node u = _bitToNode[i];
        BitSet P_cap_Nu = P & _bitNeighborhood[u];
        size_t s = P_cap_Nu.count();
        if (s >= maxBitCount)
        {
          max_u = u;
          maxBitCount = s;
        }
      }
    }

    assert(max_u != lemon::INVALID);

    BitSet P_diff_Nu = P - _bitNeighborhood[max_u];
//    std::cout << "Pivot: " << _g.id(max_u) << std::endl;
//    std::cout << "P diff NU: ";
//    printBitSet(P_diff_Nu, std::cout);
//    std::cout << std::endl;

    BitSet b(_n);

    for (size_t i = 0; i < P.size(); ++i)
    {
      if (P_diff_Nu[i])
      {
        b.set(i);
      }
      else if (P[i] && (D & (_bitNeighborhood[_bitToNode[i]] - _bitNeighborhood[max_u])).count() > 0)
      {
        b.set(i);
      }
    }


    for (size_t i = 0; i < P.size(); ++i)
    {
      if (b[i])
      {
        Node v = _bitToNode[i];
        const BitSet& N_v = _bitNeighborhood[v];
        const BitSet& Nc_v = _restrictedBitNeighborhood[v];

        BitSet P_ = P | ( D & Nc_v );
        BitSet D_ = D - Nc_v;

        BitSet R_ = R;
        R_[_nodeToBit[v]] = 1;

        BitSet X_ = X | ( S & Nc_v );
        BitSet S_ = S - Nc_v;

        // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
        bkPivot(P_ & N_v,
                D_ & N_v,
                R_,
                X_ & N_v,
                S_ & N_v);
        P[i] = 0;
        X[i] = 1;
      }
    }
  }

  /* Invariants:
   * - R is a c-clique
   * - Each node v in P is adjacent to all nodes in R
   *   and c-adjacent to a node in R
   *   (nodes in P are used to extend R, upon usage it's moved to X)
   * - Each node v in D is adjacent to all nodes in R
   *   but not c-adjacent to any node in R
   * - Each node v in X is adjacent to all nodes in R
   *   (maximal c-cliques containing R \cup v have already been reported)
   */
}


} // namespace cgp
} // namespace nina

#endif // BRONKERBOSCHRELAXED_H
