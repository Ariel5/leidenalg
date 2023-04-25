#include "ModularityVertexPartition.h"

#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

ModularityVertexPartition::ModularityVertexPartition(Graph* graph,
      vector<size_t> const& membership) :
        MutableVertexPartition(graph,
        membership)
{ }

ModularityVertexPartition::ModularityVertexPartition(Graph* graph) :
        MutableVertexPartition(graph)
{ }

ModularityVertexPartition::~ModularityVertexPartition()
{ }

ModularityVertexPartition* ModularityVertexPartition::create(Graph* graph)
{
  return new ModularityVertexPartition(graph);
}

ModularityVertexPartition* ModularityVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
  return new ModularityVertexPartition(graph, membership);
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
  TODO Ariel: do we try this for all neighbors? Where is this called?
*****************************************************************************/
double ModularityVertexPartition::diff_move(size_t v, size_t new_comm)
{
  #ifdef DEBUG
    cerr << "double ModularityVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
  #endif
  size_t old_comm = this->_membership[v]; // Current community vertex v belongs to
  double diff = 0.0;

  // So if is_directed, multiply x1. If undirected, multiply all edges x2 - makes sense
  double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());
  if (total_weight == 0.0)
    return 0.0;
  if (new_comm != old_comm)
  {
    // How "heavy" the connection is to each of the 2 communities
    // TODO Ariel is this computed everytime? Hopefully it's stored - source stays same in loop
    double w_to_old = this->weight_to_comm(v, old_comm);
    double w_from_old = this->weight_from_comm(v, old_comm);
    double w_to_new = this->weight_to_comm(v, new_comm);
    double w_from_new = this->weight_from_comm(v, new_comm);

    // Sum of weights of outgoing edges TODO for this node or community?
    double k_out = this->graph->strength(v, IGRAPH_OUT);
    // -''- incoming edges for this node
    double k_in = this->graph->strength(v, IGRAPH_IN);

    // ? I think this is for after communities are aggregated
    // Shows the sum of weights of aggregated inner edges
    double self_weight = this->graph->node_self_weight(v);

    // What's the total edge weught going out of this comm.? Probably to tell how much an impact
    //  our current node has vs. total
    double K_out_old = this->total_weight_from_comm(old_comm);
    double K_in_old = this->total_weight_to_comm(old_comm);

    // Also how much will it be affected by us adding this node
    double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
    double K_in_new = this->total_weight_to_comm(new_comm) + k_in;

    // Modularity delta? TODO Ariel
    // Where is this formula?
    double diff_old = (w_to_old - k_out*K_in_old/total_weight) + \
               (w_from_old - k_in*K_out_old/total_weight);
    double diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) + \
               (w_from_new + self_weight - k_in*K_out_new/total_weight);

    diff = diff_new - diff_old;
  }
  #ifdef DEBUG
    cerr << "exit double ModularityVertexPartition::diff_move((" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
  #endif

  double m; // Is this Volume? Why not use total_weight that they defined above? Seems identical
  if (this->graph->is_directed())
    m = this->graph->total_weight();
  else
    m = 2*this->graph->total_weight();

  return diff/m;
}


/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double ModularityVertexPartition::quality()
{
  #ifdef DEBUG
    cerr << "double ModularityVertexPartition::quality()" << endl;
  #endif
  double mod = 0.0;

  double m;
  if (this->graph->is_directed())
    m = this->graph->total_weight();
  else
    m = 2*this->graph->total_weight(); // If Undirected, count each edge weight twice (both ways)

  if (m == 0)
    return 0.0;

  // Ariel - Modularity computed here for all communities
  for (size_t c = 0; c < this->n_communities(); c++)
  {
    // Ariel Seems like a push-pull + self-update
    double w = this->total_weight_in_comm(c); // TODO What do these 3 lines mean
    double w_out = this->total_weight_from_comm(c); // push-modularity?
    double w_in = this->total_weight_to_comm(c); // pull-modularity?
    #ifdef DEBUG
      double csize = this->csize(c);
      cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
    #endif

    // Ariel - This is the crucial line
    // TODO Is graph.total_weight global? Seems so. Can we make calc. local?
      // TODO Why divide InWeight by 4 when Undirected?
    mod += w - w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
  }
  double q = (2.0 - this->graph->is_directed())*mod; // <=> if directed: return mod ; else return 2*mod
  #ifdef DEBUG
    cerr << "exit double ModularityVertexPartition::quality()" << endl;
    cerr << "return " << q/m << endl << endl;
  #endif
  return q/m;
}
