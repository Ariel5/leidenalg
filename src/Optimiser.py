from . import _c_louvain
from .VertexPartition import LinearResolutionParameterVertexPartition
from collections import namedtuple
from math import log, sqrt
import sys

# Check if working with Python 3
PY3 = (sys.version > '3');

class Optimiser(object):
  """ Class for doing community detection using the Louvain algorithm.

  The optimiser class provides a number of different methods for optimising a
  given partition. The overall optimise procedure :func:`optimise_partition`
  calls either :func:`move_nodes` or :func:`merge_nodes` (which is controlled
  by :attr:`optimise_routine`) then aggregates the graph and repeats the same
  procedure. Possible, indicated by :attr:`refine_partition` the partition is
  refined before aggregating, meaning that subsets of communities are
  considered for moving around. Which routine is used for the refinement is
  indicated by :attr:`refine_routine`.  For calculating the actual improvement
  of moving a node (corresponding a subset of nodes in the aggregate graph),
  the code relies on :func:`~VertexPartition.MutableVertexPartition.diff_move`
  which provides different values for different methods (e.g. modularity or
  CPM). Finally, the Optimiser class provides a routine to construct a
  :func:`resolution_profile` on a resolution parameter.
  """
  def __init__(self):
    """ Create a new Optimiser object """
    self._optimiser = _c_louvain._new_Optimiser();

  #########################################################3
  # consider_comms
  @property
  def consider_comms(self):
    """ Determine how alternative communities are considered for moving
    a node for *optimising* a partition. 
    
    Nodes will only move to alternative communities that improve the given
    quality function.

    Notes
    -------
    This attribute should be set to one of the following values

    * :attr:`louvain.ALL_NEIGH_COMMS`
      Consider all neighbouring communities for moving.

    * :attr:`louvain.ALL_COMMS`
      Consider all communities for moving. This is especially useful in the
      case of negative links, in which case it may be better to move a node to
      a non-neighbouring community.

    * :attr:`louvain.RAND_NEIGH_COMM`
      Consider a random neighbour community for moving. The probability to
      choose a community is proportional to the number of neighbours a node has
      in that community.

    * :attr:`louvain.RAND_COMM`
      Consider a random community for moving. The probability to choose a
      community is proportional to the number of nodes in that community.
    """
    return _c_louvain._Optimiser_get_consider_comms(self._optimiser);

  @consider_comms.setter
  def consider_comms(self, value):
    _c_louvain._Optimiser_set_consider_comms(self._optimiser, value);

  #########################################################3
  # refine consider_comms
  @property
  def refine_consider_comms(self):
    """ Determine how alternative communities are considered for moving
    a node when *refining* a partition. 
    
    Nodes will only move to alternative communities that improve the given
    quality function.

    Notes
    -------
    This attribute should be set to one of the following values

    * :attr:`louvain.ALL_NEIGH_COMMS`
      Consider all neighbouring communities for moving.

    * :attr:`louvain.ALL_COMMS`
      Consider all communities for moving. This is especially useful in the
      case of negative links, in which case it may be better to move a node to
      a non-neighbouring community.

    * :attr:`louvain.RAND_NEIGH_COMM`
      Consider a random neighbour community for moving. The probability to
      choose a community is proportional to the number of neighbours a node has
      in that community.

    * :attr:`louvain.RAND_COMM`
      Consider a random community for moving. The probability to choose a
      community is proportional to the number of nodes in that community.
    """
    return _c_louvain._Optimiser_get_refine_consider_comms(self._optimiser);

  @refine_consider_comms.setter
  def refine_consider_comms(self, value):
    _c_louvain._Optimiser_set_refine_consider_comms(self._optimiser, value);

  #########################################################3
  # optimise routine
  @property
  def optimise_routine(self):
    """ Determine the routine to use for *optimising* a partition.

    Notes
    -------
    This attribute should be set to one of the following values

    * :attr:`louvain.MOVE_NODES`
      Use :func:`~Optimiser.move_nodes`.

    * :attr:`louvain.MERGE_NODES`
      Use :func:`~Optimiser.merge_nodes`.
    """
    return _c_louvain._Optimiser_get_optimise_routine(self._optimiser);

  @optimise_routine.setter
  def optimise_routine(self, value):
    _c_louvain._Optimiser_set_optimise_routine(self._optimiser, value);

  #########################################################3
  # optimise routine
  @property
  def refine_routine(self):
    """ Determine the routine to use for *refining* a partition.

    Notes
    -------
    This attribute should be set to one of the following values

    * :attr:`louvain.MOVE_NODES`
      Use :func:`~Optimiser.move_nodes`.

    * :attr:`louvain.MERGE_NODES`
      Use :func:`~Optimiser.merge_nodes`.
    """
    return _c_louvain._Optimiser_get_refine_routine(self._optimiser);

  @refine_routine.setter
  def refine_routine(self, value):
    _c_louvain._Optimiser_set_refine_routine(self._optimiser, value);

  #########################################################3
  # refine_partition
  @property
  def refine_partition(self):
    """ boolean: if ``True`` refine partition before aggregation. """
    return _c_louvain._Optimiser_get_refine_partition(self._optimiser);

  @refine_partition.setter
  def refine_partition(self, value):
    _c_louvain._Optimiser_set_refine_partition(self._optimiser, value);

  def optimise_partition(self, partition):
    """ Optimise the given partition.

    Parameters
    ----------
    partition
      The :class:`~louvain.VertexPartition.MutableVertexPartition` to optimise.

    Returns
    -------
    float
      Improvement in quality function.

    Examples
    --------

    >>> G = ig.Graph.Famous('Zachary');
    >>> optimiser = louvain.Optimiser();
    >>> partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.optimise_partition(partition);

    """
    # Perhaps we
    diff = _c_louvain._Optimiser_optimise_partition(self._optimiser, partition._partition);
    partition._update_internal_membership();
    return diff;

  def optimise_partition_multiplex(self, partitions, layer_weights=None):
    """ Optimise the given partitions simultaneously.

    Parameters
    ----------
    partitions
      List of :class:`~louvain.VertexPartition.MutableVertexPartition` layers to optimise.

    layer_weights
      List of weights of layers.

    Returns
    -------
    float
      Improvement in quality of combined partitions, see `Notes`_.

    Notes
    -----
    .. _Notes:

    This method assumes that the partitions are defined for graphs with the
    same vertices. The connections between the vertices may be different, but
    the vertices themselves should be identical. In other words, all vertices
    should have identical indices in all graphs (i.e. node `i` is assumed to be
    the same node in all graphs). The quality of the overall partition is
    simply the sum of the individual qualities for the various partitions,
    weighted by the layer_weight. If we denote by :math:`Q_k` the quality of
    layer :math:`k` and the weight by :math:`\\lambda_k`, the overall quality is then

    .. math:: Q = \sum_k \\lambda_k Q_k.

    This is particularly useful for graphs containing negative links. When
    separating the graph in two graphs, the one containing only the positive
    links, and the other only the negative link, by supplying a negative weight
    to the latter layer, we try to find relatively many positive links within a
    community and relatively many negative links between communities. Note that
    in this case it may be better to assign a node to a community to which it
    is not connected so that :attr:`consider_comms` may be better set to
    :attr:`louvain.ALL_COMMS`.

    Besides multiplex graphs where each node is assumed to have a single
    community, it is also useful in the case of for example multiple time
    slices, or in situations where nodes can have different communities in
    different slices. The package includes some special helper functions for
    using :func:`optimise_partition_multiplex` in such cases, where there is a
    conversion required from (time) slices to layers suitable for use in this
    function.

    See Also
    --------
    louvain.slices_to_layers : Convert slices to layers.
    louvain.time_slices_to_layers : Convert time slices to layers.
    louvain.find_partition_temporal : Detect partition for time slices.

    Examples
    --------
    >>> G_pos = ig.Graph.SBM(100, pref_matrix=[[0.5, 0.1], [0.1, 0.5]], block_sizes=[50, 50]);
    >>> G_neg = ig.Graph.SBM(100, pref_matrix=[[0.1, 0.5], [0.5, 0.1]], block_sizes=[50, 50]);
    >>> optimiser = louvain.Optimiser();
    >>> partition_pos = louvain.ModularityVertexPartition(G_pos);
    >>> partition_neg = louvain.ModularityVertexPartition(G_neg);
    >>> optimiser.optimise_partition_multiplex(partitions=[partition_pos, partition_neg],
    ...                                        layer_weights=[1,-1]);

    """
    if not layer_weights:
      layer_weights = [1]*len(partitions);
    diff = _c_louvain._Optimiser_optimise_partition_multiplex(
      self._optimiser,
      [partition._partition for partition in partitions],
      layer_weights);
    for partition in partitions:
      partition._update_internal_membership();
    return diff;

  def move_nodes(self, partition, consider_comms=None):
    """ Move nodes to alternative communities for *optimising* the partition.

    Parameters
    ----------
    partition
      The partition for which to move nodes.

    consider_comms
      If ``None`` uses :attr:`~Optimiser.consider_comms`, but can be set to
      something else.

    Returns
    -------
    float
      Improvement in quality function.

    Notes
    -----
    When moving nodes, the function loops over nodes and considers moving the
    node to an alternative community. Which community depends on
    ``consider_comms``. The function terminates when no more nodes can be moved
    to an alternative community.

    See Also
    --------
    move_nodes_constrained : Move nodes for *refining* a partition.
    merge_nodes : Merge nodes rather than moving nodes.

    Examples
    --------
    >>> G = ig.Graph.Famous('Zachary');
    >>> optimiser = louvain.Optimiser();
    >>> partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.move_nodes(partition);

    """
    if (consider_comms is None):
      consider_comms = self.consider_comms;
    diff =  _c_louvain._Optimiser_move_nodes(self._optimiser, partition._partition, consider_comms);
    partition._update_internal_membership();
    return diff;

  def move_nodes_constrained(self, partition, constrained_partition, consider_comms=None):
    """ Move nodes to alternative communities for *refining* the partition.

    Parameters
    ----------
    partition
      The partition for which to move nodes.

    constrained_partition
      The partition within which we may move nodes.

    consider_comms
      If ``None`` uses :attr:`~Optimiser.refine_consider_comms`, but can be set
      to something else.

    Returns
    -------
    float
      Improvement in quality function.

    Notes
    -----
    The idea is constrain the movement of nodes to alternative communities to
    another partition. In other words, if there is a partition ``P`` which we
    want to refine, we can then initialize a new singleton partition, and move
    nodes in that partition constrained to ``P``.

    See Also
    --------
    move_nodes : Move nodes for *optimising* a partition.
    merge_nodes_constrained : Merge nodes rather than moving nodes.

    Examples
    --------
    >>> G = ig.Graph.Famous('Zachary');
    >>> optimiser = louvain.Optimiser();
    >>> partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.optimise_partition(partition);
    >>> refine_partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.move_nodes_constrained(refine_partition, partition);

    """
    if (consider_comms is None):
      consider_comms = self.refine_consider_comms;
    diff =  _c_louvain._Optimiser_move_nodes_constrained(self._optimiser, partition._partition, constrained_partition._partition, consider_comms);
    partition._update_internal_membership();
    return diff;

  def merge_nodes(self, partition, consider_comms=None):
    """ Merge nodes for *optimising* the partition.

    Parameters
    ----------
    partition
      The partition for which to merge nodes.

    consider_comms
      If ``None`` uses :attr:`~Optimiser.consider_comms`, but can be set to
      something else.

    Returns
    -------
    float
      Improvement in quality function.

    Notes
    -----
    This function loop over all nodes once and tries to merge them with another
    community.  Merging in this case implies that a node will never be removed
    from a community, only merged with other communities.

    See Also
    --------
    move_nodes_constrained : Move nodes for *refining* a partition.
    merge_nodes : Merge nodes rather than moving nodes.

    Examples
    --------
    >>> G = ig.Graph.Famous('Zachary');
    >>> optimiser = louvain.Optimiser();
    >>> partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.merge_nodes(partition);

    """
    if (consider_comms is None):
      consider_comms = self.consider_comms;
    diff =  _c_louvain._Optimiser_merge_nodes(self._optimiser, partition._partition, consider_comms);
    partition._update_internal_membership();
    return diff;

  def merge_nodes_constrained(self, partition, constrained_partition, consider_comms=None):
    """ Merge nodes for *refining* the partition.

    Parameters
    ----------
    partition
      The partition for which to merge nodes.

    constrained_partition
      The partition within which we may merge nodes.

    consider_comms
      If ``None`` uses :attr:`~Optimiser.refine_consider_comms`, but can be set
      to something else.

    Returns
    -------
    float
      Improvement in quality function.

    Notes
    -----
    The idea is constrain the merging of nodes to another partition. In other
    words, if there is a partition ``P`` which we want to refine, we can then
    initialize a new singleton partition, and move nodes in that partition
    constrained to ``P``.

    See Also
    --------
    move_nodes : Move nodes for *optimising* a partition.
    merge_nodes_constrained : Merge nodes rather than moving nodes.

    Examples
    --------
    >>> G = ig.Graph.Famous('Zachary');
    >>> optimiser = louvain.Optimiser();
    >>> partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.optimise_partition(partition);
    >>> refine_partition = louvain.ModularityVertexPartition(G);
    >>> optimiser.move_nodes_constrained(refine_partition, partition);

    """
    if (consider_comms is None):
      consider_comms = self.refine_consider_comms;
    diff =  _c_louvain._Optimiser_merge_nodes_constrained(self._optimiser, partition._partition, constrained_partition._partition, consider_comms);
    partition._update_internal_membership();
    return diff;

  def resolution_profile(self,
        graph,
        partition_type,
        resolution_range,
        weights=None,
        bisect_func=lambda p: p.bisect_value(),
        min_diff_bisect_value=1,
        min_diff_resolution=1e-3,
        linear_bisection=False,
        until_stable=False
        ):
    """ Use bisectioning on the resolution parameter in order to construct a
    resolution profile.

    Parameters
    ----------
    graph
      The graph for which to construct a resolution profile.

    partition_type
      The type of :class:`~louvain.VertexPartition.MutableVertexPartition` used
      to find a partition (must support resolution parameters obviously).

    resolution_range
      The range of resolution values that we would like to scan.

    weights
      If provided, indicates the edge attribute to use as a weight.

    Returns
    -------
    list of :class:`~louvain.VertexPartition.MutableVertexPartition` 
      A list of partitions for different resolutions.

    Other Parameters
    ----------------
    bisect_func
      The function used for bisectioning. For the methods currently
      implemented, this should usually not be altered.

    min_diff_bisect_value
      The difference in the value returned by the bisect_func below which the
      bisectioning stops (i.e. by default, a difference of a single edge does
      not trigger further bisectioning).

    min_diff_resolution
      The difference in resolution below which the bisectioning stops. For
      positive differences, the logarithmic difference is used by default, i.e.
      ``diff = log(res_1) - log(res_2) = log(res_1/res_2)``, for which ``diff >
      min_diff_resolution`` to continue bisectioning. Set the linear_bisection
      to true in order to use only linear bisectioning (in the case of negative
      resolution parameters for example, which can happen with negative
      weights).

    linear_bisection
      Whether the bisectioning will be done on a linear or on a logarithmic
      basis (if possible).

    until_stable
      If ``True`` iterate ``optimise_partition`` until no more improvement can
      be found.
    """

    # Helper function for cleaning values to be a stepwise function
    def clean_stepwise(bisect_values):
      # Check best partition for each resolution parameter
      for res, bisect in bisect_values.iteritems():
        best_bisect = bisect;
        best_quality = bisect.partition.quality(res);
        for res2, bisect2 in bisect_values.iteritems():
          if bisect2.partition.quality(res) > best_quality:
            best_bisect = bisect2;
            best_quality = bisect2.partition.quality(res);
        if best_bisect != bisect:
          bisect_values[res] = best_bisect;

      # We only need to keep the changes in the bisection values
      bisect_list = sorted([(res, part.bisect_value) for res, part in
        bisect_values.iteritems()], key=lambda x: x[0]);
      for (res1, v1), (res2, v2) \
          in zip(bisect_list,
                 bisect_list[1:]):
        # If two consecutive bisection values are the same, remove the second
        # resolution parameter
        if v1 == v2:
          del bisect_values[res2];

      for res, bisect in bisect_values.iteritems():
        bisect.partition.resolution_parameter = res;

    # We assume here that the bisection values are
    # monotonically decreasing with increasing resolution
    # parameter values.
    def ensure_monotonicity(bisect_values, new_res):
      # First check if this partition improves on any other partition
      for res, bisect_part in bisect_values.iteritems():
        if bisect_values[new_res].partition.quality(res) > bisect_part.partition.quality(res):
          bisect_values[res] = bisect_values[new_res];
      # Then check what is best partition for the new_res
      current_quality = bisect_values[new_res].partition.quality(new_res)
      best_res = new_res;
      for res, bisect_part in bisect_values.iteritems():
        if bisect_part.partition.quality(new_res) > current_quality:
          best_res = new_res;
      bisect_values[new_res] = bisect_values[best_res];

    def find_partition(self, graph, partition_type, weights=None, until_stable=False,**kwargs):
      partition = partition_type(graph,
                             weights=weights,
                             **kwargs);
      if until_stable:
        while self.optimise_partition(partition) > 0:
          pass;
      else:
        self.optimise_partition(partition)
      return partition;
    assert(issubclass(partition_type, LinearResolutionParameterVertexPartition),
        "Bisectioning only works on partitions with a linear resolution parameter.");
    # Start actual bisectioning
    bisect_values = {};
    stack_res_range = [];
    # Push first range onto the stack
    stack_res_range.append(resolution_range);
    # Make sure the bisection values are calculated
    # The namedtuple we will use in the bisection function
    BisectPartition = namedtuple('BisectPartition',
        ['partition', 'bisect_value']);
    partition = find_partition(self, graph=graph, partition_type=partition_type,
        weights=weights,until_stable=until_stable,resolution_parameter=resolution_range[0]);
    bisect_values[resolution_range[0]] = BisectPartition(partition=partition,
                                bisect_value=bisect_func(partition));
    partition = find_partition(self, graph=graph, partition_type=partition_type,
        weights=weights,until_stable=until_stable, resolution_parameter=resolution_range[1]);
    bisect_values[resolution_range[1]] = BisectPartition(partition=partition,
                                bisect_value=bisect_func(partition));
    # While stack of ranges not yet empty
    while stack_res_range:
      # Get the current range from the stack
      current_range = stack_res_range.pop();
      # Get the difference in bisection values
      diff_bisect_value = abs(bisect_values[current_range[0]].bisect_value -
                              bisect_values[current_range[1]].bisect_value);
      # Get the difference in resolution parameter (in log space if 0 is not in
      # the interval (assuming only non-negative resolution parameters).
      if current_range[0] > 0 and current_range[1] > 0 and not linear_bisection:
        diff_resolution = log(current_range[1]/current_range[0]);
      else:
        diff_resolution = abs(current_range[1] - current_range[0]);
      # Check if we still want to scan a smaller interval
      # If we would like to bisect this interval
      if diff_bisect_value > min_diff_bisect_value and \
         diff_resolution > min_diff_resolution:
        # Determine new resolution value
        if current_range[0] > 0 and current_range[1] > 0 and not linear_bisection:
          new_res = sqrt(current_range[1]*current_range[0]);
        else:
          new_res = sum(current_range)/2.0;
        # Bisect left (push on stack)
        stack_res_range.append((current_range[0], new_res));
        # Bisect right (push on stack)
        stack_res_range.append((new_res, current_range[1]));
        # If we haven't scanned this resolution value yet,
        # do so now
        if not bisect_values.has_key(new_res):
          partition = find_partition(self, graph, partition_type=partition_type,
              weights=weights, until_stable=until_stable,resolution_parameter=new_res);
          bisect_values[new_res] = BisectPartition(partition=partition,
                                      bisect_value=bisect_func(partition));
          # Because of stochastic differences in different runs, the monotonicity
          # of the bisection values might be violated, so check for any
          # inconsistencies
          ensure_monotonicity(bisect_values, new_res);

    # Ensure we only keep those resolution values for which
    # the bisection values actually changed, instead of all of them
    clean_stepwise(bisect_values);
    # Use an ordered dict so that when iterating over it, the results appear in
    # increasing order based on the resolution value.
    return sorted((bisect.partition for res, bisect in
      bisect_values.iteritems()), key=lambda x: x.resolution_parameter);
