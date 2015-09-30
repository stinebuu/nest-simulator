/*
 *  node_manager.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef NODE_MANAGER_H
#define NODE_MANAGER_H

#include <omp.h>

#include "nest.h"
#include "manager_interface.h"
#include "dict.h"
#include "communicator.h"
#include "arraydatum.h"
#include "sparse_node_array.h"

namespace nest
{

class SiblingContainer;
class Node;

class NodeManager : ManagerInterface
{
public:
  NodeManager();
  ~NodeManager();

  virtual void init();
  virtual void reset();
  virtual void reinit_nodes();

  virtual void set_status( const Dictionary& );
  virtual void get_status( Dictionary& );

  /**
   * Get properties of a node. The specified node must exist.
   * @throws nest::UnknownNode       Target does not exist in the network.
   */
  DictionaryDatum get_status( index );

  /**
   * Add a number of nodes to the network.
   * This function creates n Node objects of Model m and adds them
   * to the Network at the current position.
   * @param m valid Model ID.
   * @param n Number of Nodes to be created. Defaults to 1 if not
   * specified.
   * @throws nest::UnknownModelID
   */
  index add_node( index m, long_t n = 1 );


  /**
   * Restore nodes from an array of status dictionaries.
   * The following entries must be present in each dictionary:
   * /model - with the name or index of a neuron mode.
   *
   * The following entries are optional:
   * /parent - the node is created in the parent subnet
   *
   * Restore nodes uses the current working node as root. Thus, all
   * GIDs in the status dictionaties are offset by the GID of the current
   * working node. This allows entire subnetworks to be copied.
   */
  void restore_nodes( ArrayDatum& );


  /**
   * Set the state (observable dynamic variables) of a node to model defaults.
   * @see Node::init_state()
   */
  void init_state( index );

  /**
   * Return total number of network nodes.
   * The size also includes all Subnet objects.
   */
  index size() const;

  Subnet* get_root() const; ///< return root subnet.
  Subnet* get_cwn() const;  ///< current working node.

  /**
   * Change current working node. The specified node must
   * exist and be a subnet.
   * @throws nest::IllegalOperation Target is no subnet.
   */
  void go_to( index );

  void print( index, int );

  /**
   * Return true, if the given Node is on the local machine
   */
  bool is_local_node( Node* ) const;

  /**
   * Return true, if the given gid is on the local machine
   */
  bool is_local_gid( index gid ) const;

  /**
   * Return pointer of the specified Node.
   * @param i Index of the specified Node.
   * @param thr global thread index of the Node.
   *
   * @throws nest::UnknownNode       Target does not exist in the network.
   *
   * @ingroup net_access
   */
  Node* get_node( index, thread thr = 0 );

/**
   * Set properties of a Node. The specified node must exist.
   * @throws nest::UnknownNode       Target does not exist in the network.
   * @throws nest::UnaccessedDictionaryEntry  Non-proxy target did not read dict entry.
   * @throws TypeMismatch            Array is not a flat & homogeneous array of integers.
   */
  void set_status( index, const DictionaryDatum& );

  /**
   * Return the Subnet that contains the thread siblings.
   * @param i Index of the specified Node.
   *
   * @throws nest::NoThreadSiblingsAvailable     Node does not have thread siblings.
   *
   * @ingroup net_access
   */
  const SiblingContainer* get_thread_siblings( index n ) const;

  /**
   * Ensure that all nodes in the network have valid thread-local IDs.
   * Create up-to-date vector of local nodes, nodes_vec_.
   * This method also sets the thread-local ID on all local nodes.
   */
  void ensure_valid_thread_local_ids();

  Node* thread_lid_to_node( thread t, targetindex thread_local_id ) const;

  /**
   * Increment total number of global spike detectors by 1
   */
  void increment_n_gsd();

  /**
   * Get total number of global spike detectors
   */
  index get_n_gsd();

  /**
   * Get node lists for unfrozen nodes
   */
  vector< vector< Node* > > get_nodes_vec();


  /**
   * Prepare nodes for simulation and register nodes in node_list.
   * Calls prepare_node_() for each pertaining Node.
   * @see prepare_node_()
   */
  void prepare_nodes();

  /**
   * Invoke finalize() on nodes registered for finalization.
   */
  void finalize_nodes();

private:
  /**
   * Initialize the network data structures.
   * init_() is used by the constructor and by reset().
   * @see reset()
   */
  void init_();
  void destruct_nodes_();

  /**
   * Helper function to set properties on single node.
   * @param node to set properties for
   * @param dictionary containing properties
   * @param if true (default), access flags are called before
   *        each call so Node::set_status_()
   * @throws UnaccessedDictionaryEntry
   */
  void set_status_single_node_( Node&, const DictionaryDatum&, bool clear_flags = true );

  /**
   * Initialized buffers, register in list of nodes to update/finalize.
   * @see prepare_nodes_()
   */
  void prepare_node_( Node* );

private:
  SparseNodeArray local_nodes_; //!< The network as sparse array of local nodes
  Subnet* root_;    //!< Root node.
  Subnet* current_; //!< Current working node (for insertion).

  Model* siblingcontainer_model; //!< The model for the SiblingContainer class

  index n_gsd_; //!< Total number of global spike detectors, used for distributing them over
                //!< recording processes

  vector< vector< Node* > > nodes_vec_; //!< Nodelists for unfrozen nodes
  index nodes_vec_network_size_;        //!< Network size when nodes_vec_ was last updated




};

inline index
NodeManager::size() const
{
  return local_nodes_.get_max_gid() + 1;
}

inline Subnet*
NodeManager::get_root() const
{
  return root_;
}

inline Subnet*
NodeManager::get_cwn( void ) const
{
  return current_;
}

inline bool
NodeManager::is_local_gid( index gid ) const
{
  return local_nodes_.get_node_by_gid( gid ) != 0;
}

inline Node*
NodeManager::thread_lid_to_node( thread t, targetindex thread_local_id ) const
{
  return nodes_vec_[ t ][ thread_local_id ];
}

inline void
NodeManager::prepare_node_( Node* n )
{
  // Frozen nodes are initialized and calibrated, so that they
  // have ring buffers and can accept incoming spikes.
  n->init_buffers();
  n->calibrate();
}

inline void
NodeManager::increment_n_gsd()
{
  ++n_gsd_;
}

inline index
NodeManager::get_n_gsd()
{
  return n_gsd_;
}

inline vector< vector< Node* > >
NodeManager::get_nodes_vec()
{
  return nodes_vec_;
}

} // namespace 

#endif /* NODE_MANAGER_H */
