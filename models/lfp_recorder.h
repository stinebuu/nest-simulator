/*
 *  lfp_recorder.h
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

#ifndef LFP_RECORDER_H
#define LFP_RECORDER_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


/*
 * TODO: Documentation
 *
 */

namespace nest
{
/**
 * TODO: Documentation, for all functions
 */

class lfp_recorder : public Archiving_Node
{

public:
  lfp_recorder();
  lfp_recorder( const lfp_recorder& );
  virtual ~lfp_recorder();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& );
  void handle( DataLoggingRequest& ); // TODO: trenger vi denne?

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();
  void update( Time const&, const long, const long );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< lfp_recorder >;
  friend class UniversalDataLogger< lfp_recorder >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    std::vector< double > tau_rise;  //!< Rise time of synaptic conductance
                                     //!< in ms.
    std::vector< double > tau_decay; //!< Decay time of synaptic conductance
                                     //!< in ms.

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dictionary

    //! Return the number of receptor ports
    size_t
    n_receptors() const
    {
      return tau_rise.size();
    }
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor and assignment operator required because
   *       of C-style arrays.
   */
  struct State_
  {

    /**
     * Enumeration identifying elements in state vector State_::y_.
     * This enum identifies the elements of the vector. It must be public to be
     * accessible from the iteration function. The last two elements of this
     * enum (DG, G) will be repeated
     * n times at the end of the state vector State_::y with n being the number
     * of synapses.
     */
    enum StateVecElems
    {
      DG, // 1
      G,  // 2
      STATE_VECTOR_MIN_SIZE
    };

    static const size_t NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR = 2; // DG, G

    std::vector< double > y_; //!< neuron state

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum& );

  }; // State_

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( lfp_recorder& );
    Buffers_( const Buffers_&, lfp_recorder& );

    //! Logger for all analog data
    UniversalDataLogger< lfp_recorder > logger_;

    /** buffers and sums up incoming spikes */
    std::vector< RingBuffer > spikes_;

  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    std::vector< double > normalizer_;

    std::vector< double > P11_syn_;
    std::vector< double > P21_syn_;
    std::vector< double > P22_syn_;

    unsigned int receptor_types_size_;

  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger and RecordablesMap
  template < State_::StateVecElems elem >
  double
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // Data members -----------------------------------------------------------

  /**
   * @defgroup aeif_cond_beta_multisynapse
   * Instances of private data structures for the different types
   * of data pertaining to the model.
   * @note The order of definitions is important for speed.
   * @{
   */
  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;
  /** @} */

  //! Mapping of recordables names to access functions
  static RecordablesMap< lfp_recorder > recordablesMap_;
};

inline port
lfp_recorder::send_test_event( Node& target,
  rport receptor_type,
  synindex,
  bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
lfp_recorder::handles_test_event( DataLoggingRequest& dlr,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
lfp_recorder::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
lfp_recorder::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d );         // throws if BadProperty

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // LFP_RECORDER_H //
