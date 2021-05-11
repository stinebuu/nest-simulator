/**
 *  wang.h
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
**/
#ifndef WANG
#define WANG

#include "config.h"

#ifndef HAVE_GSL
#error "The GSL library is required for neurons that require a numerical solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// C++ includes
#include <valarray>

// Includes from sli:
//#include "dictdatum.h"

namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
**/
extern "C" inline int wang_dynamics( double, const double y[], double f[], void* pnode );


/* BeginDocumentation
  Name: wang.

  Description:



  Parameters:
  The following parameters can be set in the status dictionary.
E_L [mV]  Resting Potential
V_th [mV]  Threshold Potential
V_reset [mV]  Reset Potential
C_m [nF]  Membrane Capacitance
g_L [nS]  Leak Conductance
t_ref [ms]  Refractory period
conc_Mg2 [real]  can be wrong, should be mM
w [real]  This is wrong
 E_ex mV = 0 mV         # Excitatory reversal Potential
 E_in mV = -85.0 mV     # Inhibitory reversal Potential
 E_L mV = -70.0 mV      # Leak reversal Potential (aka resting potential)
 tau_syn_ex ms = 0.2 ms # Synaptic Time Constant Excitatory Synapse
 tau_syn_in ms = 2.0 ms # Synaptic Time Constant for Inhibitory Synapse


  Dynamic state variables:
r [integer]  counts number of tick during the refractory period
V_m [mV]  membrane potential


  Sends: SpikeEvent

  Receives: Spike,  DataLoggingRequest
*/
class wang : public ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  wang();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c calibrate().
  **/
  wang(const wang &);

  /**
   * Destructor.
  **/
  ~wang();

  // -------------------------------------------------------------------------
  //   Import sets of overloaded virtual functions.
  //   See: Technical Issues / Virtual Functions: Overriding, Overloading,
  //        and Hiding
  // -------------------------------------------------------------------------

  using Node::handles_test_event;
  using Node::handle;

  /**
   * Used to validate that we can send SpikeEvent to desired target:port.
  **/
  port send_test_event( Node& target, rport receptor_type, synindex, bool );

  // -------------------------------------------------------------------------
  //   Functions handling incoming events.
  //   We tell NEST that we can handle incoming events of various types by
  //   defining handle() for the given event.
  // -------------------------------------------------------------------------


  void handle( SpikeEvent & );        //! accept spikes
  void handle( DataLoggingRequest & );//! allow recording with multimeter
  port handles_test_event( SpikeEvent&, port );
  port handles_test_event( DataLoggingRequest&, port );

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:
  /**
   * Synapse types to connect to
   * @note Excluded upper and lower bounds are defined as INF_, SUP_.
   *       Excluding port 0 avoids accidental connections.
  **/
  enum SynapseTypes
  {
    INF_SPIKE_RECEPTOR = 0,
    AMPA,
    GABA,
    NMDA,
    SUP_SPIKE_RECEPTOR
  };

  /**
   * Reset state of neuron.
  **/
  void init_state_( const Node& proto );

  /**
   * Reset internal buffers of neuron.
  **/
  void init_buffers_();

  /**
   * Initialize auxiliary quantities, leave parameters and state untouched.
  **/
  void calibrate();

  /**
   * Take neuron through given time interval
  **/
  void update( Time const &, const long, const long );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< wang >;
  friend class UniversalDataLogger< wang >;

  /**
   * Free parameters of the neuron.
   *
   * These are the parameters that can be set by the user through @c `node.set()`.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update() and are not reset by
   * @c ResetNetwork.
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct Parameters_
  {
    double E_L;     //! Resting Potential
    double E_ex;    //! Excitatory resting Potential
    double E_in;    //! Inhibitory resting Potential
    double V_th;    //! Threshold Potential
    double V_reset; //! Reset Potential
    double C_m;     //! Membrane Capacitance
    double g_L;     //! Leak Conductance
    double t_ref;   //! Refractory period
    double tau_AMPA;
    double tau_GABA;
    double tau_rise_NMDA;
    double tau_decay_NMDA;
    double alpha;
    double g_ext_AMPA;
    double g_rec_AMPA;
    double g_GABA;
    double g_NMDA;
    //!  can be wrong, should be mM
    double conc_Mg2;

    double gsl_error_tol;

    /**
     * Initialize parameters to their default values.
    **/
    Parameters_();

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Node* node ); //!< Set values from dicitonary
  };

  /**
   * Dynamic state of the neuron.
   *
   *
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c `node.set()`. The state variables are initialized from the model
   * prototype when the node is created. State variables are reset by @c ResetNetwork.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct State_
  {
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      G_AMPA,
      G_GABA,
      G_NMDA,
      x_NMDA,
      STATE_VEC_SIZE
    };

    //! state vector, must be C-array for GSL solver
    double ode_state_[ STATE_VEC_SIZE ];
    int r_; //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  /**
   * Internal variables of the neuron.
   *
   *
   *
   * These variables must be initialized by @c calibrate, which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c calibrate(). If Variables_ has members that
   *       cannot destroy themselves, Variables_ will need a destructor.
  **/
  struct Variables_
  {
    //!  refractory time in steps
    long RefractoryCounts;
  };

  /**
   * Buffers of the neuron.
   * Usually buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c calibrate() on the first call to @c Simulate after the start of NEST,
   * ResetKernel or ResetNetwork.
   * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
   *       cannot destroy themselves, Buffers_ will need a destructor.
  **/
  struct Buffers_
  {
    Buffers_(wang &);
    Buffers_(const Buffers_ &, wang &);

    /**
     * Logger for all analog data
    **/
    UniversalDataLogger< wang > logger_;

    // -----------------------------------------------------------------------
    //   Buffers and sums of incoming spikes per timestep
    // -----------------------------------------------------------------------
    //std::vector< RingBuffer > spike_inputs_;
    RingBuffer spike_AMPA_;
    RingBuffer spike_GABA_;
    std::valarray< RingBuffer > spike_array_NMDA_;

    // -----------------------------------------------------------------------
    //   GSL ODE solver data structures
    // -----------------------------------------------------------------------

    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // integration_step_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double step_;             //!< step size in ms
    double integration_step_; //!< current integration time step, updated by GSL
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------

  inline RingBuffer& get_AMPA_spikes() { return B_.spike_AMPA_; };
  inline RingBuffer& get_GABA_spikes() { return B_.spike_GABA_; };
  inline RingBuffer& get_NMDA_spikes() { return B_.spike_array_NMDA_; };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double
  get_ode_state_elem_() const
  {
    return S_.ode_state_[ elem ];
  }

  // -------------------------------------------------------------------------
  //   Member variables of neuron model.
  //   Each model neuron should have precisely the following four data members,
  //   which are one instance each of the parameters, state, buffers and variables
  //   structures. Experience indicates that the state and variables member should
  //   be next to each other to achieve good efficiency (caching).
  //   Note: Devices require one additional data member, an instance of the
  //   ``Device`` child class they belong to.
  // -------------------------------------------------------------------------

  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static RecordablesMap< wang > recordablesMap_;
  friend int wang_dynamics( double, const double y[], double f[], void* pnode );

}; /* neuron wang */

inline port
wang::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
wang::handles_test_event( SpikeEvent&, port receptor_type )
{
  //assert( B_.spike_inputs_.size() == 3 );

  if ( !( INF_SPIKE_RECEPTOR < receptor_type && receptor_type < SUP_SPIKE_RECEPTOR ) )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
    return 0;
  }
  else
  {
    return receptor_type - 1;
  }
}

inline port
wang::handles_test_event( DataLoggingRequest& dlr, port receptor_type )
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void
wang::get_status( DictionaryDatum & d ) const
{
  P_.get( d );
  S_.get( d );
  ArchivingNode::get_status( d );

  DictionaryDatum receptor_type = new Dictionary();

  ( *receptor_type )[ names::AMPA ] = AMPA;
  ( *receptor_type )[ names::GABA ] = GABA;
  ( *receptor_type )[ names::NMDA ] = NMDA;

  ( *d )[ names::receptor_types ] = receptor_type;

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
wang::set_status( const DictionaryDatum & d )
{
  Parameters_ ptmp = P_;     // temporary copy in case of errors
  ptmp.set( d, this );       // throws if BadProperty
  State_ stmp = S_;          // temporary copy in case of errors
  stmp.set( d, ptmp, this ); // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
};
} // namespace

#endif // WANG
