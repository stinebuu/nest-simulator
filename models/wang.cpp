/*
 *  wang.cpp
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

#include "wang.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "dictdatum.h"
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap< nest::wang > nest::wang::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap< wang >::create()
  {
    // add state variables to recordables map
    insert_( names::V_m, &wang::get_ode_state_elem_< wang::State_::V_m > );

    //TODO SHOULD THESE BE NAMES???

    insert_( "G_AMPA", &wang::get_ode_state_elem_< wang::State_::G_AMPA > );
    insert_( "G_GABA", &wang::get_ode_state_elem_< wang::State_::G_GABA > );
    insert_( "G_NMDA", &wang::get_ode_state_elem_< wang::State_::G_NMDA > );
  }
}
// ---------------------------------------------------------------------------
//  Default constructors defining default parameters and state
// ---------------------------------------------------------------------------

nest::wang::Parameters_::Parameters_()
  : E_L( -70.0 )          // as mV
  , E_ex( 0.0 )           // as mV
  , E_in(-70.0)           // as mV
  , V_th(-50.0)           // as mV
  , V_reset(-55.0)        // as mV
  , C_m( 0.5 )            // as nF
  , g_L( 25.0 )           // as nS
  , t_ref( 2.0 )          // as ms
  , tau_AMPA( 2.0 )       // as ms
  , tau_GABA( 5.0 )       // as ms
  , tau_rise_NMDA( 2.0 )  // as ms
  , tau_decay_NMDA( 100 ) // as ms
  , alpha( 0.5 )          // as 1 / ms
  , g_ext_AMPA( 2.1 )     // as ns
  , g_rec_AMPA( 0.05 )    // as ns
  , g_GABA( 1.3 )         // as ns
  , g_NMDA( 0.165 )       // as ns
  , conc_Mg2( 1 )         // as mM
  , gsl_error_tol( 1e-3 )
{
}

nest::wang::State_::State_(  const Parameters_& p )
  : r_( 0 )
{
  // initial values for state variables
  ode_state_[ V_m ] = p.E_L; // as mV
  ode_state_[ G_NMDA ] = 0.0; // as real
  ode_state_[ G_AMPA ] = 0.0; // as real
  ode_state_[ G_GABA ] = 0.0; // as real
}

nest::wang::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    ode_state_[ i ] = s.ode_state_[ i ];
  }
}

nest::wang::State_& nest::wang::State_::operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    ode_state_[ i ] = s.ode_state_[ i ];
  }
  r_ = s.r_;
  return *this;
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

void
nest::wang::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::E_ex, E_ex );
  def< double >( d, names::E_in, E_in );
  def< double >( d, names::V_th, V_th );
  def< double >( d, names::V_reset, V_reset );
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::t_ref, t_ref );
  def< double >( d, names::tau_AMPA, tau_AMPA );
  def< double >( d, names::tau_GABA, tau_GABA );
  def< double >( d, names::tau_rise_NMDA, tau_rise_NMDA );
  def< double >( d, names::tau_decay_NMDA, tau_decay_NMDA );
  def< double >( d, names::alpha, alpha );
  def< double >( d, names::g_ext_AMPA, g_ext_AMPA );
  def< double >( d, names::g_rec_AMPA, g_rec_AMPA );
  def< double >( d, names::g_GABA, g_GABA );
  def< double >( d, names::g_NMDA, g_NMDA );
  def< double >( d, names::conc_Mg2, conc_Mg2 );
  def< double >( d, names::gsl_error_tol, gsl_error_tol );
}

void
nest::wang::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  // allow setting the membrane potential
  updateValueParam< double >( d, names::V_th, V_th, node );
  updateValueParam< double >( d, names::V_reset, V_reset, node );
  updateValueParam< double >( d, names::t_ref, t_ref, node );
  updateValueParam< double >( d, names::E_L, E_L, node );

  updateValueParam< double >( d, names::E_ex, E_ex, node );
  updateValueParam< double >( d, names::E_in, E_in, node );

  updateValueParam< double >( d, names::C_m, C_m, node );
  updateValueParam< double >( d, names::g_L, g_L, node );

  updateValueParam< double >( d, names::tau_AMPA, tau_AMPA, node );
  updateValueParam< double >( d, names::tau_GABA, tau_GABA, node );
  updateValueParam< double >( d, names::tau_rise_NMDA, tau_rise_NMDA, node );
  updateValueParam< double >( d, names::tau_decay_NMDA, tau_decay_NMDA, node );

  updateValueParam< double >( d, names::alpha, alpha, node );
  updateValueParam< double >( d, names::g_ext_AMPA, g_ext_AMPA, node );
  updateValueParam< double >( d, names::g_rec_AMPA, g_rec_AMPA, node );
  updateValueParam< double >( d, names::g_GABA, g_GABA, node );
  updateValueParam< double >( d, names::g_NMDA, g_NMDA, node );
  updateValueParam< double >( d, names::conc_Mg2, conc_Mg2, node );
  updateValueParam< double >( d, names::gsl_error_tol, gsl_error_tol, node );

  if ( V_reset >= V_th )
  {
    throw BadProperty( "Reset potential must be smaller than threshold." );
  }
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  if ( t_ref < 0 )
  {
    throw BadProperty( "Refractory time cannot be negative." );
  }
  if ( tau_AMPA <= 0 or tau_GABA <= 0 or tau_rise_NMDA <= 0 or tau_decay_NMDA <= 0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }
  if ( alpha <= 0 )
  {
    throw BadProperty( "alpha > 0 required." );
  }
  if ( g_ext_AMPA <= 0 or g_rec_AMPA <= 0 or g_GABA <= 0 or g_NMDA <= 0 )
  {
    throw BadProperty( "All synaptic conductances must be strictly positive." );
  }
  if ( conc_Mg2 <= 0 )
  {
    throw BadProperty( "Mg2 concentration must be strictly positive." );
  }
  if ( gsl_error_tol <= 0.0 )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

void
nest::wang::State_::get( DictionaryDatum& d ) const
{
  //TODO: SHOULD THESE BE NAMES????
  def< double >( d, names::V_m, ode_state_[ V_m ] ); // Membrane potential
  def< double >( d, "G_NMDA", ode_state_[ G_NMDA ] );
  def< double >( d, "G_AMPA", ode_state_[ G_AMPA ] );
  def< double >( d, "G_GABA", ode_state_[ G_GABA ] );
}

void
nest::wang::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
{
  updateValueParam< double >( d, names::V_m, ode_state_[ V_m ], node );
  updateValueParam< double >( d, "G_NMDA", ode_state_[ G_NMDA ], node );
  updateValueParam< double >( d, "G_AMPA", ode_state_[ G_AMPA ], node );
  updateValueParam< double >( d, "G_GABA", ode_state_[ G_GABA ], node );
}


nest::wang::Buffers_::Buffers_( wang &n )
  : logger_( n )
  //, spike_inputs_( std::vector< nest::RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
  , step_( Time::get_resolution().get_ms() )
  , integration_step_( step_ )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

nest::wang::Buffers_::Buffers_( const Buffers_ &, wang &n )
  : logger_( n )
  //, spike_inputs_( std::vector< nest::RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
  , step_( Time::get_resolution().get_ms() )
  , integration_step_( step_ )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

nest::wang::wang()
  :ArchivingNode()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();

  calibrate();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

nest::wang::wang( const wang& n_ )
  : ArchivingNode( n_ )
  , P_( n_.P_ )
  , S_( n_.S_ )
  , B_( n_.B_, *this )
{
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

nest::wang::~wang()
{
  // GSL structs may not have been allocated, so we need to protect destruction

  if ( B_.s_ )
  {
    gsl_odeiv_step_free( B_.s_ );
  }

  if ( B_.c_ )
  {
    gsl_odeiv_control_free( B_.c_ );
  }

  if ( B_.e_ )
  {
    gsl_odeiv_evolve_free( B_.e_ );
  }
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------

void
nest::wang::init_state_( const Node& proto )
{
  const wang& pr = downcast< wang >( proto );
  S_ = pr.S_;
}

void
nest::wang::init_buffers_()
{
  get_AMPA_spikes().clear(); //includes resize
  get_GABA_spikes().clear(); //includes resize

  //TODO: denne blir feil.
  //get_NMDA_spikes().clear(); //includes resize
  get_NMDA_spikes().clear(); //includes resize

  B_.logger_.reset(); // includes resize
  ArchivingNode::clear_history();

  if ( B_.s_ == 0 )
  {
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( B_.c_ == 0 )
  {
    B_.c_ = gsl_odeiv_control_y_new( P_.gsl_error_tol, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, P_.gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = wang_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );
  B_.step_ = nest::Time::get_resolution().get_ms();
  B_.integration_step_ = nest::Time::get_resolution().get_ms();
}

void
nest::wang::calibrate()
{
  B_.logger_.init();

  // internals V_
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps();
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int
nest::wang_dynamics(double, const double ode_state[], double f[], void* pnode)
{
  // a shorthand
  typedef nest::wang::State_ State_;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::wang& node = *( reinterpret_cast< nest::wang* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].

  const double I_ext_AMPA =
    node.P_.g_ext_AMPA * ( ode_state[ State_::V_m ] - node.P_.E_ex ) * ode_state[ State_::G_AMPA ];

  // MÅ GJØRE NOE MED VEKT, tror kanskje det er i update?
  const double I_rec_AMPA =
    node.P_.g_rec_AMPA * ( ode_state[ State_::V_m ] - node.P_.E_ex ) * ode_state[ State_::G_AMPA ];

  const double I_rec_GABA =
    node.P_.g_GABA * ( ode_state[ State_::V_m ] - node.P_.E_in ) * ode_state[ State_::G_GABA ];

  // TRENGER VEKT, OG SUMMERING(?)
  const double I_rec_NMDA =
      node.P_.g_NMDA * ( ode_state[ State_::V_m ] - node.P_.E_ex ) / ( 1 + node.P_.conc_Mg2 * std::exp( -0.062 * ode_state[ State_::V_m ] / 3.57 ) ) * ode_state[ State_::G_NMDA ];

  const double I_syn = I_ext_AMPA + I_rec_AMPA + I_rec_GABA + I_rec_NMDA;

  f[ State_::V_m ] = ( -node.P_.g_L * ( ode_state[ State_::V_m ] - node.P_.E_L ) - I_syn ) / node.P_.C_m;

  f[ State_::G_AMPA ] = -ode_state[ State_::G_AMPA ] / node.P_.tau_AMPA;
  f[ State_::G_GABA ] = -ode_state[ State_::G_GABA ] / node.P_.tau_GABA;
  f[ State_::G_NMDA ] = -ode_state[ State_::G_NMDA ] / node.P_.tau_decay_NMDA + node.P_.alpha * ode_state[ State_::x_NMDA ] * ( 1 - ode_state[ State_::G_NMDA ] );
  f[ State_::x_NMDA ] = -ode_state[ State_::x_NMDA ] / node.P_.tau_rise_NMDA;

  return GSL_SUCCESS;
}

void
nest::wang::update(nest::Time const & origin,const long from, const long to)
{
  for ( long lag = from ; lag < to ; ++lag )
  {
    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
                                                 B_.c_,
                                                 B_.s_,
                                                 &B_.sys_,              // system of ODE
                                                 &t,                    // from t
                                                 B_.step_,              // to t <= step
                                                 &B_.integration_step_, // integration step size
                                                 S_.ode_state_ );       // neuronal state

      if ( status != GSL_SUCCESS )
      {
        throw nest::GSLSolverFailure( get_name(), status );
      }
    }

    /* replace analytically solvable variables with precisely integrated values  */
    S_.ode_state_[ State_::G_AMPA ] += get_AMPA_spikes().get_value( lag );
    S_.ode_state_[ State_::G_GABA ] += get_GABA_spikes().get_value( lag );
    S_.ode_state_[ State_::x_NMDA ] += get_NMDA_spikes().get_value( lag );

    // absolute refractory period
    if ( S_.r_ != 0 )
    { // neuron is absolute refractory
      --S_.r_;
      S_.ode_state_[ State_::V_m ] = P_.V_reset;
    }
    else if ( S_.ode_state_[ State_::V_m ] >= P_.V_th )
    { // neuron is not absolute refractory
      S_.r_ = V_.RefractoryCounts;
      S_.ode_state_[ State_::V_m ] = P_.V_reset;

      set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ) );

      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send( *this, se, lag );
    }

    // voltage logging
    B_.logger_.record_data( origin.get_steps() + lag );
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void
nest::wang::handle( nest::DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
nest::wang::handle( nest::SpikeEvent &e )
{
  assert( e.get_delay_steps() > 0 );
  //assert( e.get_rport() < static_cast< int >( B_.spike_inputs_.size() ) );
  double steps = e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() );
  double weight_mul = e.get_weight() * e.get_multiplicity();

  if ( e.get_rport() == 0 )
  {
    B_.spike_AMPA_.add_value( steps, weight_mul );
  }
  else if ( e.get_rport() == 1 )
  {
    B_.spike_GABA_.add_value( steps, weight_mul );
  }
  else if ( e.get_rport() == 2 )
  {
    B_.spike_array_NMDA_= B_.spike_array_NMDA_.apply( []( RingBuffer rb ){ return rb.add_value( steps, weight_mul ); } );
  }
  else
  {
    throw BadProperty( "Specified receptor not defined." );
  }


  //B_.spike_inputs_[ e.get_rport() ].add_value(
  //  e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
  //  e.get_weight() * e.get_multiplicity() );
}
