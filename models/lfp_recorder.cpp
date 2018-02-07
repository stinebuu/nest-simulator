/*
 *  lfp_recorder.cpp
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

#include "lfp_recorder.h"

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"


/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::lfp_recorder >
  nest::lfp_recorder::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< lfp_recorder >::create()
{
  // use standard names wherever you can for consistency!
  insert_( Name("lfp"), &lfp_recorder:: get_y_elem_< lfp_recorder::State_::G > );
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

lfp_recorder::Parameters_::Parameters_()
  : tau_rise( 1, 2.0 )   // ms
  , tau_decay( 1, 20.0 ) // ms
{
}

lfp_recorder::State_::State_( const Parameters_& p )
  : y_( STATE_VECTOR_MIN_SIZE, 0.0 )
{
}

lfp_recorder::State_::State_( const State_& s )
{
  y_ = s.y_;
}

lfp_recorder::State_& lfp_recorder::State_::
operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program

  y_ = s.y_;
  return *this;
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
lfp_recorder::Parameters_::get( DictionaryDatum& d ) const
{
  def< size_t >( d, names::n_receptors, n_receptors() );
  ArrayDatum tau_rise_ad( tau_rise );
  ArrayDatum tau_decay_ad( tau_decay );
  def< ArrayDatum >( d, names::tau_rise, tau_rise_ad );
  def< ArrayDatum >( d, names::tau_decay, tau_decay_ad );
}

void
lfp_recorder::Parameters_::set( const DictionaryDatum& d )
{
  const size_t old_n_receptors = n_receptors();
  bool taur_flag =
    updateValue< std::vector< double > >( d, names::tau_rise, tau_rise );
  bool taud_flag =
    updateValue< std::vector< double > >( d, names::tau_decay, tau_decay );
  if ( taur_flag || taud_flag )
  { // receptor arrays have been modified
    if ( ( tau_rise.size() != old_n_receptors
           || tau_decay.size() != old_n_receptors )
      && ( not taur_flag || not taud_flag ) )
    {
      throw BadProperty(
        "If the number of receptor ports is changed, the two arrays "
        "tau_rise and tau_decay must be provided." );
    }
    for ( size_t i = 0; i < tau_rise.size(); ++i )
    {
      if ( tau_rise[ i ] <= 0 || tau_decay[ i ] <= 0 )
      {
        throw BadProperty(
          "All synaptic time constants must be strictly positive" );
      }
      if ( tau_decay[ i ] < tau_rise[ i ] )
      {
        throw BadProperty(
          "Synaptic rise time must be smaller than or equal to decay time." ); // TODO: Denne slår ut noen ganger, ikke relevant for oss?
      }
    }
  }
}

void
lfp_recorder::State_::get( DictionaryDatum& d ) const
{
  std::vector< double >* dg = new std::vector< double >();
  std::vector< double >* g = new std::vector< double >();

  for ( size_t i = 0;
        i < ( y_.size() / State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR );
        ++i )
  {
    dg->push_back( y_[ State_::DG
      + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] );
    g->push_back( y_[ State_::G
      + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] );
  }

  ( *d )[ names::dg ] = DoubleVectorDatum( dg );
  ( *d )[ names::g ] = DoubleVectorDatum( g );
}

void
lfp_recorder::State_::set( const DictionaryDatum& d )
{
}

lfp_recorder::Buffers_::Buffers_(
  lfp_recorder& n )
  : logger_( n )
{
}

lfp_recorder::Buffers_::Buffers_( const Buffers_& b,
  lfp_recorder& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

lfp_recorder::lfp_recorder()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

lfp_recorder::lfp_recorder(
  const lfp_recorder& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

lfp_recorder::~lfp_recorder()
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
lfp_recorder::init_state_( const Node& proto )
{
  const lfp_recorder& pr =
    downcast< lfp_recorder >( proto );
  S_ = pr.S_;
}

void
lfp_recorder::init_buffers_()
{
  B_.spikes_.clear(); // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();
}

void
lfp_recorder::calibrate()
{
  // Ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  const double h = Time::get_resolution().get_ms();

  V_.P11_syn_.resize( P_.n_receptors() );
  V_.P21_syn_.resize( P_.n_receptors() );
  V_.P22_syn_.resize( P_.n_receptors() );

  S_.y_.resize( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * P_.n_receptors(), 0.0 );

  B_.spikes_.resize( P_.n_receptors() );

  // normalizer_ will be initialized in the loop below.
  V_.normalizer_.resize( P_.n_receptors() );

  for ( size_t i = 0; i < P_.n_receptors(); i++ )
  {
    // Set matrix components
    V_.P11_syn_[ i ] = std::exp( -h / P_.tau_decay[ i ] );
    V_.P22_syn_[ i ] = std::exp( -h / P_.tau_rise[ i ] );
    V_.P21_syn_[ i ] = ( ( P_.tau_decay[ i ] * P_.tau_rise[ i ] ) / ( P_.tau_decay[ i ] - P_.tau_rise[ i ] ) ) * ( V_.P11_syn_[ i ] - V_.P22_syn_[ i ] );

    // The denominator (denom1) that appears in the expression of the peak time
    // is computed here to check that it is != 0.
    // Another denominator denom2 appears in the expression of the
    // normalization factor normalizer_.
    // Both denom1 and denom2 are null if tau_decay = tau_rise, but they
    // can also be null if tau_decay and tau_rise are not equal but very
    // close to each other, due to the numerical precision limits.
    // In such case the beta function reduces to the alpha function,
    // and the normalization factor for the alpha function should be used.
    double denom1 = P_.tau_decay[ i ] - P_.tau_rise[ i ];
    double denom2 = 0;
    if ( denom1 != 0 )
    {
      // peak time
      const double t_p = P_.tau_decay[ i ] * P_.tau_rise[ i ]
        * std::log( P_.tau_decay[ i ] / P_.tau_rise[ i ] ) / denom1;

      // Another denominator is computed here to check that it is != 0
      denom2 = std::exp( -t_p / P_.tau_decay[ i ] )
        - std::exp( -t_p / P_.tau_rise[ i ] );
    }
    if ( denom2 == 0 ) // if rise time == decay time use alpha function
    {                  // use normalization for alpha function in this case
      V_.normalizer_[ i ] = 1.0 * numerics::e / P_.tau_decay[ i ];
    }
    else // if rise time != decay time use beta function
    {
      // normalization factor for beta function
      V_.normalizer_[ i ] = ( ( 1. / P_.tau_rise[ i ] ) - ( 1. / P_.tau_decay[ i ] ) ) / denom2;
    }

    B_.spikes_[ i ].resize();
  }
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */
void
lfp_recorder::update( Time const& origin,
  const long from,
  const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag ) // proceed by stepsize B_.step_
  {
    for ( size_t i = 0; i < P_.n_receptors(); ++i )
    {
      S_.y_[ State_::G + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] =
          V_.P21_syn_[ i ] * S_.y_[ State_::DG + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] +
          V_.P22_syn_[ i ] * S_.y_[ State_::G + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ];

      S_.y_[ State_::DG
        + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] *= V_.P11_syn_[ i ];

      S_.y_[ State_::DG
              + ( State_::NUMBER_OF_STATES_ELEMENTS_PER_RECEPTOR * i ) ] +=
                  V_.normalizer_[ i ] * B_.spikes_[ i ].get_value( lag );
    }

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );

  } // for-loop
}

port
lfp_recorder::handles_test_event( SpikeEvent&,
  rport receptor_type )
{
  if ( receptor_type <= 0
    || receptor_type > static_cast< port >( P_.n_receptors() ) )
  {
    throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
  }
  return receptor_type;
}

void
lfp_recorder::handle( SpikeEvent& e )
{
  if ( e.get_weight() < 0 )
  {
    throw BadProperty(
      "Synaptic weights for conductance-based multisynapse models "
      "must be positive." );
  }
  assert( e.get_delay() > 0 );
  assert(
    ( e.get_rport() > 0 ) && ( ( size_t ) e.get_rport() <= P_.n_receptors() ) );

  B_.spikes_[ e.get_rport() - 1 ].add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}

void
lfp_recorder::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

} // namespace nest

