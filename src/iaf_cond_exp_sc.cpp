/*
 *  iaf_cond_exp_sc.cpp
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

#include "iaf_cond_exp_sc.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< stepcurrentmodule::iaf_cond_exp_sc > stepcurrentmodule::iaf_cond_exp_sc::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< stepcurrentmodule::iaf_cond_exp_sc >::create()
{
  // use standard names whereever you can for consistency!
  insert_( nest::names::V_m, &stepcurrentmodule::iaf_cond_exp_sc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc::State_::V_M > );
  insert_( nest::names::g_ex, &stepcurrentmodule::iaf_cond_exp_sc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc::State_::G_EXC > );
  insert_( nest::names::g_in, &stepcurrentmodule::iaf_cond_exp_sc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc::State_::G_INH > );
}
}

extern "C" inline int
stepcurrentmodule::iaf_cond_exp_sc_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef stepcurrentmodule::iaf_cond_exp_sc::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const stepcurrentmodule::iaf_cond_exp_sc& node = *( reinterpret_cast< stepcurrentmodule::iaf_cond_exp_sc* >( pnode ) );

  const bool is_refractory = node.S_.r_ > 0;

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // Clamp membrane potential to V_reset while refractory, otherwise bound
  // it to V_th.
  const double V = is_refractory ? node.P_.V_reset_ : std::min( y[ S::V_M ], node.P_.V_th_ );

  const double I_syn_exc = y[ S::G_EXC ] * ( V - node.P_.E_ex );
  const double I_syn_inh = y[ S::G_INH ] * ( V - node.P_.E_in );
  const double I_L = node.P_.g_L * ( V - node.P_.E_L );

  // V dot
  f[ 0 ] = is_refractory ? 0.0 : ( -I_L + node.B_.I_stim_ + node.P_.I_e - I_syn_exc - I_syn_inh ) / node.P_.C_m;

  f[ 1 ] = -y[ S::G_EXC ] / node.P_.tau_synE;
  f[ 2 ] = -y[ S::G_INH ] / node.P_.tau_synI;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

stepcurrentmodule::iaf_cond_exp_sc::Parameters_::Parameters_()
  : V_th_( -55.0 )    // mV
  , V_reset_( -60.0 ) // mV
  , t_ref_( 2.0 )     // ms
  , g_L( 16.6667 )    // nS
  , C_m( 250.0 )      // pF
  , E_ex( 0.0 )       // mV
  , E_in( -85.0 )     // mV
  , E_L( -70.0 )      // mV
  , tau_synE( 0.2 )   // ms
  , tau_synI( 2.0 )   // ms
  , I_e( 0.0 )        // pA
{
}

stepcurrentmodule::iaf_cond_exp_sc::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ V_M ] = p.E_L;
  y_[ G_EXC ] = y_[ G_INH ] = 0;
}

stepcurrentmodule::iaf_cond_exp_sc::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

stepcurrentmodule::iaf_cond_exp_sc::State_&
stepcurrentmodule::iaf_cond_exp_sc::State_::operator=( const State_& s )
{
  r_ = s.r_;
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
  return *this;
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
stepcurrentmodule::iaf_cond_exp_sc::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::V_th, V_th_ );
  def< double >( d, nest::names::V_reset, V_reset_ );
  def< double >( d, nest::names::t_ref, t_ref_ );
  def< double >( d, nest::names::g_L, g_L );
  def< double >( d, nest::names::E_L, E_L );
  def< double >( d, nest::names::E_ex, E_ex );
  def< double >( d, nest::names::E_in, E_in );
  def< double >( d, nest::names::C_m, C_m );
  def< double >( d, nest::names::tau_syn_ex, tau_synE );
  def< double >( d, nest::names::tau_syn_in, tau_synI );
  def< double >( d, nest::names::I_e, I_e );
}

void
stepcurrentmodule::iaf_cond_exp_sc::Parameters_::set( const DictionaryDatum& d, Node* node )
{
  // allow setting the membrane potential
  nest::updateValueParam< double >( d, nest::names::V_th, V_th_, node );
  nest::updateValueParam< double >( d, nest::names::V_reset, V_reset_, node );
  nest::updateValueParam< double >( d, nest::names::t_ref, t_ref_, node );
  nest::updateValueParam< double >( d, nest::names::E_L, E_L, node );

  nest::updateValueParam< double >( d, nest::names::E_ex, E_ex, node );
  nest::updateValueParam< double >( d, nest::names::E_in, E_in, node );

  nest::updateValueParam< double >( d, nest::names::C_m, C_m, node );
  nest::updateValueParam< double >( d, nest::names::g_L, g_L, node );

  nest::updateValueParam< double >( d, nest::names::tau_syn_ex, tau_synE, node );
  nest::updateValueParam< double >( d, nest::names::tau_syn_in, tau_synI, node );

  nest::updateValueParam< double >( d, nest::names::I_e, I_e, node );
  if ( V_reset_ >= V_th_ )
  {
    throw BadProperty( "Reset potential must be smaller than threshold." );
  }
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  if ( t_ref_ < 0 )
  {
    throw BadProperty( "Refractory time cannot be negative." );
  }
  if ( tau_synE <= 0 || tau_synI <= 0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }
}

void
stepcurrentmodule::iaf_cond_exp_sc::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::V_m, y_[ V_M ] ); // Membrane potential
  def< double >( d, nest::names::g_ex, y_[ G_EXC ] );
  def< double >( d, nest::names::g_in, y_[ G_INH ] );
}

void
stepcurrentmodule::iaf_cond_exp_sc::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
{
  nest::updateValueParam< double >( d, nest::names::V_m, y_[ V_M ], node );
  nest::updateValueParam< double >( d, nest::names::g_ex, y_[ G_EXC ], node );
  nest::updateValueParam< double >( d, nest::names::g_in, y_[ G_INH ], node );
}

stepcurrentmodule::iaf_cond_exp_sc::Buffers_::Buffers_( iaf_cond_exp_sc& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

stepcurrentmodule::iaf_cond_exp_sc::Buffers_::Buffers_( const Buffers_&, iaf_cond_exp_sc& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

stepcurrentmodule::iaf_cond_exp_sc::iaf_cond_exp_sc()
  : ArchivingNode()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

stepcurrentmodule::iaf_cond_exp_sc::iaf_cond_exp_sc( const iaf_cond_exp_sc& n )
  : ArchivingNode( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

stepcurrentmodule::iaf_cond_exp_sc::~iaf_cond_exp_sc()
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

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
stepcurrentmodule::iaf_cond_exp_sc::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  ArchivingNode::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

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
    B_.c_ = gsl_odeiv_control_y_new( 1e-3, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-3, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = iaf_cond_exp_sc_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
stepcurrentmodule::iaf_cond_exp_sc::pre_run_hook()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
stepcurrentmodule::iaf_cond_exp_sc::update( Time const& origin, const long from, const long to )
{

  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag )
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
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state
      if ( status != GSL_SUCCESS )
      {
        throw GSLSolverFailure( get_name(), status );
      }
    }

    S_.y_[ State_::G_EXC ] += B_.spike_exc_.get_value( lag );
    S_.y_[ State_::G_INH ] += B_.spike_inh_.get_value( lag );

    // absolute refractory period
    if ( S_.r_ )
    { // neuron is absolute refractory
      --S_.r_;
      S_.y_[ State_::V_M ] = P_.V_reset_;
    }
    else
      // neuron is not absolute refractory
      if ( S_.y_[ State_::V_M ] >= P_.V_th_ )
      {
        S_.r_ = V_.RefractoryCounts_;
        S_.y_[ State_::V_M ] = P_.V_reset_;

        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

        SpikeEvent se;
        kernel().event_delivery_manager.send( *this, se, lag );
      }

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
stepcurrentmodule::iaf_cond_exp_sc::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spike_inh_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() );
  }
}

void
stepcurrentmodule::iaf_cond_exp_sc::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
}

void
stepcurrentmodule::iaf_cond_exp_sc::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
