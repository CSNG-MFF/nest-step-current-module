/*
 *  iaf_cond_exp_sc_nc.cpp
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

#include "iaf_cond_exp_sc_nc.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "compose.hpp"
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "booldatum.h"
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< stepcurrentmodule::iaf_cond_exp_sc_nc > stepcurrentmodule::iaf_cond_exp_sc_nc::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< stepcurrentmodule::iaf_cond_exp_sc_nc >::create()
{
  // use standard names whereever you can for consistency!
  insert_( nest::names::V_m, &stepcurrentmodule::iaf_cond_exp_sc_nc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc_nc::State_::V_M > );
  insert_( nest::names::g_ex, &stepcurrentmodule::iaf_cond_exp_sc_nc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc_nc::State_::G_EXC > );
  insert_( nest::names::g_in, &stepcurrentmodule::iaf_cond_exp_sc_nc::get_y_elem_< stepcurrentmodule::iaf_cond_exp_sc_nc::State_::G_INH > );
}
}

extern "C" inline int
stepcurrentmodule::iaf_cond_exp_sc_nc_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef stepcurrentmodule::iaf_cond_exp_sc_nc::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const stepcurrentmodule::iaf_cond_exp_sc_nc& node = *( reinterpret_cast< stepcurrentmodule::iaf_cond_exp_sc_nc* >( pnode ) );

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
  f[ 0 ] = is_refractory ? 0.0 : ( -I_L + node.B_.I_stim_ + node.B_.I_step_amp_ + node.B_.I_noise_amp_ + node.P_.I_e - I_syn_exc - I_syn_inh ) / node.P_.C_m; // ADDING NOISE

  f[ 1 ] = -y[ S::G_EXC ] / node.P_.tau_synE;
  f[ 2 ] = -y[ S::G_INH ] / node.P_.tau_synI;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

stepcurrentmodule::iaf_cond_exp_sc_nc::Parameters_::Parameters_()
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
  , amp_time_stamps_()
  , amp_values_()     // pA
  , allow_offgrid_amp_times_( false )
  , noise_mean_( 0.0 ) // pA NOISE
  , noise_std_( 0.0 )     // pA NOISE
  , noise_std_mod_( 0.0 ) // pA NOISE
  , noise_freq_( 0.0 )    // Hz NOISE
  , noise_phi_deg_( 0.0 )  // degree NOISE
  , noise_dt_( get_default_dt() ) // NOISE
{
}

stepcurrentmodule::iaf_cond_exp_sc_nc::State_::State_( const Parameters_& p )
  : r_( 0 )
  , y_0_( 0.0 )   // pA NOISE
  , y_1_( 0.0 )   // pA NOISE
{
  y_[ V_M ] = p.E_L;
  y_[ G_EXC ] = y_[ G_INH ] = 0;
}

stepcurrentmodule::iaf_cond_exp_sc_nc::State_::State_( const State_& s )
  : r_( s.r_ )
  , y_0_( s.y_0_ )   // pA NOISE
  , y_1_( s.y_1_ )   // pA NOISE
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

stepcurrentmodule::iaf_cond_exp_sc_nc::State_&
stepcurrentmodule::iaf_cond_exp_sc_nc::State_::operator=( const State_& s )
{
  r_ = s.r_;
  y_0_ = s.y_0_; // NOISE
  y_1_ = s.y_1_; // NOISE
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
stepcurrentmodule::iaf_cond_exp_sc_nc::Parameters_::get( DictionaryDatum& d ) const
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

  std::vector< double >* times_ms = new std::vector< double >();
  times_ms->reserve( amp_time_stamps_.size() );
  for ( auto amp_time_stamp : amp_time_stamps_ )
  {
    times_ms->push_back( amp_time_stamp.get_ms() );
  }
  ( *d )[ nest::names::amplitude_times ] = DoubleVectorDatum( times_ms );
  ( *d )[ nest::names::amplitude_values ] = DoubleVectorDatum( new std::vector< double >( amp_values_ ) );
  ( *d )[ nest::names::allow_offgrid_times ] = BoolDatum( allow_offgrid_amp_times_ );
 
  ( *d )[ nest::names::mean ] = noise_mean_; // NOISE
  ( *d )[ nest::names::std ] = noise_std_; // NOISE
  ( *d )[ nest::names::std_mod ] = noise_std_mod_; // NOISE
  ( *d )[ nest::names::dt ] = noise_dt_.get_ms(); // NOISE
  ( *d )[ nest::names::phase ] = noise_phi_deg_; // NOISE
  ( *d )[ nest::names::frequency ] = noise_freq_; // NOISE
}

nest::Time
stepcurrentmodule::iaf_cond_exp_sc_nc::Parameters_::validate_time_( double t, const Time& t_previous )
{
  // Force the amplitude change time to the grid
  // First, convert the time to tics, may not be on grid
  Time t_amp = Time::ms( t );
  if ( not t_amp.is_grid_time() )
  {
    if ( allow_offgrid_amp_times_ )
    {
      // In this case, we need to round to the end of the step
      // in which t lies, ms_stamp does that for us.
      t_amp = Time::ms_stamp( t );
    }
    else
    {
      std::stringstream msg;
      msg << "iaf_cond_exp_sc_nc: Time point " << t << " is not representable in current resolution.";
      throw BadProperty( msg.str() );
    }
  }

  const Time now = kernel().simulation_manager.get_time();
  if ( t_amp < now )
  {
    throw BadProperty(
		      String::compose("Amplitude can only change for the future (t >= %1).", now.get_ms() ) );
  }

  // t_amp is now the correct time stamp given the chosen options
  if ( t_amp <= t_previous )
  {
    throw BadProperty(
      "step_current_generator: amplitude "
      "times must be at strictly increasing "
      "time steps." );
  }

  // when we get here, we know that the spike time is valid
  return t_amp;
}


void
stepcurrentmodule::iaf_cond_exp_sc_nc::Parameters_::set( const DictionaryDatum& d, Buffers_& b, Node* node )
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

  std::vector< double > new_times;
  const bool times_changed = updateValue< std::vector< double > >( d, nest::names::amplitude_times, new_times );
  const bool values_changed = updateValue< std::vector< double > >( d, nest::names::amplitude_values, amp_values_ );
  const bool allow_offgrid_changed = updateValue< bool >( d, nest::names::allow_offgrid_times, allow_offgrid_amp_times_ );

  if ( times_changed xor values_changed )
  {
	throw BadProperty( "Amplitude times and values must be reset together." );
  }

  if ( allow_offgrid_changed and not( times_changed or amp_time_stamps_.empty() ) )
  {
	// times_changed implies values_changed
	throw BadProperty(
	  "allow_offgrid_times can only be changed before "
      "amplitude_times have been set, or together with "
      "amplitude_times and amplitude_values." );
  }

  const size_t times_size = times_changed ? new_times.size() : amp_time_stamps_.size();

  if ( times_size != amp_values_.size() )
  {
    throw BadProperty( "Amplitude times and values have to be the same size." );
  }

  if ( times_changed )
  {
	std::vector< Time > new_stamps;
	new_stamps.reserve( times_size );

	if ( not new_times.empty() )
	{
      // insert first change, we are sure we have one
      new_stamps.push_back( validate_time_( new_times[ 0 ], Time::neg_inf() ) );

      // insert all others
      for ( size_t idx = 1; idx < times_size; ++idx )
      {
    	new_stamps.push_back( validate_time_( new_times[ idx ], new_stamps[ idx - 1 ] ) );
      }
    }

    // if we get here, all times have been successfully converted
    amp_time_stamps_.swap( new_stamps );
  }

  if ( times_changed or values_changed )
  {
    b.I_step_idx_ = 0; // reset if we got new data
  }

  nest::updateValueParam< double >( d, nest::names::mean, noise_mean_, node ); // NOISE
  nest::updateValueParam< double >( d, nest::names::std, noise_std_, node ); // NOISE
  nest::updateValueParam< double >( d, nest::names::std_mod, noise_std_mod_, node ); // NOISE
  nest::updateValueParam< double >( d, nest::names::frequency, noise_freq_, node ); // NOISE
  nest::updateValueParam< double >( d, nest::names::phase, noise_phi_deg_, node ); // NOISE
  double noise_dt; // NOISE

  // NOISE
  if ( nest::updateValueParam< double >( d, nest::names::dt, noise_dt, node ) )
  {
    noise_dt_ = Time::ms( noise_dt );
  }
  if ( noise_std_ < 0 )
  {
    throw BadProperty( "The standard deviation cannot be negative." );
  }
  if ( noise_std_mod_ < 0 )
  {
    throw BadProperty( "The standard deviation cannot be negative." );
  }
  if ( noise_std_mod_ > noise_std_ )
  {
    throw BadProperty(
      "The modulation apmlitude must be smaller or equal to the baseline "
      "amplitude." );
  }

  if ( not noise_dt_.is_step() )
  {
    throw BadProperty("The timestep should be a multiple of the resolution.");
  }
}

void
stepcurrentmodule::iaf_cond_exp_sc_nc::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::V_m, y_[ V_M ] ); // Membrane potential
  def< double >( d, nest::names::g_ex, y_[ G_EXC ] );
  def< double >( d, nest::names::g_in, y_[ G_INH ] );

  def< double >( d, nest::names::g_ex, y_[ G_EXC ] );
  def< double >( d, nest::names::g_in, y_[ G_INH ] );
  ( *d )[ nest::names::y_0 ] = y_0_; // NOISE
  ( *d )[ nest::names::y_1 ] = y_1_; // NOISE
}

void
stepcurrentmodule::iaf_cond_exp_sc_nc::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
{
  nest::updateValueParam< double >( d, nest::names::V_m, y_[ V_M ], node );
  nest::updateValueParam< double >( d, nest::names::g_ex, y_[ G_EXC ], node );
  nest::updateValueParam< double >( d, nest::names::g_in, y_[ G_INH ], node );
}

stepcurrentmodule::iaf_cond_exp_sc_nc::Buffers_::Buffers_( iaf_cond_exp_sc_nc& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
  , next_step_( 0 ) // NOISE
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

stepcurrentmodule::iaf_cond_exp_sc_nc::Buffers_::Buffers_( const Buffers_& b, iaf_cond_exp_sc_nc& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
  , next_step_( b.next_step_ ) // NOISE
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

stepcurrentmodule::iaf_cond_exp_sc_nc::iaf_cond_exp_sc_nc()
  : ArchivingNode()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

stepcurrentmodule::iaf_cond_exp_sc_nc::iaf_cond_exp_sc_nc( const iaf_cond_exp_sc_nc& n )
  : ArchivingNode( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

stepcurrentmodule::iaf_cond_exp_sc_nc::~iaf_cond_exp_sc_nc()
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
stepcurrentmodule::iaf_cond_exp_sc_nc::init_buffers_()
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

  B_.sys_.function = iaf_cond_exp_sc_nc_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );


  B_.I_step_idx_ = 0;
  B_.I_step_amp_ = 0.0;
  B_.I_stim_ = 0.0;
 
  B_.next_step_ = 0; // NOISE
  B_.I_noise_amp_ = 0.0; // NOISE
}

void
stepcurrentmodule::iaf_cond_exp_sc_nc::pre_run_hook()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );


  V_.dt_steps_ = P_.noise_dt_.get_steps(); // NOISE

  const double h = Time::get_resolution().get_ms(); // NOISE
  const double t = kernel().simulation_manager.get_time().get_ms(); // NOISE

  // scale Hz to ms NOISE
  const double omega = 2.0 * numerics::pi * P_.noise_freq_ / 1000.0;
  const double phi_rad = P_.noise_phi_deg_ * 2.0 * numerics::pi / 360.0;

  // initial state NOISE
  S_.y_0_ = std::cos( omega * t + phi_rad );
  S_.y_1_ = std::sin( omega * t + phi_rad );

  // matrix elements NOISE
  V_.A_00_ = std::cos( omega * h );
  V_.A_01_ = -std::sin( omega * h );
  V_.A_10_ = std::sin( omega * h );
  V_.A_11_ = std::cos( omega * h );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
stepcurrentmodule::iaf_cond_exp_sc_nc::update( Time const& origin, const long from, const long to )
{

  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const long start = origin.get_steps();

  for ( long lag = from; lag < to; ++lag )
  {

	if ( B_.I_step_idx_ < P_.amp_time_stamps_.size()
		 and origin.get_steps() + lag == P_.amp_time_stamps_[ B_.I_step_idx_ ].get_steps() )
	{
	  B_.I_step_amp_ = P_.amp_values_[ B_.I_step_idx_ ];
	  B_.I_step_idx_++;
	}

    double t = 0.0;
    const long now = start + lag; // NOISE

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

    // NOISE
    if ( P_.noise_std_mod_ != 0. )
    {
      const double y_0 = S_.y_0_;
      S_.y_0_ = V_.A_00_ * y_0 + V_.A_01_ * S_.y_1_;
      S_.y_1_ = V_.A_10_ * y_0 + V_.A_11_ * S_.y_1_;
    }

    // >= in case we woke from inactivity NOISE
    if ( now >= B_.next_step_ )
    {
      // compute new noise current
      B_.I_noise_amp_ = P_.noise_mean_
          + std::sqrt( P_.noise_std_ * P_.noise_std_ + S_.y_1_ * P_.noise_std_mod_ * P_.noise_std_mod_ )
            * V_.normal_dist_( get_vp_specific_rng( get_thread() ) );
      
      // use now as reference, in case we woke up from inactive period
      B_.next_step_ = now + V_.dt_steps_;
    }

  }
}

void
stepcurrentmodule::iaf_cond_exp_sc_nc::handle( SpikeEvent& e )
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
stepcurrentmodule::iaf_cond_exp_sc_nc::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
}

void
stepcurrentmodule::iaf_cond_exp_sc_nc::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
