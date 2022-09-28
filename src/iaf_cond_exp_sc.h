/*
 *  iaf_cond_exp_sc.h
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

#ifndef IAF_COND_EXP_SC_H
#define IAF_COND_EXP_SC_H

// Generated includes:
#include "config.h"

#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace stepcurrentmodule
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
 */
extern "C" int iaf_cond_exp_sc_dynamics( double, const double*, double*, void* );

/* BeginUserDocs: neuron, integrate-and-fire, conductance-based

Short description
+++++++++++++++++

Simple conductance based leaky integrate-and-fire neuron model

Description
+++++++++++

``iaf_cond_exp_sc`` is an implementation of a spiking neuron using integrate-and-fire
dynamics with conductance-based synapses and includes a built-in step current generator.
Incoming spike events induce a postsynaptic change of conductance modelled by an exponential
function. The exponential function is normalized such that an event of weight 1.0 results
in a peak conductance of 1 nS. See also [1]_.

Parameters
++++++++++

The following parameters can be set in the status dictionary.

==================== =====  =======================================================
 V_m                 mV      Membrane potential
 E_L                 mV      Leak reversal potential
 C_m                 pF      Capacity of the membrane
 t_ref               ms      Duration of refractory period
 V_th                mV      Spike threshold
 V_reset             mV      Reset potential of the membrane
 E_ex                mV      Excitatory reversal potential
 E_in                mV      Inhibitory reversal potential
 g_L                 nS      Leak conductance
 tau_syn_ex          ms      Exponential decay time constant of excitatory synaptic
                             conductance kernel
 tau_syn_in          ms      Exponential decay time constant of inhibitory synaptic
                             conductance kernel
 I_e                 pA      Constant input current

 amplitude_times     ms      Times at which current changes (list)
 amplitude_values    pA      Amplitudes of step current current (list)
 allow_offgrid_times         If True, allow off-grid times (default: False)
==================== ======  =======================================================

.. note::

	If ``allow_offgrid_spikes`` is set false, times will be rounded to the
	nearest step if they are less than tic/2 from the step, otherwise NEST
	reports an error. If true, times are rounded to the nearest step if
	within tic/2 from the step, otherwise they are rounded up to the *end*
	of the step.

	Times of amplitude changes must be strictly increasing after conversion
	to simulation time steps. The option ``allow_offgrid_times`` may be
	useful, for example, if you are using randomized times for current changes
	which typically would not fall onto simulation time steps.


Sends
+++++

SpikeEvent

Receives
++++++++

SpikeEvent, CurrentEvent, DataLoggingRequest

References
++++++++++

.. [1] Meffin H, Burkitt AN, Grayden DB (2004). An analytical
       model for the large, fluctuating synaptic conductance state typical of
       neocortical neurons in vivo. Journal of Computational Neuroscience,
       16:159-175.
       DOI: https://doi.org/10.1023/B:JCNS.0000014108.03012.81

See also
++++++++

iaf_psc_delta, iaf_psc_exp, iaf_cond_exp_sc

EndUserDocs*/

using nest::SpikeEvent;
using nest::CurrentEvent;
using nest::DataLoggingRequest;

using nest::rport;
using nest::synindex;
using nest::port;
using nest::delay;
using nest::Time;

using nest::RecordablesMap;
using nest::UniversalDataLogger;
using nest::RingBuffer;
using nest::kernel;

using nest::BadProperty;
using nest::GSLSolverFailure;
using nest::UnknownReceptorType;


class iaf_cond_exp_sc : public nest::ArchivingNode
{

public:
  iaf_cond_exp_sc();
  iaf_cond_exp_sc( const iaf_cond_exp_sc& );
  ~iaf_cond_exp_sc();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( DataLoggingRequest& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_buffers_();
  void pre_run_hook();
  void update( Time const&, const long, const long );

  struct Buffers_;

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int iaf_cond_exp_sc_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< iaf_cond_exp_sc >;
  friend class UniversalDataLogger< iaf_cond_exp_sc >;

private:
  // ----------------------------------------------------------------

  //! Model parameters
  struct Parameters_
  {
    double V_th_;    //!< Threshold Potential in mV
    double V_reset_; //!< Reset Potential in mV
    double t_ref_;   //!< Refractory period in ms
    double g_L;      //!< Leak Conductance in nS
    double C_m;      //!< Membrane Capacitance in pF
    double E_ex;     //!< Excitatory reversal Potential in mV
    double E_in;     //!< Inhibitory reversal Potential in mV
    double E_L;      //!< Leak reversal Potential (aka resting potential) in mV
    double tau_synE; //!< Time constant for excitatory synaptic kernel in ms
    double tau_synI; //!< Time constant for inhibitory synaptic kernel in ms
    double I_e;      //!< Constant Current in pA

    std::vector< Time > amp_time_stamps_;  //!< Times of amplitude changes
    std::vector< double > amp_values_;     //!< Amplitude values activated at given times
    bool allow_offgrid_amp_times_;         //!< Allow and round up amplitude times not on steps

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Buffers_& b, Node* node ); //!< Set values from dicitonary

    /**
     * Return time as Time object if valid, otherwise throw BadProperty
     *
     * @param amplitude time, ms
     * @param previous time stamp
     */
    Time validate_time_( double, const Time& );
  };

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor required because of C-style array.
   */
  struct State_
  {
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_M = 0,
      G_EXC,
      G_INH,
      STATE_VEC_SIZE
    };

    //! neuron state, must be C-array for GSL solver
    double y_[ STATE_VEC_SIZE ];
    int r_; //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );

    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  // ----------------------------------------------------------------

private:
  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( iaf_cond_exp_sc& );                  //!< Sets buffer pointers to 0
    Buffers_( const Buffers_&, iaf_cond_exp_sc& ); //!< Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< iaf_cond_exp_sc > logger_;

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;
    RingBuffer currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // Since IntergrationStep_ is initialized with step_, and the resolution
    // cannot change after nodes have been created, it is safe to place both
    // here.
    double step_;            //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    size_t I_step_idx_; //!< index of current amplitude
    double I_step_amp_; //!< current amplitude

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    int RefractoryCounts_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< iaf_cond_exp_sc > recordablesMap_;
};


inline port
iaf_cond_exp_sc::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
iaf_cond_exp_sc::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
iaf_cond_exp_sc::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
iaf_cond_exp_sc::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}


inline void
iaf_cond_exp_sc::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  ArchivingNode::get_status( d );

  ( *d )[ nest::names::recordables ] = recordablesMap_.get_list();
}

inline void
iaf_cond_exp_sc::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_;     // temporary copy in case of errors
  ptmp.set( d, B_, this );       // throws if BadProperty
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
}

} // namespace

#endif // HAVE_GSL
#endif // IAF_COND_EXP_SC_H
