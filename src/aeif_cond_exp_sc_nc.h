/*
 *  aeif_cond_exp_sc_nc.h
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

#ifndef AEIF_COND_EXP_SC_NC_H
#define AEIF_COND_EXP_SC_NC_H

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
#include "random_generators.h" // NOISE

namespace stepcurrentmodule
{
/**
 * Function computing right-hand side of ODE for GSL solver if Delta_T != 0.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int aeif_cond_exp_sc_nc_dynamics( double, const double*, double*, void* );

/**
 * Function computing right-hand side of ODE for GSL solver if Delta_T == 0.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int aeif_cond_exp_sc_nc_dynamics_DT0( double, const double*, double*, void* );

/* BeginUserDocs: neuron, adaptive threshold, integrate-and-fire, conductance-based

Short description
+++++++++++++++++

Conductance based exponential integrate-and-fire neuron model

Description
+++++++++++

``aeif_cond_exp_sc_nc`` is the adaptive exponential integrate and fire neuron
according to Brette and Gerstner (2005), with postsynaptic
conductances in the form of truncated exponentials;  it includes a built-in step current generator.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg
solver with adaptive stepsize to integrate the differential equation.

The membrane potential is given by the following differential equation:

.. math::

 C dV/dt= -g_L(V-E_L)+g_L \cdot \Delta_T \cdot \exp((V-V_T)/\Delta_T)-g_e(t)(V-E_e) \\
                                                     -g_i(t)(V-E_i)-w + I_e

and

.. math::

 \tau_w \cdot dw/dt= a(V-E_L) -W

Note that the spike detection threshold V_peak is automatically set to
:math:`V_th+10 mV` to avoid numerical instabilites that may result from
setting V_peak too high.

For implementation details see the
`aeif_models_implementation <../model_details/aeif_models_implementation.ipynb>`_ notebook.

See also [1]_.


Parameters:
+++++++++++++
The following parameters can be set in the status dictionary.

======== ======= =======================================
**Dynamic state variables:**
--------------------------------------------------------
 V_m     mV      Membrane potential
 g_ex    nS      Excitatory synaptic conductance
 g_in    nS      Inhibitory synaptic conductance
 w       pA      Spike-adaptation current
======== ======= =======================================

======== ======= =======================================
**Membrane Parameters**
--------------------------------------------------------
 C_m     pF      Capacity of the membrane
 t_ref   ms      Duration of refractory period
 V_reset mV      Reset value for V_m after a spike
 E_L     mV      Leak reversal potential
 g_L     nS      Leak conductance
 I_e     pA      Constant external input current
======== ======= =======================================

======== ======= ==================================
**Spike adaptation parameters**
---------------------------------------------------
 a       nS      Subthreshold adaptation
 b       pA      Spike-triggered adaptation
 Delta_T mV      Slope factor
 tau_w   ms      Adaptation time constant
 V_th    mV      Spike initiation threshold
 V_peak  mV      Spike detection threshold
======== ======= ==================================

=========== ======= ===========================================================
**Synaptic parameters**
-------------------------------------------------------------------------------
 E_ex       mV      Excitatory reversal potential
 tau_syn_ex ms      Exponential decay time constant of excitatory synaptic
                    conductance kernel
 E_in       mV      Inhibitory reversal potential
 tau_syn_in ms      Exponential decay time constant of inhibitory synaptic
                    conductance kernel
=========== ======= ===========================================================

============= ======= =========================================================
**Integration parameters**
-------------------------------------------------------------------------------
gsl_error_tol real    This parameter controls the admissible error of the
                      GSL integrator. Reduce it if NEST complains about
                      numerical instabilities.
============= ======= =========================================================

==================== ======  =======================================================
**Step current parameters**
-------------------------------------------------------------------------------
 amplitude_times     ms      Times at which current changes (list)
 amplitude_values    pA      Amplitudes of step current current (list)
 allow_offgrid_times         If True, allow off-grid times (default: False)
==================== ======  =======================================================

==================== ======  =======================================================
**Noise current parameters**
-------------------------------------------------------------------------------
 noise_mean          pA      The mean value :math:`\mu` of the noise current NOISE
 noise_std           pA      The standard deviation :math:`\sigma` of the noise current NOISE
 noise_dt            ms      The interval :math:`\delta` between changes in current (default: 10 * resolution) NOISE
 noise_std_mod       pA      The modulation :math:`\sigma_{\text{mod}}` of the standard deviation of the noise current NOISE
 noise_frequency     Hz      The frequency of the sine modulation NOISE
 noise_phase         deg     The phase of sine modulation (0–360) NOISE
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

.. [1] Brette R and Gerstner W (2005). Adaptive Exponential
       Integrate-and-Fire Model as an Effective Description of Neuronal
       Activity. J Neurophysiol 94:3637-3642.
       DOI: https://doi.org/10.1152/jn.00686.2005


See also
++++++++

iaf_cond_exp_sc_nc, aeif_cond_alpha

EndUserDocs */

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
using nest::NumericalInstability;
using nest::UnknownReceptorType;

using nest::normal_distribution; // NOISE
using nest::get_vp_specific_rng; // NOISE
using nest::get_vp_specific_rng; // NOISE


class aeif_cond_exp_sc_nc : public nest::ArchivingNode
{

public:
  aeif_cond_exp_sc_nc();
  aeif_cond_exp_sc_nc( const aeif_cond_exp_sc_nc& );
  ~aeif_cond_exp_sc_nc();

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
  void update( const Time&, const long, const long );

  struct Buffers_;

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int aeif_cond_exp_sc_nc_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< aeif_cond_exp_sc_nc >;
  friend class UniversalDataLogger< aeif_cond_exp_sc_nc >;

private:
  // ----------------------------------------------------------------

  //! Independent parameters
  struct Parameters_
  {
    double V_peak_;  //!< Spike detection threshold in mV
    double V_reset_; //!< Reset Potential in mV
    double t_ref_;   //!< Refractory period in ms

    double g_L;        //!< Leak Conductance in nS
    double C_m;        //!< Membrane Capacitance in pF
    double E_ex;       //!< Excitatory reversal Potential in mV
    double E_in;       //!< Inhibitory reversal Potential in mV
    double E_L;        //!< Leak reversal Potential (aka resting potential) in mV
    double Delta_T;    //!< Slope factor in mV
    double tau_w;      //!< Adaptation time-constant in ms
    double a;          //!< Subthreshold adaptation in nS
    double b;          //!< Spike-triggered adaptation in pA
    double V_th;       //!< Spike threshold in mV
    double tau_syn_ex; //!< Excitatory synaptic kernel decay time in ms
    double tau_syn_in; //!< Inhibitory synaptic kernel decay time in ms
    double I_e;        //!< Intrinsic current in pA

    double gsl_error_tol; //!< Error bound for GSL integrator

    std::vector< Time > amp_time_stamps_;  //!< Times of amplitude changes
    std::vector< double > amp_values_;     //!< Amplitude values activated at given times
    bool allow_offgrid_amp_times_;         //!< Allow and round up amplitude times not on steps

    double noise_mean_;    //!< mean current, in pA NOISE
    double noise_std_;     //!< standard deviation of current, in pA NOISE
    double noise_std_mod_; //!< standard deviation of current modulation, in pA NOISE
    double noise_freq_;    //!< Standard frequency in Hz NOISE
    double noise_phi_deg_; //!< Phase of sinusodial noise modulation (0-360 deg) NOISE
    Time noise_dt_;        //!< time interval between updates NOISE


    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Buffers_& b, Node* node ); //!< Set values from dictionary

    /**
     * Return time as Time object if valid, otherwise throw BadProperty
     *
     * @param amplitude time, ms
     * @param previous time stamp
     */
    Time validate_time_( double, const Time& );
    Time get_default_dt(); // NOISE
  };

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor required because of C-style array.
   */
  struct State_
  {
    /**
     * Enumeration identifying elements in state array State_::y_.
     * The state vector must be passed to GSL as a C array. This enum
     * identifies the elements of the vector. It must be public to be
     * accessible from the iteration function.
     */
    enum StateVecElems
    {
      V_M = 0,
      G_EXC, // 1
      G_INH, // 2
      W,     // 3
      STATE_VEC_SIZE
    };

    //! neuron state, must be C-array for GSL solver
    double y_[ STATE_VEC_SIZE ];
    unsigned int r_; //!< number of refractory steps remaining

    double y_0_; // NOISE
    double y_1_; // NOISE

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
    Buffers_( aeif_cond_exp_sc_nc& );                  //!< Sets buffer pointers to 0
    Buffers_( const Buffers_&, aeif_cond_exp_sc_nc& ); //!< Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< aeif_cond_exp_sc_nc > logger_;

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;
    RingBuffer currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing the GSL system

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

    long next_step_; //!< time step of next change in current NOISE
    double I_noise_amp_;   // noise amplitude NOISE
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    /**
     * Threshold detection for spike events: P.V_peak if Delta_T > 0.,
     * P.V_th if Delta_T == 0.
     */
    double V_peak;

    unsigned int refractory_counts_;

    normal_distribution normal_dist_; //!< normal distribution NOISE

    long dt_steps_;  //!< update interval in steps NOISE
    double omega_;   //!< frequency [radian/s] NOISE
    double phi_rad_; //!< phase of sine current (0-2Pi rad) NOISE

    // The exact integration matrix NOISE
    double A_00_; 
    double A_01_;
    double A_10_;
    double A_11_;

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
  static RecordablesMap< aeif_cond_exp_sc_nc > recordablesMap_;
};

inline port
aeif_cond_exp_sc_nc::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
aeif_cond_exp_sc_nc::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
aeif_cond_exp_sc_nc::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
aeif_cond_exp_sc_nc::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
aeif_cond_exp_sc_nc::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  ArchivingNode::get_status( d );

  ( *d )[ nest::names::recordables ] = recordablesMap_.get_list();
}

inline void
aeif_cond_exp_sc_nc::set_status( const DictionaryDatum& d )
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

inline Time
aeif_cond_exp_sc_nc::Parameters_::get_default_dt()
{
  return 10 * Time::get_resolution();
}

} // namespace

#endif // HAVE_GSL
#endif // AEIF_COND_EXP_SC_H
