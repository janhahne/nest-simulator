/*
 *  astrocyte.cpp
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


#include "astrocyte.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cmath> // in case we need isnan() // fabs
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
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

nest::RecordablesMap< nest::astrocyte > nest::astrocyte::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< astrocyte >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::IP3, &astrocyte::get_y_elem_< astrocyte::State_::IP3 > );
  insert_( names::Ca_astro, &astrocyte::get_y_elem_< astrocyte::State_::CALC > );
  insert_( names::F_IP3R, &astrocyte::get_y_elem_< astrocyte::State_::F_IP3R > );
}

extern "C" int
astrocyte_dynamics( double time, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::astrocyte::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::astrocyte& node = *( reinterpret_cast< nest::astrocyte* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double& ip3 = y[ S::IP3 ];
  const double& calc = y[ S::CALC ];
  const double& f_ip3r = y[ S::F_IP3R ];

  const double alpha_f_ip3r = node.P_.a2_ * node.P_.d2_ * ( ip3 + node.P_.d1_ ) / ( ip3 +
      node.P_.d3_ );
  const double beta_f_ip3r = node.P_.a2_ * calc;
  const double I_pump = node.P_.v3_ * std::pow(calc, 2) / (std::pow(node.P_.k3_, 2) + std::pow(calc, 2));
  const double m_inf = ip3 / (ip3 + node.P_.d1_);
  const double n_inf = calc / (calc + node.P_.d5_);
  const double calc_ER = (node.P_.c0_ - calc) / node.P_.c1_;
  const double I_leak = node.P_.c1_ * node.P_.v2_ * (calc_ER - calc);
  const double I_channel = node.P_.c1_ * node.P_.v1_ * std::pow(m_inf, 3) * std::pow(n_inf, 3) *
    std::pow(f_ip3r, 3) * (calc_ER - calc);


  // set I_gap depending on interpolation order
  double gap = 0.0;

  const double t = time / node.B_.step_;

  switch ( kernel().simulation_manager.get_wfr_interpolation_order() )
  {
  case 0:
    gap = -node.B_.sumj_g_ij_ * ip3 + node.B_.interpolation_coefficients[ node.B_.lag_ ];
    break;

  case 1:
    gap = -node.B_.sumj_g_ij_ * ip3 + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 1 ] * t;
    break;

  case 3:
    gap = -node.B_.sumj_g_ij_ * ip3 + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 1 ] * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 2 ] * t * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 3 ] * t * t * t;
    break;

  default:
    throw BadProperty( "Interpolation order must be 0, 1, or 3." );
  }

  const double I_gap = gap;

  f[ S::IP3 ] = ( node.P_.ip3_eq_ - ip3 ) / node.P_.tau_ip3_;
  f[ S::CALC ] = I_channel - I_pump + I_leak;
  f[ S::F_IP3R ] = alpha_f_ip3r * ( 1.0 - f_ip3r ) - beta_f_ip3r * f_ip3r;

  return GSL_SUCCESS;
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::astrocyte::Parameters_::Parameters_()
  : tau_ip3_( 7142.0 )   // ms
  , r_ip3_( 5.0 )      // uM / ms
  , d1_( 0.13 )    // uM
  , d2_( 1.049 )    // uM
  , d3_( 0.9434 )    // uM
  , d5_( 0.08234 )    // uM
  , v1_( 0.006 )    // 1 / ms
  , v2_( 0.00011 )    // 1 / ms
  , v3_( 0.0009 )    // uM / ms
  , k3_( 0.1 )    // uM
  , a2_( 0.0002 )    // 1 / (uM*ms)
  , c0_( 2.0 )   // uM
  , c1_( 0.185 ) 
  , ip3_eq_( 0.16 )    // uM
{
}

nest::astrocyte::State_::State_( const Parameters_& p )
{
  y_[ IP3 ] = p.ip3_eq_;
  y_[ CALC ] = 0.073;   // uM
  y_[ F_IP3R ] = 0.793;
}

nest::astrocyte::State_::State_( const State_& s )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

nest::astrocyte::State_& nest::astrocyte::State_::operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program
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
nest::astrocyte::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::tau_ip3, tau_ip3_ );
  def< double >( d, names::r_ip3, r_ip3_ );
  // TODO: add other parameters (also in nestkernel/nestnames.h and .cpp)
}

void
nest::astrocyte::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::tau_ip3, tau_ip3_ );
  updateValue< double >( d, names::r_ip3, r_ip3_ );
  // TODO: add other parameters and check range
  /*
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  */
}

void
nest::astrocyte::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::IP3, y_[ IP3 ] );
  def< double >( d, names::Ca_astro, y_[ CALC ] );
  def< double >( d, names::F_IP3R, y_[ F_IP3R ] );
}

void
nest::astrocyte::State_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::IP3, y_[ IP3 ] );
  updateValue< double >( d, names::Ca_astro, y_[ CALC ] );
  updateValue< double >( d, names::F_IP3R, y_[ F_IP3R ] );
  // TODO: check range of parameters
  /*
  if ( y_[ HH_M ] < 0 || y_[ HH_H ] < 0 || y_[ HH_N ] < 0 || y_[ HH_P ] < 0 )
  {
    throw BadProperty( "All (in)activation variables must be non-negative." );
  }
  */
}

nest::astrocyte::Buffers_::Buffers_( astrocyte& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::astrocyte::Buffers_::Buffers_( const Buffers_&, astrocyte& n )
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

nest::astrocyte::astrocyte()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

nest::astrocyte::astrocyte( const astrocyte& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

nest::astrocyte::~astrocyte()
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
nest::astrocyte::init_state_( const Node& proto )
{
  const astrocyte& pr = downcast< astrocyte >( proto );
  S_ = pr.S_;
}

void
nest::astrocyte::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.currents_.clear();  // includes resize

  // allocate strucure for gap events here
  // function is called from Scheduler::prepare_nodes() before the
  // first call to update
  // so we already know which interpolation scheme to use according
  // to the properties of this neurons
  // determine size of structure depending on interpolation scheme
  // and unsigned int Scheduler::min_delay() (number of simulation time steps
  // per min_delay step)

  // resize interpolation_coefficients depending on interpolation order
  const size_t buffer_size =
    kernel().connection_manager.get_min_delay() * ( kernel().simulation_manager.get_wfr_interpolation_order() + 1 );

  B_.interpolation_coefficients.resize( buffer_size, 0.0 );

  B_.last_y_values.resize( kernel().connection_manager.get_min_delay(), 0.0 );

  B_.sumj_g_ij_ = 0.0;

  Archiving_Node::clear_history();

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
    B_.c_ = gsl_odeiv_control_y_new( 1e-6, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-6, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = astrocyte_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::astrocyte::calibrate()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

bool
nest::astrocyte::update_( Time const& origin, const long from, const long to, const bool called_from_wfr_update )
{

  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const size_t interpolation_order = kernel().simulation_manager.get_wfr_interpolation_order();
  const double wfr_tol = kernel().simulation_manager.get_wfr_tol();
  bool wfr_tol_exceeded = false;

  // allocate memory to store the new interpolation coefficients
  // to be sent by gap event
  const size_t buffer_size = kernel().connection_manager.get_min_delay() * ( interpolation_order + 1 );
  std::vector< double > new_coefficients( buffer_size, 0.0 );
  std::vector< double > sic_values( kernel().connection_manager.get_min_delay(), 0.0 );


  // parameters needed for piecewise interpolation
  double y_i = 0.0, y_ip1 = 0.0, hf_i = 0.0, hf_ip1 = 0.0;
  double f_temp[ State_::STATE_VEC_SIZE ];

  for ( long lag = from; lag < to; ++lag )
  {

    // B_.lag is needed by astrocyte_dynamics to
    // determine the current section
    B_.lag_ = lag;

    if ( called_from_wfr_update )
    {
      y_i = S_.y_[ State_::IP3 ];
      if ( interpolation_order == 3 )
      {
        astrocyte_dynamics( 0, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_i = B_.step_ * f_temp[ State_::IP3 ];
      }
    }

    double t = 0.0;
    // TODO: deleted to much?
    const double U_old = S_.y_[ State_::IP3 ];

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

    if ( not called_from_wfr_update )
    {
      // log state data
      B_.logger_.record_data( origin.get_steps() + lag );

      // set new input current
      B_.I_stim_ = B_.currents_.get_value( lag );
    }
    else // if(called_from_wfr_update)
    {
      // check if deviation from last iteration exceeds wfr_tol
      wfr_tol_exceeded = wfr_tol_exceeded or fabs( S_.y_[ State_::IP3 ] - B_.last_y_values[ lag ] ) > wfr_tol;
      B_.last_y_values[ lag ] = S_.y_[ State_::IP3 ];

      // update different interpolations

      // constant term is the same for each interpolation order
      new_coefficients[ lag * ( interpolation_order + 1 ) + 0 ] = y_i;

      switch ( interpolation_order )
      {
      case 0:
        break;

      case 1:
        y_ip1 = S_.y_[ State_::IP3 ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = y_ip1 - y_i;
        break;

      case 3:
        y_ip1 = S_.y_[ State_::IP3 ];
        astrocyte_dynamics( B_.step_, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_ip1 = B_.step_ * f_temp[ State_::IP3 ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = hf_i;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 2 ] = -3 * y_i + 3 * y_ip1 - 2 * hf_i - hf_ip1;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 3 ] = 2 * y_i - 2 * y_ip1 + hf_i + hf_ip1;
        break;

      default:
        throw BadProperty( "Interpolation order must be 0, 1, or 3." );
      }
    }

    if ( not called_from_wfr_update )
    {
      // this is to add the incoming spikes to the state variable
      S_.y_[ State_::IP3 ] += P_.r_ip3_ * B_.spike_exc_.get_value( lag );
    }
    else
    {
      // this is to add the incoming spikes to the state variable
      S_.y_[ State_::IP3 ] += P_.r_ip3_ * B_.spike_exc_.get_value_wfr_update( lag );
    }
    double calc_thr = S_.y_[ State_::CALC ] * 1000.0 - 196.69;
    if ( calc_thr > 1.0 )
    {
      /* The SIC is converted to pA from uA/cm2 in the original publication */
      sic_values[ lag ] = std::pow(25, 2) * 3.14 * std::pow(10, -2) * std::log( calc_thr );
    }

  } // end for-loop

  // if not called_from_wfr_update perform constant extrapolation
  // and reset last_y_values
  if ( not called_from_wfr_update )
  {
    for ( long temp = from; temp < to; ++temp )
    {
      // TODO: this is for the connection between two astrocytes
      new_coefficients[ temp * ( interpolation_order + 1 ) + 0 ] = S_.y_[ State_::IP3 ];
    }

    std::vector< double >( kernel().connection_manager.get_min_delay(), 0.0 ).swap( B_.last_y_values );
  }

  // Send gap-event
  GapJunctionEvent ge;
  ge.set_coeffarray( new_coefficients );
  kernel().event_delivery_manager.send_secondary( *this, ge );
  
  // Send sic-event
  SICEvent sic;
  sic.set_coeffarray( sic_values );
  kernel().event_delivery_manager.send_secondary( *this, sic );

  // Reset variables
  B_.sumj_g_ij_ = 0.0;
  std::vector< double >( buffer_size, 0.0 ).swap( B_.interpolation_coefficients );

  return wfr_tol_exceeded;
}

void
nest::astrocyte::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
}

void
nest::astrocyte::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
}

void
nest::astrocyte::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
nest::astrocyte::handle( GapJunctionEvent& e )
{
  const double weight = e.get_weight();

  B_.sumj_g_ij_ += weight;

  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    B_.interpolation_coefficients[ i ] += weight * e.get_coeffvalue( it );
    ++i;
  }
}

#endif // HAVE_GSL
