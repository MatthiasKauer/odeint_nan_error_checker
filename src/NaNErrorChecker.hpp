/*
 * NaNErrorChecker.cpp
 *
 *  Created on: Apr 13, 2013
 *      Author: matthias
 */
#ifndef NANERRORCHECKER_HPP_
#define NANERRORCHECKER_HPP_

#include <algorithm>
#include <cmath>
#include <cmath>

#include <boost/config.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/copy.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>

#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

//called in controlled_runge_kutta.hpp line 710 as:
//value_type max_rel_err = m_error_checker.error( m_stepper.algebra() , in , dxdt_in , m_xerr.m_v , dt );

//MATTHIAS' ADDITIONS

//called in controlled_runge_kutta.hpp line 710 as:
//value_type max_rel_err = m_error_checker.error( m_stepper.algebra() , in , dxdt_in , m_xerr.m_v , dt );
namespace boost {
namespace numeric {
namespace odeint {

template
<
class Value ,
class Algebra = range_algebra ,
class Operations = default_operations
>
class nan_error_checker
{
public:

    typedef Value value_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;

    nan_error_checker(
            value_type eps_abs = static_cast< value_type >( 1.0e-6 ) ,
            value_type eps_rel = static_cast< value_type >( 1.0e-6 ) ,
            value_type a_x = static_cast< value_type >( 1 ) ,
            value_type a_dxdt = static_cast< value_type >( 1 ) )
    : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
    { }


    template< class State , class Deriv , class Err , class Time >
    value_type error( const State &x_old , const Deriv &dxdt_old , Err &x_err , Time dt ) const
    {
        return error( algebra_type() , x_old , dxdt_old , x_err , dt );
    }

    template<class State, class Deriv, class Err, class Time>
    value_type error(algebra_type &algebra, const State &x_old,
    		const Deriv &dxdt_old, Err &x_err, Time dt) const {
    	// this overwrites x_err !
    	algebra.for_each3(x_err, x_old, dxdt_old,
    			typename operations_type::template rel_error<value_type>(m_eps_abs,
    					m_eps_rel, m_a_x, m_a_dxdt * get_unit_value(dt)));

    	bool nonan=std::none_of(boost::begin(x_err), boost::begin(x_err), [](value_type v) {return std::isnan(v);} );

    	if(nonan) {
    		value_type res = algebra.reduce( x_err ,
    				typename operations_type::template maximum< value_type >() , static_cast< value_type >( 0 ) );
    		return res;
    	}
    	else {
    		return 2.0; //if( max_rel_err > 1.0 ) is checked in l.712 controlled_runge_kutta.hpp
    	}
    }
private:

    value_type m_eps_abs;
    value_type m_eps_rel;
    value_type m_a_x;
    value_type m_a_dxdt;

};

}}}
#endif
