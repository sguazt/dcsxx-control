/**
 * \file pid_controller.hpp
 *
 * \brief PID controller.
 *
 * A Proportial-Integral-Derivative (PID) controller attempts to correct the
 * error between a measured process variable and a desired setpoint by
 * calculating and then outputting a corrective action that can adjust the
 * process accordingly and rapidly, to keep the error minimal.
 *
 * The control law is defined as:
 * \f{align*}{
 * u(t) &= u_P(t)+u_I(t)+u_D(t) \\
 *      &= K_P e(t) + K_I \int_0^t{q(e(t)) e(\tau)}\,{d\tau} +
 *         K_D\frac{de}{dt}(t),
 *         \quad q(e(t)) = \begin{cases}
 *                           1,& |e(t)|\le E\\
 *                           0,& \text{otherwise}
 *                         \end{cases}
 * \f}
 * where:
 * - \f$u_P(t)\f$ is the proportional term, which determines the reaction to the
 *   current error;
 * - \f$u_I(t)\f$ is the integral term, which determines the reaction based on
 *   the sum of recent errors;
 * - \f$u_D(t)\f$ is the derivative term, which determines the reaction based on
 *   the rate at which the error has been changing.
 * - \f$K_P\f$ is the proportional gain: larger values typically mean faster
 *   response since the larger the error, the larger the proportional term
 *   compensation; however, an excessively large proportional gain will lead to
 *   system instability and oscillation.
 * - \f$K_I\f$ is the integral gain: larger values imply steady state errors are
 *   eliminated more quickly; however, the trade-off is larger overshoot, since
 *   any negative error integrated during transient response must be integrated
 *   away by positive error before we reach steady state.
 * - \f$K_D\f$ is the derivative gain: larger values decrease overshoot, but
 *   slows down transient response and may lead to instability due to signal
 *   noise amplification in the differentiation of the error.
 * - \f$q(\cdot)\f$ is a function used for turning on or off the integration
 *   when the amplitude (i.e., the absolute value) of the error signal is above
 *   a cutoff level \f$E\f$; its purpose is to mitigate the integrator
 *   windup (also called reset windup) problem, resulting whenever actuator
 *   saturation occurs, that is when integrating a large error signal during
 *   periods when the actuator is unable to respond fully to its commanded
 *   behavior; the problem is that the integrated error term can reach a very
 *   large value during these periods, which results in significant overshoot in
 *   the system response; oscillation is also a possibility, with the actuator
 *   banging from one end of its range of motion to the other.
 * .
 *
 * \note
 *  Actually, this implementation is based on the discretized version of the
 *  above control law, that is:
 *  \f[
 *    u(t+\Delta\,t) = K_P e(t+\Delta\,t) + K_I \Delta\,t \sum_{i=0}^{t+\Delta\,t}{e(t+i\Delta\,t)} + K_D\frac{e(t+\Delta\,t)-e(t)}{\Delta\,t}
 *  \f]
 *  where \f$\Delta\,t\f$ is the sampling time</em> (also known as <em>control
 *  time</em> or <em>step time</em>).
 *
 * Copyright (C) 2009-2010  Distributed Computing System (DCS) Group, Computer
 * Science Department - University of Piemonte Orientale, Alessandria (Italy).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#ifndef DCS_CONTROL_PID_CONTROLLER_HPP
#define DCS_CONTROL_PID_CONTROLLER_HPP


#include <cmath>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <functional>
#include <limits>
#include <stdexcept>


namespace dcs { namespace control {


/**
 * \brief Class for single-loop PID controllers.
 *
 * A single-loop PID controller can be used to control Single-Input Single-Output
 * (SISO) systems.
 *
 * The control law is thus given by:
 * \f[
 *   u(k) = K_Pe(k)+K_I\sum_{i=0}^k{e(k)}+K_D(e(k)-e(k-1))
 * \f]
 * where:
 * - \f$u_j(k)\f$ is the <em>control signal</em> at time \f$k\f$,
 * - \f$e_j(k)\f$ is the <em>control error</em> at time \f$k\f$,
 * - \f$K_{P,j}\f$ is the <em>proportional gain</em>,
 * - \f$K_{I,j}\f$ is the <em>integral gain</em>,
 * - \f$K_{D,j}\f$ is the <em>derivative gain</em>.
 * .
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */
template <typename RealT = double>
class pid_controller
{
	public: typedef RealT real_type;

	/**
	 * \brief A constructor.
	 *
	 * \param kp The proportional gain.
	 * \param ki The integral gain.
	 * \param kd The derivative gain.
	 * \param step_time The control time.
	 * \param err_thresh The error cutoff level used for mitigating the
	 *  integrator windup problem; set to \f$\infty\f$ for not using it.
	 */
	public: pid_controller(real_type kp, real_type ki, real_type kd, real_type step_time, real_type err_thresh = ::std::numeric_limits<real_type>::infinity())
		: kp_(kp),
		  ki_(ki),
		  kd_(kd),
		  ts_(step_time),
		  err_thresh_(err_thresh),
		  prev_err_(0),
		  integral_(0),
		  started_(false)
	{
		// Empty
	}


	public: real_type sampling_time() const
	{
		return ts_;
	}


	public: real_type proportional_gain() const
	{
		return kp_;
	}


	public: real_type integral_gain() const
	{
		return ki_;
	}


	public: real_type derivative_gain() const
	{
		return kd_;
	}


	/**
	 * \brief Update the control signal according to the given control error.
	 * \param error The current control error.
	 * \return The new control signal.
	 *
	 * \note
	 *  This version uses the positional algorithm:
	 *  \f{align*}{
	 *    u(k) &= u_P(k) + u_I(k) + u_D(k) \\
	 *         &= K_P e(k) + K_I \sum_{i=0}^{k} e(i) + K_D (e(k)-e(k-1))
	 *  \f}
	 */
	public: real_type control(real_type error)
	{
		// Update if the error magnitude is below the threshold.
		if (
				err_thresh_ == ::std::numeric_limits<real_type>::infinity()
				|| ::std::abs(error) < err_thresh_
		) {
			// Update the error integral according to the forward difference
			// equation:
			//   (u_I(k+h) - u_I(k))h = K_I ==> u_I(k+h) = u_I(k) + h*e(k+h)*K_I
			integral_ += ts_*error;
		}

		// Compute the error derivative according to the backward difference
		// equation:
		// K_D*(u_D(k)-u_D(k-1))/h+u_D(k) = -K_D*e(k)/h
		real_type deriv(0);
		if (started_)
		{
			deriv = (error - prev_err_) / ts_;
		}

		prev_err_ = error;

		if (!started_)
		{
			started_ = true;
		}

		// Return the PID controller actuator command
		return kp_*error + ki_*integral_ + kd_*deriv;
	}


	// Use the incremental algorithm (aka velocity algorithm)
	// u(t_k)=u(t_{k-1})+K_p\left[\left(1+\dfrac{\Delta t}{T_i}+\dfrac{T_d}{\Delta t}\right)e(t_k)+\left(-1-\dfrac{2T_d}{\Delta t}\right)e(t_{k-1})+\dfrac{T_d}{\Delta t}e(t_{k-2})\right] 
	//real_type control(real_type error);


	/// Proportional gain.
	private: real_type kp_;
	/// Integral gain.
	private: real_type ki_;
	/// Differential gain.
	private: real_type kd_;
	/// Controller step time.
	private: real_type ts_;
	/// Error cutoff level for preventing the integrator windup problem.
	private: real_type err_thresh_;
	/// The control error at the previous control interval.
	private: real_type prev_err_;
	/// Integrator state.
	private: real_type integral_;
	/// Tell if the controller is started.
	private: bool started_;
};


/**
 * \brief Class for multi-loop PID controllers.
 *
 * A multi-loop PID controller can be used to control Multi-Input Multi-Output
 * (MIMO) systems.
 * In a multi-loop PID controller, control is performed by considering each
 * single-loop PID controller independently.
 *
 * The control law is thus given by:
 * \f[
 *   u_j(k) = K_{P,j}e_j(k)+K_{I,j}\sum_{i=0}^k{e_j(k)}+K_{D,j}(e_j(k)-e_j(k-1))
 * \f]
 * for \f$j=1,\ldots,n_u\f$, where:
 * - \f$u_j(k)\f$ is the <em>control signal</em> at time \f$k\f$ for the
 *   \f$j\f$-th single-loop PID controller,
 * - \f$e_j(k)\f$ is the <em>control error</em> at time \f$k\f$ for the
 *   \f$j\f$-th single-loop PID controller,
 * - \f$K_{P,j}\f$ is the <em>proportional gain</em> for the \f$j\f$-th
 *   single-loop PID controller,
 * - \f$K_{I,j}\f$ is the <em>integral gain</em> for the \f$j\f$-th single-loop
 *   PID controller,
 * - \f$K_{D,j}\f$ is the <em>derivative gain</em> for the \f$j\f$-th
 *   single-loop PID controller.
 * .
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */
template <typename VectorT, typename RealT = double>
class multiloop_pid_controller
{
	public: typedef VectorT vector_type;
	public: typedef RealT real_type;

	/**
	 * \brief A constructor.
	 *
	 * \param kp The proportional gain matrix.
	 * \param ki The integral gain matrix.
	 * \param kd The derivative gain matrix.
	 * \param step_time The control (step) time.
	 * \param err_thresh The error cutoff level used for mitigating the
	 *  integrator windup problem; set to \f$\infty\f$ for not using it.
	 */
	public: multiloop_pid_controller(vector_type kp, vector_type ki, vector_type kd, real_type step_time, vector_type err_thresh = vector_type(1, ::std::numeric_limits<real_type>::infinity()))
		: kp_(kp),
		  ki_(ki),
		  kd_(kd),
		  ts_(step_time),
		  err_thresh_(err_thresh),
		  prev_err_(::boost::numeric::ublasx::size(kp), 0),
		  integral_(::boost::numeric::ublasx::size(ki), 0),
		  started_(false)
	{
		// preconditions
		DCS_ASSERT(
			::boost::numeric::ublasx::size(kp) == ::boost::numeric::ublasx::size(ki)
			&&
			::boost::numeric::ublasx::size(kp) == ::boost::numeric::ublasx::size(kd),
			throw ::std::invalid_argument("Non conformant dimensions for K_P, K_I and K_D gains.")
		);

		// Make sure the error threshold vector has a compliant size.

		typedef typename ::boost::numeric::ublas::vector_traits<vector_type>::size_type size_type;

		size_type n = ::boost::numeric::ublasx::size(kp);
		size_type m = ::boost::numeric::ublasx::size(err_thresh_);
		if (m < n)
		{
			err_thresh_.resize(n, true);
			real_type val = (m > 0) ? err_thresh_(m-1) : ::std::numeric_limits<real_type>::infinity();
			while (m < n)
			{
				err_thresh_[m++] = val;
			}
		}
	}


	public: real_type sampling_time() const
	{
		return ts_;
	}


	public: vector_type proportional_gain() const
	{
		return kp_;
	}


	public: vector_type integral_gain() const
	{
		return ki_;
	}


	public: vector_type derivative_gain() const
	{
		return kd_;
	}


	/**
	 * \brief Compute the control signal according to the given control error.
	 * \param error The current control error.
	 * \return The new control signal.
	 *
	 * \note
	 *  This version uses the positional algorithm:
	 *  \f{align*}{
	 *    u(k) &= u_P(k) + u_I(k) + u_D(k) \\
	 *         &= K_P e(k) + K_I \sum_{i=0}^{k} e(i) + K_D (e(k)-e(k-1))
	 *  \f}
	 */
	public: vector_type control(vector_type error)
	{
		// preconditions
		DCS_ASSERT(
			::boost::numeric::ublasx::size(error) == ::boost::numeric::ublasx::size(err_thresh_),
			throw ::std::invalid_argument("Wrong length of error vector.")
		);

		typedef typename vector_type::size_type size_type;

		// Update if the error magnitude is below the threshold.
		size_type err_len = ::boost::numeric::ublasx::size(error);
		for (size_type i = 0; i < err_len; ++i)
		{
			if (
				(err_thresh_(i) == ::std::numeric_limits<real_type>::infinity())
				||
				(::std::abs(error(i)) < err_thresh_(i))
			) {
				// Update the error integral according to the forward difference
				// equation:
				//   (u_I(k+h) - u_I(k))/h = K_I ==> u_I(k+h) = u_I(k) + h*e(k+h)*K_I
				integral_(i) += ts_*error(i);
			}
		}

		// Compute the error derivative according to the backward difference
		// equation:
		// K_D*(u_D(k)-u_D(k-h))/h+u_D(k) = -K_D*e(k)/h
		vector_type deriv(err_len, 0);
		if (started_)
		{
			deriv = (error - prev_err_) / ts_;
		}

		prev_err_ = error;

		if (!started_)
		{
			started_ = true;
		}

		// Return the PID controller actuator command
		return ::boost::numeric::ublas::element_prod(kp_, error)
			   + ::boost::numeric::ublas::element_prod(ki_, integral_)
			   + ::boost::numeric::ublas::element_prod(kd_, deriv);
	}


	// Use the incremental algorithm (aka velocity algorithm)
	// u(k) = u(k-1) + K_P \left[\left(1+\frac{\Delta t}{T_I}+\frac{T_D}{\Delta t}\right)e(k)+\left(-1-\frac{2T_D}{\Delta t}e(k-1)+\frac{T_D}{\Delta t}e(k-2)\right)\right]
	//vector_type control(vector_type error);


	/// Proportional gain.
	private: vector_type kp_;
	/// Integral gain.
	private: vector_type ki_;
	/// Differential gain.
	private: vector_type kd_;
	/// Control (step) time.
	private: real_type ts_;
	/// Error cutoff level for preventing the integrator windup problem.
	private: vector_type err_thresh_;
	/// The control error at the previous control interval.
	private: vector_type prev_err_;
	/// Integrator state.
	private: vector_type integral_;
	/// Tell if the controller is started.
	private: bool started_;
};


/**
 * \brief Class for Multi-Input Multi-Output (MIMO) PID controllers.
 *
 * A MIMO PID controller can be used to control MIMO systems.
 * A MIMO PID controller differs from a multi-loop PID controller since it takes
 * care of possible interactions between control signals. 
 * The control law is thus given by:
 * \f[
 *   \mathbf{u}(k) = \mathbf{K}_P \mathbf{e}(k) + \mathbf{K}_I \sum_{i=0}^k{\mathbf{e}(k)} + \mathbf{K}_D(\mathbf{e}(k)-\mathbf{e}(k-1))
 * \f]
 * where:
 * - \f$\mathbf{u}(k)\f$ is the <em>control signal vector</em> at time \f$k\f$,
 * - \f$\mathbf{e}(k)\f$ is the <em>control error vector</em> at time \f$k\f$,
 * - \f$\mathbf{K}_P\f$ is the <em>proportional gain matrix</em>,
 * - \f$\mathbf{K}_I\f$ is the <em>integral gain matrix</em>,
 * - \f$\mathbf{K}_D\f$ is the <em>derivative gain matrix</em>.
 * .
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */
template <typename VectorT, typename MatrixT, typename RealT = double>
class mimo_pid_controller
{
	public: typedef VectorT vector_type;
	public: typedef MatrixT matrix_type;
	public: typedef RealT real_type;

	/**
	 * \brief A constructor.
	 *
	 * \param kp The proportional gain matrix.
	 * \param ki The integral gain matrix.
	 * \param kd The derivative gain matrix.
	 * \param step_time The control (step) time.
	 * \param err_thresh The error cutoff level used for mitigating the
	 *  integrator windup problem; set to \f$\infty\f$ for not using it.
	 */
	public: mimo_pid_controller(matrix_type kp, matrix_type ki, matrix_type kd, real_type step_time, vector_type err_thresh = vector_type(::std::numeric_limits<real_type>::infinity()))
		: kp_(kp),
		  ki_(ki),
		  kd_(kd),
		  ts_(step_time),
		  err_thresh_(err_thresh),
		  prev_err_(::boost::numeric::ublasx::num_rows(kp), 0),
		  integral_(::boost::numeric::ublasx::num_rows(ki), 0),
		  started_(false)
	{
		// preconditions
		DCS_ASSERT(
			::boost::numeric::ublasx::num_rows(kp) == ::boost::numeric::ublasx::num_rows(ki)
			&&
			::boost::numeric::ublasx::num_columns(kp) == ::boost::numeric::ublasx::num_columns(ki)
			&&
			::boost::numeric::ublasx::num_rows(kp) == ::boost::numeric::ublasx::num_rows(kd)
			&&
			::boost::numeric::ublasx::num_columns(kp) == ::boost::numeric::ublasx::num_columns(kd),
			throw ::std::invalid_argument("Non conformant dimensions for K_P, K_I and K_D gains.")
		);

		// Make sure the error threshold vector has a compliant size.

		typedef typename ::boost::numeric::ublas::vector_traits<vector_type>::size_type size_type;

		size_type n = ::boost::numeric::ublasx::num_rows(kp);
		size_type m = ::boost::numeric::ublasx::size(err_thresh_);
		if (m < n)
		{
			err_thresh_.resize(n, true);
			real_type val = (m > 0) ? err_thresh_(m-1) : ::std::numeric_limits<real_type>::infinity();
			while (m < n)
			{
				err_thresh_[m++] = val;
			}
		}
	}


	public: real_type sampling_time() const
	{
		return ts_;
	}


	public: matrix_type proportional_gain() const
	{
		return kp_;
	}


	public: matrix_type integral_gain() const
	{
		return ki_;
	}


	public: matrix_type derivative_gain() const
	{
		return kd_;
	}


	/**
	 * \brief Update the control signal according to the given control error.
	 * \param error The current control error.
	 * \return The new control signal.
	 *
	 * \note
	 *  This version uses the positional algorithm:
	 *  \f{align*}{
	 *    u(k) &= u_P(k) + u_I(k) + u_D(k) \\
	 *         &= K_P e(k) + K_I \sum_{i=0}^{k} e(i) + K_D (e(k)-e(k-1))
	 *  \f}
	 */
	public: vector_type control(vector_type error)
	{
		// preconditions
		DCS_ASSERT(
			::boost::numeric::ublasx::size(error) == ::boost::numeric::ublasx::size(err_thresh_),
			throw ::std::invalid_argument("Wrong length of error vector.")
		);

		typedef typename vector_type::size_type size_type;

		// Update if the error magnitude is below the threshold.
		size_type err_len = ::boost::numeric::ublasx::size(error);
		for (size_type i = 0; i < err_len; ++i)
		{
			if (
				(err_thresh_(i) == ::std::numeric_limits<real_type>::infinity())
				||
				(::std::abs(error(i)) < err_thresh_(i))
			) {
				// Update the error integral according to the forward difference
				// equation:
				//   (u_I(k+h) - u_I(k))/h = K_I ==> u_I(k+h) = u_I(k) + h*e(k+h)*K_I
				integral_(i) += ts_*error(i);
			}
		}

		// Compute the error derivative according to the backward difference
		// equation:
		// K_D*(u_D(k)-u_D(k-h))/h+u_D(k) = -K_D*e(k)/h
		vector_type deriv(err_len, 0);
		if (started_)
		{
			deriv = (error - prev_err_) / ts_;
		}

		prev_err_ = error;

		if (!started_)
		{
			started_ = true;
		}

		// Return the PID controller actuator command
		return ::boost::numeric::ublas::prod(kp_, error)
			   + ::boost::numeric::ublas::prod(ki_, integral_)
			   + ::boost::numeric::ublas::prod(kd_, deriv);
	}


	// Use the incremental algorithm (aka velocity algorithm)
	// u(k) = u(k-1) + K_P \left[\left(1+\frac{\Delta t}{T_I}+\frac{T_D}{\Delta t}\right)e(k)+\left(-1-\frac{2T_D}{\Delta t}e(k-1)+\frac{T_D}{\Delta t}e(k-2)\right)\right]
	//vector_type control(vector_type error);


	/// Proportional gain.
	private: matrix_type kp_;
	/// Integral gain.
	private: matrix_type ki_;
	/// Differential gain.
	private: matrix_type kd_;
	/// Control (step) time.
	private: real_type ts_;
	/// Error cutoff level for preventing the integrator windup problem.
	private: vector_type err_thresh_;
	/// The control error at the previous control interval.
	private: vector_type prev_err_;
	/// Integrator state.
	private: vector_type integral_;
	/// Tell if the controller is started.
	private: bool started_;
};


/**
 * \brief PID controller in standard form.
 *
 * A single-loop PID controller can be used to control Single-Input Single-Output
 * (SISO) systems.
 *
 * The control law is thus given by:
 * \f[
 *   u(k) = K_P\left(e(k)+\frac{1}{T_I}\sum_{i=0}^k{e(i)}+T_D(e(k)-e(k-1)\right))
 * \f]
 * where:
 * - \f$u_j(k)\f$ is the <em>control signal</em> at time \f$k\f$,
 * - \f$e_j(k)\f$ is the <em>control error</em> at time \f$k\f$,
 * - \f$K_{P,j}\f$ is the <em>proportional gain</em>,
 * - \f$T_I\f$ is the <em>integral time</em>,
 * - \f$T_D\f$ is the <em>derivative time</em>.
 * .
 *
 * The ideal and standard PID controller are related by:
 * \f{align*}{
 *   K_P &= K_P \\
 *   K_I &= \frac{K_P}{T_I} \\
 *   K_D &= K_P T_D
 * \f}
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */
template <typename RealT = double>
class std_pid_controller
{
	public: typedef RealT real_type;


	public: std_pid_controller(real_type kp, real_type int_time, real_type deriv_time, real_type step_time, real_type err_thresh = ::std::numeric_limits<real_type>::infinity())
		: pid_(kp, kp/int_time, kp*deriv_time, err_thresh)
	{
		// Empty
	}


	public: real_type sampling_time() const
	{
		return pid_.sampling_time();
	}


	public: real_type proportional_gain() const
	{
		return pid_.proportional_gain();
	}


	public: real_type integral_time() const
	{
		return pid_.proportional_gain()/pid_.integral_gain();
	}


	public: real_type derivative_time() const
	{
		return pid_.derivative_gain()/pid_.proportional_gain();
	}


	public: real_type control(real_type error)
	{
		return pid_.control(error);
	}


	private: pid_controller<real_type> pid_;
};

}} // dcs::control


#endif // DCS_CONTROL_PID_CONTROLLER_HPP
