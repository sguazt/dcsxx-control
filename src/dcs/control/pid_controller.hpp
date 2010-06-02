/**
 * \file pid_controller.hpp
 *
 * \brief PID controller.
 *
 * \_licinclude docs/LICENSE-embed.txt
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#include <cmath>
#include <limits>

namespace dcs { namespace control {

/**
 * \brief Class for PID controllers.
 *
 * A Proportial-Integral-Derivative (PID) controller attempts to correct the
 * error between a measured process variable and a desired setpoint by
 * calculating and then outputting a corrective action that can adjust the
 * process accordingly and rapidly, to keep the error minimal.
 *
 * The control law is defined as:
 * \f{align*}{
 * u(t) &= u_P(t)+u_I(t)+u_D(t) \\
 *      &= K_P e(t) + K_I \int_0^t{q(e(t)) e(\tau)}\,{d\tau} + K_D\frac{de}{dt}(t), \quad q(e(t)) = \begin{cases}1,& |e(t)|\le E \\ 0,& \text{otherwise}\end{cases}
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
	public: pid_controller(real_type kp, real_type ki, real_type kd, real_type step_time, real_type err_thresh = std::numeric_limits<real_type>::infinity());

	/**
	 * \brief Update the control input according to the given control error.
	 * \param error The current control error.
	 * \return The new control input.
	 */
	public: real_type update(real_type error);

	/// Proportional gain.
	private: real_type kp_;
	/// Integral gain.
	private: real_type ki_;
	/// Differential gain.
	private: real_type kd_;
	/// Controller step time.
	private: real_type h_;
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
 * \brief PID controller in standard form.
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */
template <typename RealT = double>
class std_pid_controller
{
	public: typedef RealT real_type;

	public: std_pid_controller(real_type kp, real_type int_time, real_type deriv_time, real_type step_time, real_type err_thresh = std::numeric_limits<real_type>::infinity());
	public: real_type update(real_type error);

	private: pid_controller<real_type> pid_;
};


//@{ pid_controller ////////////////////////////////////////////////////////////

template <typename RealT>
pid_controller<RealT>::pid_controller(RealT kp, RealT ki, RealT kd, RealT step_time, RealT err_thresh)
	: kp_(kp),
	  ki_(ki),
	  kd_(kd),
	  h_(step_time),
	  err_thresh_(err_thresh),
	  prev_err_(0),
	  integral_(0),
	  started_(false)
{
	// Empty
}

// Use the positional algorithm:
// u(k) = u_P(k) + u_I(k) + u_D(k)
//      = K_P e(k) + K_I sum_{i=1}^{k} e(i) + K_D (e(k)-e(k-1))
template <typename RealT>
RealT pid_controller<RealT>::update(RealT error)
{
	// Update if the error magnitude is below the threshold.
	if (
			err_thresh_ == std::numeric_limits<RealT>::infinity()
			|| std::fabs(error) < err_thresh_
	) {
		// Update the error integral according to the forward difference
		// equation:
		//   (u_I(k+1) - u_I(k))/h = K_I ==> u_I(k+1) = u_I(k) + h*e(k+1)*K_I
		//integral_ += h_*error;
		integral_ += error;
	}

	// Compute the error derivative according to the backward difference
	// equation:
	// K_D*(u_D(k)-u_D(k-1))/h+u_D(k) = -K_D*e(k)/h
	RealT deriv;
	if (!started_)
	{
		started_ = true;
		deriv = 0;
	}
	else
	{
		//deriv = (error - prev_err_) / h_;
		deriv = error - prev_err_;
	}

	prev_err_ = error;

	// Return the PID controller actuator command
	return kp_*error + ki_*integral_ + kd_*deriv;
}

// Use the incremental algorithm (aka velocity algorithm)
// u(k) = u(k-1) + K_P \left[\left(1+\frac{\Delta t}{T_I}+\frac{T_D}{\Delta t}\right)e(k)+\left(-1-\frac{2T_D}{\Delta t}e(k-1)+\frac{T_D}{\Delta t}e(k-2)\right)\right]
//template <typename RealT>
//RealT pid_controller<RealT>update(RealT error)
//{
//
//}

//@} pid_controller ////////////////////////////////////////////////////////////


//@{ std_pid_controller ////////////////////////////////////////////////////////

template <typename RealT>
std_pid_controller<RealT>::std_pid_controller(RealT kp, RealT int_time, RealT deriv_time, RealT step_time, RealT err_thresh)
	: pid_(kp, kp/int_time, kp*deriv_time, err_thresh)
{
	// Empty
}

template <typename RealT>
RealT std_pid_controller<RealT>::update(RealT error)
{
	return pid_.update(error);
}

//@} std_pid_controller ////////////////////////////////////////////////////////

}} // dcs::control
