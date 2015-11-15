/**
 * \file dcs/control/design/matlab_linear_mpc.hpp
 *
 * \brief Linear Model Predictive Control based on the MATLAB Model Predictive
 *  Control toolbox.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2015 Marco Guazzone (marco.guazzone@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef DCS_CONTROL_DESIGN_MATLAB_LINEAR_MPC_HPP
#define DCS_CONTROL_DESIGN_MATLAB_LINEAR_MPC_HPP


#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublasx/operation/empty.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <cmath>
#include <cstddef>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <dcs/control/detail/matlab_utility.hpp>
#include <dcs/system/posix_process.hpp>
#include <dcs/exception.hpp>
#include <dcs/macro.hpp>
#include <sstream>
#include <stdexcept>


namespace dcs { namespace control {

namespace matlab_lmpc_detail {

template <typename RealT>
class matlab_output_consumer
{
    public: matlab_output_consumer()
	: ok_(true)
    {
    }

    public: void operator()(dcs::system::posix_process& matlab_process)
    {
        ok_ = true;

        if (!matlab_process.alive())
        {
			ok_ = false;
			errmsg_ = "MATLAB is not running";
            return;
        }

        std::istream& is = matlab_process.output_stream();

//DCS_DEBUG_TRACE("BEGIN parsing MATLAB output");//XXX

        bool parse_line = false;
        while (matlab_process.alive() && is.good())
        {
//DCS_DEBUG_TRACE("Going to read another line: GOOD: " << is.good() << " - EOF: " << is.eof() << " - FAIL: " << is.fail() << " - BAD: " << is.bad() << " - !(): " << !static_cast<bool>(is) << "");//XXX
            std::string line;

            std::getline(is, line);

//DCS_DEBUG_TRACE("Read from MATLAB --> " << line);//XXX

            if (line.find("???") != std::string::npos
                || line.find("Error:") != std::string::npos)
            {
                DCS_DEBUG_TRACE("An error is occurred while executing MATLAB: " << line);//XXX
                ok_ = false;
				errmsg_ = line;
                break;
            }

            if (parse_line)
            {
                if (line.find("[/dcs::control::matlab_linear_mpc]") != std::string::npos)
                {
                    // The end of parsable lines
                    parse_line = false;
                }
                else
                {
                    typename std::string::size_type pos;

                    if ((pos = line.find("QPCode=")) != std::string::npos)
                    {
                        qpcode_ = line.substr(pos+7);
            			boost::to_lower(qpcode_);
//DCS_DEBUG_TRACE("Parsed as QPCode=" << qpcode_);//XXX
                    }
                    else if ((pos = line.find("Iterations=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+11), iterations_);
//DCS_DEBUG_TRACE("Parsed as Iterations=" << iterations_);//XXX
                    }
                    else if ((pos = line.find("cost=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+5), cost_);
//DCS_DEBUG_TRACE("Parsed as cost=" << cost_);//XXX
                    }
                    else if ((pos = line.find("uopt=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+5), uopt_);
//DCS_DEBUG_TRACE("Parsed as uopt=" << uopt_);//XXX
                    }
                    else if ((pos = line.find("Xopt=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+5), Xopt_);
//DCS_DEBUG_TRACE("Parsed as Xopt=" << Xopt_);//XXX
                    }
                    else if ((pos = line.find("Uopt=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+5), Uopt_);
//DCS_DEBUG_TRACE("Parsed as Uopt=" << Uopt_);//XXX
                    }
                    else if ((pos = line.find("Xopt=")) != std::string::npos)
                    {
                        detail::parse_matlab_str(line.substr(pos+5), Yopt_);
//DCS_DEBUG_TRACE("Parsed as Yopt=" << Yopt_);//XXX
                    }
                }
            }
            else
            {
                if (line.find("[dcs::control::matlab_linear_mpc]") != std::string::npos)
                {
                    // The beginning of parsable lines
                    parse_line = true;
                }
            }
//DCS_DEBUG_TRACE("Looping: GOOD: " << is.good() << " - EOF: " << is.eof() << " - FAIL: " << is.fail() << " - BAD: " << is.bad() << " - !(): " << !static_cast<bool>(is) << "");//XXX
        }

//DCS_DEBUG_TRACE("DONE WITH LOOPING");//XXX

//DCS_DEBUG_TRACE("END parsing MATLAB output");//XXX
//DCS_DEBUG_TRACE("IS state: " << is.good() << " - " << is.eof() << " - " << is.fail() << " - " << is.bad());//XXX
    }

	public: bool ok_;
	public: std::string errmsg_;
	public: boost::numeric::ublas::vector<RealT> uopt_;
    public: boost::numeric::ublas::matrix<RealT> Xopt_;
    public: boost::numeric::ublas::matrix<RealT> Uopt_;
    public: boost::numeric::ublas::matrix<RealT> Yopt_;
	public: std::string qpcode_;
	public: std::size_t iterations_;
	public: RealT cost_;
}; // matlab_output_consumer

} // Namespace matlab_lmpc_detail


/**
 * \brief Model Predictive Control (MPC) for linear systems (Maciejowski,2000).
 *
 * The linear Model Predictive Control (MPC) framework can be stated as follows.
 * Given a linear time-invariant (LTI) system:
 * \f[
 *  x(k+1) = Ax(k) + Bu(k),\\
 *  y(k) = Cx(k).
 * \f]
 * the MPC controller finds the best sequence of control actions \f$u(k|k),\ldots,u(k+H_c-1|k)\f$ by solving the following optimization problem:
 * \f{align}
 *  \text{minimize} \quad J(k) =\sum_{i=1}^{H_p} \lVert y(k+i|k)-y_{\text{ref}}(k+i) \rVert^2_{W_y(i)} + \sum_{i=0}^{H_c-1} \lVert \Delta u(k+i|k) \rVert^2_{W_{\Delta u}(i)},\\
 *  \text{subject to} \quad y_{\text{min}} \le y \ge y_{\text{max}},\\
 *                    \quad \Delta y_{\text{min}} \ge \Delta y \le \Delta y_{\text{max}},\\
 *                    \quad u_{\text{min}} \le u \ge u_{\text{max}},\\
 *                    \quad \Delta u_{\text{min}} \le \Delta u \ge \Delta u_{\text{max}}.
 * \f}
 * with initial condition:
 * \f[
 *  u(k-1|k) = u(k-1),\\
 *  x(k|k) = x(k).
 * \f]
 * where:
 * - The notation \f$\lVert x \rVert^2_Q\f$ represents the quadratic form \f$x^T Q x\f$, where \f$x\f$ is a vector and \f$Q\f$ is a symmetric matrix.
 * - \f$H_p\f$ is the prediction horizon
 * - \f$H_c\f$ is the control horizon
 * - \f$x(k+1|k),u(k+1|k),y(k+i|k)\f$ are the predicted system state, system input and system output at the \f$(k+i)\f$-th prediction horizon step, respectively, which are based on measurements up to time \f$k\f$. Note that \f$u(k-1|k)=u(k-1)\f$.
 * Note, we always assume that \f$H_c \le H_p\f$ and that \f$\Delta u(k+i|k)=0\f$ for \f$i \ge H_c\f$, so that \f$u(k+i|k)=u(k+H_c-1|k)\f$ for all \f$i \ge H_c\f$.
 *
 * \note Soft constraints and manipulated variable tracking are not supported yet.
 *
 * References
 * -# J.M. Maciejowski. "Predictive Control with Constraints," Prentice-Hall, 2000.
 * .
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 */
template <typename RealT>
class matlab_linear_mpc_controller
{
	public: typedef RealT real_type;
	public: typedef boost::numeric::ublas::matrix<real_type, boost::numeric::ublas::column_major> matrix_type;
	public: typedef boost::numeric::ublas::vector<real_type> vector_type;


	public: matlab_linear_mpc_controller()
	{
	}

	public: template <typename WyMatrixT,
					  typename WduMatrixT,
					  typename YMinVectorT,
					  typename YMaxVectorT,
					  typename DeltaYMinVectorT,
					  typename DeltaYMaxVectorT,
					  typename UMinVectorT,
					  typename UMaxVectorT,
					  typename DeltaUMinVectorT,
					  typename DeltaUMaxVectorT>
			matlab_linear_mpc_controller(boost::numeric::ublas::matrix_expression<WyMatrixT> const& Wy,
										 boost::numeric::ublas::matrix_expression<WduMatrixT> const& Wdu,
										 boost::numeric::ublas::vector_expression<YMinVectorT> const& ymin,
										 boost::numeric::ublas::vector_expression<YMaxVectorT> const& ymax,
										 boost::numeric::ublas::vector_expression<DeltaYMinVectorT> const& dymin,
										 boost::numeric::ublas::vector_expression<DeltaYMaxVectorT> const& dymax,
										 boost::numeric::ublas::vector_expression<UMinVectorT> const& umin,
										 boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
										 boost::numeric::ublas::vector_expression<DeltaUMinVectorT> const& dumin,
										 boost::numeric::ublas::vector_expression<DeltaUMaxVectorT> const& dumax,
										 std::size_t Hp,
										 std::size_t Hc)
	: Wy_(Wy),
	  Wdu_(Wdu),
	  ymin_(ymin),
	  ymax_(ymax),
	  dymin_(dymin),
	  dymax_(dymax),
	  umin_(umin),
	  umax_(umax),
	  dumin_(dumin),
	  dumax_(dumax),
	  Hp_(Hp),
	  Hc_(Hc)
	{
		namespace ublasx = boost::numeric::ublasx;

		const std::size_t ny = ublasx::size(ymin);
		const std::size_t nu = ublasx::size(umin);

		// pre: Hc_ > 0
		DCS_ASSERT(Hc_ > 0,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Control horizon cannot be zero"));
		// pre: Hp_ >= Hc_
		DCS_ASSERT(Hp_ >= Hc_,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Prediction horizon cannot be less than control horizon"));
		// pre: num_rows(Wy_) >= ny
		DCS_ASSERT(ublasx::num_rows(Wy_) >= ny,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too few rows for output weight matrix Wy"));
		// pre: num_rows(Wy_) <= ny*Hp_
		DCS_ASSERT(ublasx::num_rows(Wy_) <= ny*Hp_,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too many rows for output weight matrix Wy"));
		// pre: num_columns(Wy_) >= ny
		DCS_ASSERT(ublasx::num_columns(Wy_) >= ny,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too few columns for output weight matrix Wy"));
		// pre: num_columns(Wy_) <= ny*Hp_
		DCS_ASSERT(ublasx::num_columns(Wy_) <= ny*Hp_,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too many columns for output weight matrix Wy"));
		// pre: num_rows(Wdu_) >= nu_
		DCS_ASSERT(ublasx::num_rows(Wdu_) >= nu,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too few rows for manipulated variable weight matrix Wdu"));
		// pre: num_rows(Wdu_) <= nu*Hc_
		DCS_ASSERT(ublasx::num_rows(Wdu_) <= nu*Hc_,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too many rows for manipulated variable weight matrix Wdu"));
		// pre: num_columns(Wdu_) >= nu
		DCS_ASSERT(ublasx::num_columns(Wdu_) >= nu,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too few columns for manipulated variable weight matrix Wdu"));
		// pre: num_columns(Wdu_) <= nu*Hc_
		DCS_ASSERT(ublasx::num_columns(Wdu_) <= nu*Hc_,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too many columns for manipulated variable weight matrix Wdu"));
		// pre: size(ymin) == size(ymax)
		DCS_ASSERT(ublasx::size(ymin) == ublasx::size(ymax),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors ymin and ymax have incompatible size"));
		// pre: size(ymin) == size(dymin)
		DCS_ASSERT(ublasx::size(ymin) == ublasx::size(dymin),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors ymin and dymin have incompatible size"));
		// pre: size(ymin) == size(dymax)
		DCS_ASSERT(ublasx::size(ymin) == ublasx::size(dymax),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors ymin and dymax have incompatible size"));
		// pre: size(umin) == size(umax)
		DCS_ASSERT(ublasx::size(umin) == ublasx::size(umax),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors umin and umax have incompatible size"));
		// pre: size(umin) == size(dumin)
		DCS_ASSERT(ublasx::size(umin) == ublasx::size(dumin),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors umin and dumin have incompatible size"));
		// pre: size(umin) == size(dumax)
		DCS_ASSERT(ublasx::size(umin) == ublasx::size(dumax),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Vectors umin and dumax have incompatible size"));
	}

	public: void reset()
	{
		Wy_.clear();
		Wdu_.clear();
		ymin_.clear();
		ymax_.clear();
		dymin_.clear();
		dymax_.clear();
		umin_.clear();
		umax_.clear();
		dumin_.clear();
		dumax_.clear();
		Hp_ = 0;
		Hc_ = 0;
		A_.clear();
		B_.clear();
		C_.clear();
		Xopt_.clear();
		Uopt_.clear();
		Yopt_.clear();
		//J_ = std::numeric_limits<RealT>::infinity();
	}

	public: template <typename AMatrixT,
					  typename BMatrixT,
					  typename CMatrixT>
			void solve(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					   boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					   boost::numeric::ublas::matrix_expression<CMatrixT> const& C)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		const std::size_t nx = ublasx::num_rows(A);
		const std::size_t nu = ublasx::num_columns(B);
		const std::size_t ny = ublasx::num_rows(C);

		// pre: A is a square matrix
		DCS_ASSERT(nx == ublasx::num_columns(A),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A must be square"));
		// pre: num_rows(B) == num_rows(A)
		DCS_ASSERT(nx == ublasx::num_rows(B),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix B has a wrong number of rows"));
		// pre: num_columns(C) == num_rows(A)
		DCS_ASSERT(nx == ublasx::num_columns(C),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix C has a wrong number of columns"));
		// pre: num_rows(C) == ny
		DCS_ASSERT(ny == ublasx::num_rows(C),
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix C has a wrong number of rows"));

		A_ = A;
		B_ = B;
		C_ = C;

		if (ublasx::empty(Wy_))
		{
			Wy_ = ublas::identity_matrix<real_type>(ny);
		}
		if (ublasx::empty(Wdu_))
		{
			Wdu_= ublas::zero_matrix<real_type>(nu);
		}
//		if (ublasx::empty(Wdu_))
//		{
//			Wdu_= ublas::zero_matrix<real_type>(nu);
//		}
		if (ublasx::empty(umin_))
		{
			umin_= ublas::scalar_vector<real_type>(nu, -std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(umax_))
		{
			umax_= ublas::scalar_vector<real_type>(nu, +std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(dumin_))
		{
			dumin_= ublas::scalar_vector<real_type>(nu, -std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(dumax_))
		{
			dumax_= ublas::scalar_vector<real_type>(nu, +std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(ymin_))
		{
			ymin_= ublas::scalar_vector<real_type>(ny, -std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(ymax_))
		{
			ymax_= ublas::scalar_vector<real_type>(ny, +std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(dymin_))
		{
			dymin_= ublas::scalar_vector<real_type>(ny, -std::numeric_limits<real_type>::infinity());
		}
		if (ublasx::empty(dymax_))
		{
			dymax_= ublas::scalar_vector<real_type>(ny, +std::numeric_limits<real_type>::infinity());
		}

		// Spans the weighting matrices over the horizons
		if (ublasx::num_rows(Wdu_) == nu)
		{
			ublas::matrix<real_type> WduHc = ublas::zero_matrix<real_type>(nu*Hc_);
			for (std::size_t i = 0; i < Hc_; ++i)
			{
				const std::size_t i0 = i*nu;
				const std::size_t i1 = i0+nu;
				ublas::subrange(WduHc, i0, i1, i0, i1) = Wdu_;
			}
			Wdu_ = WduHc;
		}
		if (ublasx::num_rows(Wy_) == ny)
		{
			ublas::matrix<real_type> WyHp = ublas::zero_matrix<real_type>(ny*Hp_);
			for (std::size_t i = 0; i < Hp_; ++i)
			{
				const std::size_t i0 = i*ny;
				const std::size_t i1 = i0+ny;
				ublas::subrange(WyHp, i0, i1, i0, i1) = Wy_;
			}
			Wy_ = WyHp;
		}
	}

	public: template <typename XVectorT,
					  typename UVectorT,
					  typename YRefMatrixT>
			vector_type control(boost::numeric::ublas::vector_expression<XVectorT> const& x,
								boost::numeric::ublas::vector_expression<UVectorT> const& u,
								boost::numeric::ublas::matrix_expression<YRefMatrixT> const& Yref)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		const std::size_t nx = ublasx::num_rows(A_);
		const std::size_t nu = ublasx::num_columns(B_);
		const std::size_t ny = ublasx::num_rows(C_);

		DCS_ASSERT(ublasx::size(x) == nx,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Bad size for vector x"));
		DCS_ASSERT(ublasx::size(u) == nu,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Bad size for vector u"));
		DCS_ASSERT(ublasx::num_columns(Yref) == ny,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Bad number of columns for matrix Yref"));
		//if (ublasx::num_rows(Yref) < Hp_)
		DCS_ASSERT(ublasx::num_rows(Yref) > 0,
				   DCS_EXCEPTION_THROW(std::invalid_argument, "Too few rows for matrix Yref"));

		//FIXME: in MATLAB it's not possible to set constraints on the output variables rate (variables dymin_ and dymax_)

		Xopt_.resize(0, 0, false);
		Uopt_.resize(0, 0, false);
		Yopt_.resize(0, 0, false);
		//J_ = std::numeric_limits<RealT>::infinity();

		std::vector<std::string> matlab_args;

		matlab_args.push_back("-nodisplay");
		//matlab_args.push_back("-nojvm");
		matlab_args.push_back("-nodesktop");

		std::ostringstream oss;
		oss << "-r \""
			<< " use_custom_state_est = 1;" //FIXME: this disable the MATLAB's default state estimation (i.e., steady-state kalman filter) in favor of custom state estimation
			<< " try"
			<< "  nu = " << nu << ";"
			<< "  ny = " << ny << ";"
			<< "  A = " << detail::to_matlab_str(A_) << ";"
			<< "  B = " << detail::to_matlab_str(B_) << ";"
			<< "  C = " << detail::to_matlab_str(C_) << ";"
			<< "  D = zeros(ny,nu);"
			<< "  sys = ss(A,B,C,D,1);"
			<< "  Hp = " << Hp_ << ";"
			<< "  Hc = " << Hc_ << ";"
			<< "  ctrl = mpc(sys, sys.Ts, Hp, Hc);"
			<< "  if use_custom_state_est,"
			<< "   setEstimator(ctrl,'custom');"
			<< "  end;"
			<< "  ctrl.Weights.MVRate = reshape(diag(" << detail::to_matlab_str(Wdu_) << "), nu, Hc)';"
			<< "  ctrl.Weights.OV = reshape(diag(" << detail::to_matlab_str(Wy_) << "), ny, Hp)';";
		for (std::size_t i = 0; i < nu; ++i)
		{
			oss << "  ctrl.MV(" << (i+1) << ").Min = " << umin_(i) << ";"
				<< "  ctrl.MV(" << (i+1) << ").Max = " << umax_(i) << ";"
				<< "  ctrl.MV(" << (i+1) << ").RateMin = " << dumin_(i) << ";"
				<< "  ctrl.MV(" << (i+1) << ").RateMax = " << dumax_(i) << ";";
		}
		for (std::size_t i = 0; i < ny; ++i)
		{
			oss << "  ctrl.OV(" << (i+1) << ").Min = " << ymin_(i) << ";"
				<< "  ctrl.OV(" << (i+1) << ").Max = " << ymax_(i) << ";";
		}
		oss << "  x = " << detail::to_matlab_str(x) << ";"
			<< "  u = " << detail::to_matlab_str(u) << ";"
			<< "  y = sys.C*x+sys.D*u;"
			<< "  r = " << detail::to_matlab_str(Yref) << ";"
			<< "  ctrlstate = mpcstate(ctrl, x, [], [], u);"
			<< "  if use_custom_state_est,"
			<< "   ctrlstate.Plant = x;"
			<< "  end;"
			<< "  [uopt, info] = mpcmove(ctrl, ctrlstate, y, r);"
			<< "  format long;"
			<< "  disp('--- [dcs::control::matlab_linear_mpc] ---');"
			<< "  disp(['uopt=', mat2str(uopt)]);"
			<< "  disp(['Uopt=', mat2str(info.Uopt)]);"
			<< "  disp(['Yopt=', mat2str(info.Yopt)]);"
			<< "  disp(['Xopt=', mat2str(info.Xopt)]);"
			<< "  disp(['QPCode=', info.QPCode]);"
			<< "  disp(['Iterations=', num2str(info.Iterations)]);"
			<< "  disp(['Cost=', num2str(info.Cost)]);"
			<< "  disp('--- [/dcs::control::matlab_linear_mpc] ---');"
			<< " catch me,"
			<< "  disp(['??? Error: ', me.message]);"
			<< " end;"
			<< " quit force;"
			<< "\"";

		matlab_args.push_back(oss.str());

		matlab_lmpc_detail::matlab_output_consumer<RealT> matlab_consumer;
		detail::run_matlab("matlab", matlab_args.begin(), matlab_args.end(), matlab_consumer);

		if (!matlab_consumer.ok_)
		{
			DCS_EXCEPTION_THROW(std::runtime_error, matlab_consumer.errmsg_);
		}

		Xopt_ = matlab_consumer.Xopt_;
		Uopt_ = matlab_consumer.Uopt_;
		Yopt_ = matlab_consumer.Yopt_;
		//J_ = matlab_consumer.cost_;

		//DCS_DEBUG_TRACE("CONTROL -- predicted X: " << Xopt_);//XXX
		//DCS_DEBUG_TRACE("CONTROL -- optimal U: " << Uopt_);//XXX
		//DCS_DEBUG_TRACE("CONTROL -- predicted Y: " << Yopt_);//XXX

		return ublas::row(Uopt_, 0);
	}

	public: template <typename XVectorT,
					  typename UVectorT,
					  typename YRefVectorT>
			vector_type control(boost::numeric::ublas::vector_expression<XVectorT> const& x,
								boost::numeric::ublas::vector_expression<UVectorT> const& u,
								boost::numeric::ublas::vector_expression<YRefVectorT> const& yref)
	{
		namespace ublas = boost::numeric::ublas;
		namespace ublasx = boost::numeric::ublasx;

		matrix_type Yref(Hp_, boost::numeric::ublasx::size(yref));
		for (std::size_t i = 0; i < Hp_; ++i)
		{
			ublas::row(Yref, i) = yref;
		}
		return this->control(x, u, Yref);
	}

	// Mimics the MATLAB's sim function
	public: template <typename XVectorT,
					  typename UVectorT,
					  typename YRefMatrixT>
			matrix_type simulate(boost::numeric::ublas::vector_expression<XVectorT> const& x,
								 boost::numeric::ublas::vector_expression<UVectorT> const& u,
								 boost::numeric::ublas::matrix_expression<YRefMatrixT> const& Yref)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		const std::size_t nx = ublasx::num_rows(A_);
		const std::size_t nu = ublasx::num_columns(B_);
		const std::size_t ny = ublasx::num_rows(C_);
		const std::size_t nsim = ublasx::num_rows(Yref);

		if (ublasx::size(x) != nx)
		{
			DCS_EXCEPTION_THROW(std::invalid_argument, "Bad size for vector x");
		}
		if (ublasx::size(u) != nu)
		{
			DCS_EXCEPTION_THROW(std::invalid_argument, "Bad size for vector u");
		}
		if (ublasx::num_columns(Yref) != ny)
		{
			DCS_EXCEPTION_THROW(std::invalid_argument, "Bad number of columns for matrix Yref");
		}
		//if (ublasx::num_rows(Yref) < Hp_)
		//{
		//	DCS_EXCEPTION_THROW(std::invalid_argument, "Too few rows for matrix Yref");
		//}

		ublas::vector<real_type> xx = x;
		ublas::vector<real_type> uu = u;
		ublas::matrix<real_type> U(nsim, nu, 0);
#ifdef DCS_DEBUG
		ublas::matrix<real_type> X(nsim, nx, 0);//XXX
		ublas::matrix<real_type> Y(nsim, ny, 0);//XXX
#endif // DCS_DEBUG
		for (std::size_t k = 0; k < nsim; ++k)
		{
			this->control(xx, uu, ublas::row(Yref(), k));
			xx = ublas::row(Xopt_, 1);
			uu = ublas::row(Uopt_, 0);
			ublas::row(U, k) = uu;
#ifdef DCS_DEBUG
			//ublas::row(X, k) = xx; //XXX
			//ublas::row(Y, k) = ublas::row(Yopt_, 1); //XXX
			ublas::row(X, k) = ublas::row(Xopt_, 0); //XXX
			ublas::row(Y, k) = ublas::row(Yopt_, 0); //XXX
#endif // DCS_DEBUG
		}

		//DCS_DEBUG_TRACE("SIMULATION -- Yref: " << Yref);//XXX
		//DCS_DEBUG_TRACE("SIMULATION -- X: " << X);//XXX
		//DCS_DEBUG_TRACE("SIMULATION -- U: " << U);//XXX
		//DCS_DEBUG_TRACE("SIMULATION -- Y: " << Y);//XXX

		return U;
	}

	public: matrix_type Wy() const
	{
		return Wy_;
	}

	public: matrix_type Wdu() const
	{
		return Wdu_;
	}

	public: vector_type ymin() const
	{
		return ymin_;
	}

	public: vector_type ymax() const
	{
		return ymax_;
	}

	public: vector_type dymin() const
	{
		return dymin_;
	}

	public: vector_type dymax() const
	{
		return dymax_;
	}

	public: vector_type umin() const
	{
		return umin_;
	}

	public: vector_type umax() const
	{
		return umax_;
	}

	public: vector_type dumin() const
	{
		return dumin_;
	}

	public: vector_type dumax() const
	{
		return dumax_;
	}

	public: std::size_t prediction_horizon() const
	{
		return Hp_;
	}

	public: std::size_t control_horizon() const
	{
		return Hc_;
	}

	public: matrix_type predicted_states() const
	{
		return Xopt_;
	}

	public: matrix_type predicted_outputs() const
	{
		return Yopt_;
	}

	public: matrix_type optimal_inputs() const
	{
		return Uopt_;
	}

//TODO: we need to collect the objective values of the QP problem during each iteration in the control method
//	private: real_type cost() const
//	{
//		return J_;
//	}


	private: matrix_type Wy_; ///< Weight matrix for the output variables
	private: matrix_type Wdu_; ///< Weight matrix for the manipulated variable adjustments (moves)
	private: vector_type ymin_; ///< Contraint vector definining the lower bounds of the output variables
	private: vector_type ymax_; ///< Contraint vector definining the upper bounds of the output variables
	private: vector_type dymin_; ///< Contraint vector definining the lower bounds of the manipulated variable adjustments
	private: vector_type dymax_; ///< Contraint vector definining the upper bounds of the manipulated variable adjustments
	private: vector_type umin_; ///< Contraint vector definining the lower bounds of the manipulated variables
	private: vector_type umax_; ///< Contraint vector definining the upper bounds of the manipulated variables
	private: vector_type dumin_; ///< Contraint vector definining the lower bounds of the manipulated variable adjustments
	private: vector_type dumax_; ///< Contraint vector definining the upper bounds of the manipulated variable adjustments
	private: std::size_t Hp_; ///< Prediction horizon
	private: std::size_t Hc_; ///< Control horizon
	private: matrix_type A_; ///< The state matrix of the linear system
	private: matrix_type B_; ///< The input matrix of the linear system
	private: matrix_type C_; ///< The output matrix of the linear system
	//private: matrix_type Ae_; ///< The state matrix of the extended linear system
	//private: matrix_type Be_; ///< The input matrix of the extended linear system
	//private: matrix_type Ce_; ///< The output matrix of the extended linear system
	//private: matrix_type Iu_;
	//private: matrix_type Ru_;
	//private: matrix_type Rx_;
	//private: matrix_type dRx_;
	//private: matrix_type Rx1_;
	//private: vector_type umax1_;
	//private: vector_type umin1_;
	//private: vector_type dumax1_;
	//private: vector_type dumin1_;
	//private: vector_type ymax1_;
	//private: vector_type ymin1_;
	//private: vector_type dymax1_;
	//private: vector_type dymin1_;
	//private: matrix_type H_; 
	//private: matrix_type Alpha_; 
	private: matrix_type Xopt_;
	private: matrix_type Uopt_;
	private: matrix_type Yopt_;
	//private: real_type J_;
}; // matlab_linear_mpc_controller


template <typename RealT,
		  typename AMatrixT,
		  typename BMatrixT,
		  typename CMatrixT,
		  typename WyMatrixT,
		  typename WduMatrixT,
		  typename YMinVectorT,
		  typename YMaxVectorT,
		  typename DeltaYMinVectorT,
		  typename DeltaYMaxVectorT,
		  typename UMinVectorT,
		  typename UMaxVectorT,
		  typename DeltaUMinVectorT,
		  typename DeltaUMaxVectorT,
		  typename XVectorT,
		  typename UVectorT,
		  typename YRefMatrixT>
boost::numeric::ublas::vector<RealT> matlab_linear_mpc(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
												boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
												boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
												boost::numeric::ublas::matrix_expression<WyMatrixT> const& Wy,
												boost::numeric::ublas::matrix_expression<WduMatrixT> const& Wdu,
												boost::numeric::ublas::vector_expression<YMinVectorT> const& ymin,
												boost::numeric::ublas::vector_expression<YMaxVectorT> const& ymax,
												boost::numeric::ublas::vector_expression<DeltaYMinVectorT> const& dymin,
												boost::numeric::ublas::vector_expression<DeltaYMaxVectorT> const& dymax,
												boost::numeric::ublas::vector_expression<UMinVectorT> const& umin,
												boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
												boost::numeric::ublas::vector_expression<DeltaUMinVectorT> const& dumin,
												boost::numeric::ublas::vector_expression<DeltaUMaxVectorT> const& dumax,
												std::size_t Hp,
												std::size_t Hc,
												boost::numeric::ublas::vector_expression<XVectorT> const& x,
												boost::numeric::ublas::vector_expression<UVectorT> const& u,
												boost::numeric::ublas::matrix_expression<YRefMatrixT> const& Yref)
{
	matlab_linear_mpc_controller<RealT> lmpc(Wy, Wdu, ymin, ymax, dymin, dymax, umin, umax, dumin, dumax, Hp, Hc);
	lmpc.solve(A, B, C);
	return lmpc.control(x, u, Yref);
}


template <typename RealT,
		  typename AMatrixT,
		  typename BMatrixT,
		  typename CMatrixT,
		  typename WyMatrixT,
		  typename WduMatrixT,
		  typename YMinVectorT,
		  typename YMaxVectorT,
		  typename DeltaYMinVectorT,
		  typename DeltaYMaxVectorT,
		  typename UMinVectorT,
		  typename UMaxVectorT,
		  typename DeltaUMinVectorT,
		  typename DeltaUMaxVectorT,
		  typename XVectorT,
		  typename UVectorT,
		  typename YRefVectorT>
boost::numeric::ublas::vector<RealT> matlab_linear_mpc(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
												boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
												boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
												boost::numeric::ublas::matrix_expression<WyMatrixT> const& Wy,
												boost::numeric::ublas::matrix_expression<WduMatrixT> const& Wdu,
												boost::numeric::ublas::vector_expression<YMinVectorT> const& ymin,
												boost::numeric::ublas::vector_expression<YMaxVectorT> const& ymax,
												boost::numeric::ublas::vector_expression<DeltaYMinVectorT> const& dymin,
												boost::numeric::ublas::vector_expression<DeltaYMaxVectorT> const& dymax,
												boost::numeric::ublas::vector_expression<UMinVectorT> const& umin,
												boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
												boost::numeric::ublas::vector_expression<DeltaUMinVectorT> const& dumin,
												boost::numeric::ublas::vector_expression<DeltaUMaxVectorT> const& dumax,
												std::size_t Hp,
												std::size_t Hc,
												boost::numeric::ublas::vector_expression<XVectorT> const& x,
												boost::numeric::ublas::vector_expression<UVectorT> const& u,
												boost::numeric::ublas::vector_expression<YRefVectorT> const& yref)
{
	matlab_linear_mpc_controller<RealT> lmpc(Wy, Wdu, ymin, ymax, dymin, dymax, umin, umax, dumin, dumax, Hp, Hc);
	lmpc.solve(A, B, C);
	return lmpc.control(x, u, yref);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_DESIGN_MATLAB_LINEAR_MPC_HPP
