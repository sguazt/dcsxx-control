/**
 * \file dcs/control/design/linear_mpc.hpp
 *
 * \brief Linear Model Predictive Control
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

#ifndef DCS_CONTROL_DESIGN_LINEAR_MPC_HPP
#define DCS_CONTROL_DESIGN_LINEAR_MPC_HPP


#include <algorithm>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/empty.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/pow.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <cmath>
#include <cstddef>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <dcs/control/solver/qp.hpp>
#include <dcs/exception.hpp>
#include <dcs/macro.hpp>
#include <stdexcept>


namespace dcs { namespace control {

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
class linear_mpc_controller
{
	public: typedef RealT real_type;
	public: typedef boost::numeric::ublas::matrix<real_type, boost::numeric::ublas::column_major> matrix_type;
	public: typedef boost::numeric::ublas::vector<real_type> vector_type;


	public: linear_mpc_controller()
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
			linear_mpc_controller(boost::numeric::ublas::matrix_expression<WyMatrixT> const& Wy,
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
		Ae_.clear();
		Be_.clear();
		Ce_.clear();
		Iu_.clear();
		Ru_.clear();
		Rx_.clear();
		dRx_.clear();
		Rx1_.clear();
		umax1_.clear();
		umin1_.clear();
		dumax1_.clear();
		dumin1_.clear();
		ymax1_.clear();
		ymin1_.clear();
		dymax1_.clear();
		dymin1_.clear();
		H_.clear();
		Alpha_.clear();
		Xopt_.clear();
		Uopt_.clear();
		Yopt_.clear();
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

		// Forms the augmented system 
		Ae_ = ublasx::cat_columns(ublasx::cat_rows(A_, B_),
								  ublasx::cat_rows(ublas::zero_matrix<real_type>(nu, nx),
												   ublas::identity_matrix<real_type>(nu)));
		Be_ = ublasx::cat_columns(B_, ublas::identity_matrix<real_type>(nu));
		Ce_ = ublasx::cat_rows(C_, ublas::zero_matrix<real_type>(ny, nu));

		const std::size_t nxe = ublasx::num_rows(Ae_);
		const std::size_t nue = ublasx::num_columns(Be_);
		const std::size_t nye = ublasx::num_rows(Ce_);

		// Auxiliary variables
		ublas::matrix<real_type> RRR = ublas::zero_matrix<real_type>(nye*Hp_, nue);
		for (std::size_t i = 0; i < Hp_; ++i)
		{
			const std::size_t r0 = i*nye;
			const std::size_t r1 = r0 + nye;
			ublas::matrix<real_type> tmp1 = ublas::prod(Ce_, ublasx::pow(Ae_, i));
			ublas::subrange(RRR, r0, r1, 0, nue) = ublas::prod(tmp1, Be_);
		}

		// Constructs the Ru and Rx matrices
		Ru_ = ublas::matrix<real_type>(nye*Hp_, nue*Hc_, 0);
		for (std::size_t i = 0; i < Hc_; ++i)
		{
			ublas::matrix<real_type> tmp = ublasx::cat_columns(ublas::zero_matrix<real_type>(i*nye, nue), ublas::subrange(RRR, 0, nye*(Hp_-i), 0, nue));
			ublas::subrange(Ru_, 0, nye*Hp_, i*nue, (i+1)*nue) = ublasx::cat_columns(ublas::zero_matrix<real_type>(i*nye, nue), ublas::subrange(RRR, 0, nye*(Hp_-i), 0, nue));
		}
		Rx_ = ublas::zero_matrix<real_type>(nye*Hp_, nxe);
		for (std::size_t i = 0; i < Hp_; ++i)
		{
			const std::size_t r0 = i*nye;
			const std::size_t r1 = r0 + nye;
			ublas::subrange(Rx_, r0, r1, 0, nxe) = ublas::prod(Ce_, ublasx::pow(Ae_, i));
		}

		// Auxiliary matrix
		// FIXME: RRR already contains part of info that is computed in the following loop
		//RRR = ublas::zero_matrix<real_type>(nye*Hp_, nue);
		//for (std::size_t i = 0; i < Hp_; ++i)
		//{
		//	ublas::subrange(RRR, i*nye, i*nye+nye, 0, nue) = ublas::prod(ublas::prod(Ce_, ublasx::pow(Ae_, i)), Be_);
		//	if (i > 0)
		//	{
		//		ublas::subrange(RRR, i*nye, i*nye+nye, 0, nue) -= ublas::prod(ublas::prod(Ce_, ublasx::pow(Ae_, i-1)), Be_);
		//	}
		//}
		RRR.resize((Hp_+2)*nye, nue);
		for (std::size_t i = 1; i < (Hp_+2); ++i)
		{
			const std::size_t r0 = i*nye;
			const std::size_t r1 = r0 + nye;
			if (i >= Hp_)
			{
				ublas::matrix<real_type> tmp1 = ublas::prod(Ce_, ublasx::pow(Ae_, i));
				ublas::subrange(RRR, r0, r1, 0, nue) = ublas::prod(tmp1, Be_);
			}
			//ublas::subrange(RRR, r0, r1, 0, nue) -= ublas::prod(ublas::prod(Ce_, ublasx::pow(Ae_, i-1)), Be_);
			ublas::matrix<real_type> tmp1 = ublas::prod(Ce_, ublasx::pow(Ae_, i-1));
			ublas::subrange(RRR, r0, r1, 0, nue) -= ublas::prod(tmp1, Be_);
		}

		// Constructs the dRu, dRx, dRu1, and dRx1 matrices for delta y constraints
		//ublas::matrix<real_type> dRu(ublasx::num_rows(RRR), ublasx::num_columns(RRR)*Hc_);
		ublas::matrix<real_type> dRu(ublasx::num_rows(RRR)-2*nye, ublasx::num_columns(RRR)*Hc_);
		ublas::subrange(dRu, 0, ublasx::num_rows(dRu), 0, ublasx::num_columns(RRR)) = ublas::subrange(RRR, 2*nye, ublasx::num_rows(RRR), 0, ublasx::num_columns(RRR));
		if (Hc_ > 1)
		{
			const std::size_t r0 = nye;
			const std::size_t r1 = ublasx::num_rows(RRR)-nye;
			const std::size_t c0 = ublasx::num_columns(RRR);
			const std::size_t c1 = c0 + ublasx::num_columns(RRR);
			ublas::subrange(dRu, 0, ublasx::num_rows(dRu), c0, c1) = ublas::subrange(RRR, r0, r1, 0, c0);
		}
		for (std::size_t i = 2; i < Hc_; ++i)
		{
			//const std::size_t c0 = 2*(i+1)*ublasx::num_columns(RRR);
			const std::size_t c0 = i*ublasx::num_columns(RRR);
			const std::size_t c1 = c0 + ublasx::num_columns(RRR);
			//if (i > 2)
			//{
			//	ublas::subrange(dRu, 0, ublasx::num_rows(dRu), c0, c1) = ublasx::cat_columns(ublas::zero_matrix<real_type>((i-2)*nye, nue), ublas::subrange(RRR, 0, nye*(Hp_+2-i), 0, ublasx::num_columns(RRR)));
			//}
			//else
			//{
			//	ublas::subrange(dRu, 0, ublasx::num_rows(dRu), c0, c1) = ublas::subrange(RRR, 0, nye*(Hp_+2-i), 0, ublasx::num_columns(RRR));
			//}
			ublas::subrange(dRu, 0, ublasx::num_rows(dRu), c0, c1) = ublasx::cat_columns(ublas::zero_matrix<real_type>((i-2)*nye, nue), ublas::subrange(RRR, 0, nye*(Hp_+2-i), 0, ublasx::num_columns(RRR)));
		}
		dRx_ = ublas::matrix<real_type>(nye*Hp_, nxe, 0);
		ublas::subrange(dRx_, 0, nye, 0, nxe) = ublas::prod(Ce_, Ae_) - Ce_;
		for (std::size_t i = 0; i < (Hp_-1); ++i)
		{
			const std::size_t r0 = (i+1)*nye;
			const std::size_t r1 = r0 + nye;
			ublas::subrange(dRx_, r0, r1, 0, nxe) = ublas::prod(Ce_, ublasx::pow(Ae_, i+1)) - ublas::prod(Ce_, ublasx::pow(Ae_, i));
		}
		ublas::matrix<real_type> Ru1 = ublas::subrange(Ru_, 0, nye, 0, ublasx::num_columns(Ru_));
		Rx1_ = ublas::subrange(Rx_, 0, nye, 0, ublasx::num_columns(Rx_));

		// $I_u$ matrix
		Iu_ = ublas::matrix<real_type>(nue*Hc_, nue, 0);
		for (std::size_t i = 0; i < Hc_; ++i)
		{
			const std::size_t r0 = i*nue;
			const std::size_t r1 = r0 + nue;
			ublas::subrange(Iu_, r0, r1, 0, nue) = ublas::identity_matrix<real_type>(nue);
		}

		// $I_{\Delta u}$ matrix
		ublas::matrix<real_type> Idu = ublas::zero_matrix<real_type>(nue*Hc_);
		for (std::size_t i = 0; i < Hc_; ++i)
		{
			const std::size_t r0 = i*nue;
			const std::size_t r1 = r0 + nue;

			for (std::size_t j = 0; j <= i; ++j)
			{
				const std::size_t c0 = j*nue;
				const std::size_t c1 = c0 + nue;
				ublas::subrange(Idu, r0, r1, c0, c1) = ublas::identity_matrix<real_type>(nue);
			}
		}

		// Additional variables needed for constraint handling
		ublas::matrix<real_type> Hp1 = ublas::scalar_matrix<real_type>(1, Hp_, 1.0);
		ublas::matrix<real_type> Hc1 = ublas::scalar_matrix<real_type>(1, Hc_, 1.0);
		ublas::matrix<real_type> ImHc = ublas::identity_matrix<real_type>(Hc_*nue);
		umax1_ = ublas::vector<real_type>(nu*Hc_);
		umin1_ = ublas::vector<real_type>(nu*Hc_);
		dumax1_ = ublas::vector<real_type>(nu*Hc_);
		dumin1_ = ublas::vector<real_type>(nu*Hc_);
		for (std::size_t i = 0; i < Hc_; ++i)
		{
			const std::size_t i0 = i*nu;
			const std::size_t i1 = i0 + nu;
			ublas::subrange(umax1_, i0, i1) = umax_;
			ublas::subrange(umin1_, i0, i1) = umin_;
			ublas::subrange(dumax1_, i0, i1) = dumax_;
			ublas::subrange(dumin1_, i0, i1) = dumin_;
		}
		ymax1_ = ublas::vector<real_type>(ny*Hp_, 0);
		ymin1_ = ublas::vector<real_type>(ny*Hp_, 0);
		dymax1_ = ublas::vector<real_type>(ny*Hp_, 0);
		dymin1_ = ublas::vector<real_type>(ny*Hp_, 0);
		for (std::size_t i = 0; i < Hp_; ++i)
		{
			const std::size_t i0 = i*ny;
			const std::size_t i1 = i0 + ny;
			ublas::subrange(ymax1_, i0, i1) = ymax_;
			ublas::subrange(ymin1_, i0, i1) = ymin_;
			ublas::subrange(dymax1_, i0, i1) = dymax_;
			ublas::subrange(dymin1_, i0, i1) = dymin_;
		}

		// $\Alpha$ matrix
		Alpha_ = ublas::matrix<real_type>(ublasx::num_rows(Idu)*2 + ublasx::num_rows(ImHc)*2 + ublasx::num_rows(Ru_)*2 + ublasx::num_rows(Ru1)*2 + ublasx::num_rows(dRu)*2, ublasx::num_columns(Idu), 0);
		{
			std::size_t r0 = 0;
			std::size_t r1 = ublasx::num_rows(Idu);
			std::size_t c0 = 0;
			std::size_t c1 = ublasx::num_columns(Idu);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = Idu;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(Idu);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = -Idu;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(ImHc);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = ImHc;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(ImHc);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = -ImHc;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(Ru_);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = Ru_;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(Ru_);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = -Ru_;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(Ru1);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = Ru1;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(dRu);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = dRu;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(Ru1);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = -Ru1;
			r0 = r1;
			r1 = r0 + ublasx::num_rows(dRu);
			ublas::subrange(Alpha_, r0, r1, c0, c1) = -dRu;
		}

		// $H$ matrix
		//ublas::matrix<real_type> H = 2*(ublas::prod(ublas::prod(ublas::trans(Ru), Wy_), Ru) + Wdu_);
		{
			ublas::matrix<real_type> tmp1 = ublas::prod(ublas::trans(Ru_), Wy_);
			H_ = 2*(ublas::prod(tmp1, Ru_) + Wdu_);
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

		const std::size_t nxe = ublasx::num_rows(Ae_);
		const std::size_t nue = ublasx::num_columns(Be_);

		ublas::vector<real_type> xe(nxe, 0);
		ublas::subrange(xe, 0, nx) = x;
		ublas::subrange(xe, nx, nxe) = u;
		ublas::vector<real_type> ye = ublas::prod(ublas::prod(Ce_, Ae_), xe);

		// Main simulation loop
		Uopt_ = ublas::zero_matrix<real_type>(Hp_, nu);
		Xopt_ = ublas::zero_matrix<real_type>(Hp_+1, nx);
		ublas::row(Xopt_, 0) = x;
		Yopt_ = ublas::zero_matrix<real_type>(Hp_+1, ny);
		ublas::row(Yopt_, 0) = ublas::prod(C_, x);
		ublas::vector<real_type> uu = u;

		ublas::vector<real_type> w(ny*Hp_);
		for (std::size_t j = 0; j < Hp_; ++j)
		{
			ublas::subrange(w, j*ny, (j+1)*ny) = ublas::row(Yref(), std::min(j, ublasx::num_rows(Yref)-1));
		}

		// $\omega$ vector
		ublas::vector<real_type> omega(ublasx::size(umax1_)*4+ublasx::size(dymax_)*2+ublasx::size(ymax1_)*4, 0);
		{
			std::size_t r0 = 0;
			std::size_t r1 = ublasx::size(umax1_);
			ublas::subrange(omega, r0, r1) =  umax1_ - ublas::prod(Iu_, uu);
			r0 = r1;
			r1 = r0 + ublasx::size(umin1_);
			ublas::subrange(omega, r0, r1) = -umin1_ + ublas::prod(Iu_, uu);
			r0 = r1;
			r1 = r0 + ublasx::size(dumax1_);
			ublas::subrange(omega, r0, r1) =  dumax1_;
			r0 = r1;
			r1 = r0 + ublasx::size(dumin1_);
			ublas::subrange(omega, r0, r1) = -dumin1_;
			r0 = r1;
			r1 = r0 + ublasx::size(ymax1_);
			ublas::subrange(omega, r0, r1) =  ymax1_ - ublas::prod(ublas::prod(Rx_, Ae_), xe);
			r0 = r1;
			r1 = r0 + ublasx::size(ymin1_);
			ublas::subrange(omega, r0, r1) = -ymin1_ + ublas::prod(ublas::prod(Rx_, Ae_), xe);
			r0 = r1;
			r1 = r0 + ublasx::size(dymax_);
			ublas::subrange(omega, r0, r1) =  dymax_ - ublas::prod(ublas::prod(Rx1_, Ae_), xe) + ye;
			r0 = r1;
			r1 = r0 + ublasx::size(dymax1_);
			ublas::subrange(omega, r0, r1) =  dymax1_ - ublas::prod(ublas::prod(dRx_, Ae_), xe);
			r0 = r1;
			r1 = r0 + ublasx::size(dymin_);
			ublas::subrange(omega, r0, r1) = -dymin_ + ublas::prod(ublas::prod(Rx1_, Ae_), xe) - ye;
			r0 = r1;
			r1 = r0 + ublasx::size(dymax1_);
			ublas::subrange(omega, r0, r1) = -dymin1_ + ublas::prod(ublas::prod(dRx_, Ae_), xe);
		}

		// $c$ vector
		ublas::vector<real_type> c = 2*(ublas::prod(ublas::prod(ublas::trans(Ru_), ublas::trans(Wy_)), ublas::prod(ublas::prod(Rx_, Ae_), xe) - w));

		// Computes the control output
//DCS_DEBUG_TRACE("H: " << H_);//XXX
//DCS_DEBUG_TRACE("c: " << c);//XXX
//DCS_DEBUG_TRACE("Alpha: " << Alpha_);//XXX
//DCS_DEBUG_TRACE("Omega: " << omega);//XXX
		ublas::vector<real_type> du_opt = qp_solve<real_type>(H_, c, Alpha_, omega);
//DCS_DEBUG_TRACE("dU*: " << du_opt);//XXX

		// Integrates u
		for (std::size_t k = 0; k < Hp_; ++k)
		{
			ublas::vector<real_type> du = (k < Hc_) ? ublas::subrange(du_opt, k*nue, (k+1)*nue) : ublas::subrange(du_opt, (Hc_-1)*nue, Hc_*nue);
			uu += du;
			ublas::row(Uopt_, k) = uu;

			xe = ublas::prod(Ae_, xe) + ublas::prod(Be_, du);
			ye = ublas::prod(Ce_, xe);

			ublas::row(Xopt_, k+1) = ublas::prod(A_, ublas::row(Xopt_,k)) + ublas::prod(B_, uu);
			ublas::row(Yopt_, k+1) = ublas::prod(C_, ublas::row(Xopt_,k+1));
		}

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
	private: matrix_type Ae_; ///< The state matrix of the extended linear system
	private: matrix_type Be_; ///< The input matrix of the extended linear system
	private: matrix_type Ce_; ///< The output matrix of the extended linear system
	private: matrix_type Iu_;
	private: matrix_type Ru_;
	private: matrix_type Rx_;
	private: matrix_type dRx_;
	private: matrix_type Rx1_;
	private: vector_type umax1_;
	private: vector_type umin1_;
	private: vector_type dumax1_;
	private: vector_type dumin1_;
	private: vector_type ymax1_;
	private: vector_type ymin1_;
	private: vector_type dymax1_;
	private: vector_type dymin1_;
	private: matrix_type H_; 
	private: matrix_type Alpha_; 
	private: matrix_type Xopt_;
	private: matrix_type Uopt_;
	private: matrix_type Yopt_;
	//private: real_type J_;
}; // linear_mpc_controller


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
boost::numeric::ublas::vector<RealT> linear_mpc(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
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
	linear_mpc_controller<RealT> lmpc(Wy, Wdu, ymin, ymax, dymin, dymax, umin, umax, dumin, dumax, Hp, Hc);
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
boost::numeric::ublas::vector<RealT> linear_mpc(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
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
	linear_mpc_controller<RealT> lmpc(Wy, Wdu, ymin, ymax, dymin, dymax, umin, umax, dumin, dumax, Hp, Hc);
	lmpc.solve(A, B, C);
	return lmpc.control(x, u, yref);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_DESIGN_LINEAR_MPC_HPP
