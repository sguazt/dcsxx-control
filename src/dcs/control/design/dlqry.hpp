/**
 * \file dcs/control/design/dlqry.hpp
 *
 * \brief Infinite-horizon Discrete Linear Quadratic state-feedback Regulator
 *  with output weighting.
 *
 * Given a dicrete time-invariant state-space system model:
 * \f{align*}{
 *   \mathbf{x}(k+1) &= \mathbf{A}\mathbf{x}(k)+\mathbf{B}\mathbf{u}(k) \\
 *   \mathbf{y}(k+1) &= \mathbf{C}\mathbf{x}(k+1)+\mathbf{D}\mathbf{x}(k+1)
 * \f}
 * subject to initial conditions:
 * \f[
 *   \mathbf{x}(0) = \mathbf{x}_0
 * \f]
 * where \f$\mathbf{A} \in \mathbb{R}^{n \times n}\f$,
 *       \f$\mathbf{B} \in \mathbb{R}^{n \times m}\f$,
 *       \f$\mathbf{C} \in \mathbb{R}^{p \times n}\f$, and
 *       \f$\mathbf{D} \in \mathbb{R}^{p \times m}\f$,
 * and assuming that \f$(A,B)\f$ is \e stabilizable,
 * the <em>Infinite-horizon, Discrete Linear Quadratic (LQ) state-feedback
 * Regulator</em> controller problem calculates the optimal \e gain matrix
 * \f$\mathbf{K}\f$ such that the state-feedback law:
 * \f[
 *   \mathbf{u}(k) = -K\mathbf{x}(k)
 * \f]
 * minimizes the following quadratic cost function (<em>performance index</em>):
 * \f{align*}
 *   J(\mathbf{u}) &= \sum_{k=1}^{\infty}{\mathbf{y}^T(k)\mathbf{Q}\mathbf{y}(k) + \mathbf{u}^T(k)\mathbf{R}\mathbf{u}(k) + 2\mathbf{y}^T(k)\mathbf{N}\mathbf{u}(k)} \\
 *                 &= \sum_{k=1}^{\infty}{\begin{pmatrix}\mathbf{y}^T(k) & \mathbf{u}^T(k)\end{pmatrix}\begin{pmatrix}\mathbf{Q} & \mathbf{N} \\ \mathbf{N}^T & \mathbf{R}\end{pmatrix}\begin{pmatrix}\mathbf{y}(k) \\ \mathbf{u}(k)\end{pmatrix}}
 * \f}
 * where:
 * - the <em>error weighting matrix</em>
 *           \f$\mathbf{Q}\in\mathbb{R}^{n \times n}\f$ is a symmetric positive
 *           semi-definite real matrix that penalizes the output vector 
 *           \f$\mathbf{y}\f$ in the cost function,
 * - the <em>control weighting matrix</em>
 *           \f$\mathbf{R}\in\mathbb{R}^{m \times m}\f$ is a symmetric positive
 *           definite real matrix that penalizes the input vector
 *           \f$\mathbf{u}\f$ in the cost function, and
 * - the <em>cross-coupling weighting matrix</em>
 *           \f$\mathbf{N}\in\mathbb{R}^{n \times m}\f$ is a real matrix that
 *           penalizes the cross product between input and output vectors, such
 *           that \f$(\mathbf{Q}-\mathbf{N}\mathbf{R}^{-1}\mathbf{N}^{T})\f$ is
 *           positive semi-definite.
 * .
 *
 * The optimal gain matrix \f$\mathbf{K}\f$ is derived from:
 *  \f[
 *    \mathbf{K} = (\mathbf{B}^{T}\mathbf{X}\mathbf{B}+\mathbf{R})^{-1}(\mathbf{B}^{T}\mathbf{X}\mathbf{A}+\mathbf{N}^{T})
 *  \f]
 *  where \f$\mathbf{X}\f$ is the stabilizing infinite-horizon solution of the
 *  associated Discrete-time Algebraic Riccati Equation:
 *  \f[
 *    \mathbf{A}^{T}\mathbf{X}+\mathbf{X}\mathbf{A}-(\mathbf{X}\mathbf{B}+\mathbf{N})\mathbf{R}^{-1}(\mathbf{B}^{T}\mathbf{X}+\mathbf{N}^{T})+\mathbf{Q}=\mathbf{0}
 *  \f]
 *
 *
 * \note
 *  -# Matrices \f$\mathbf{Q}\f$ and \f$\mathbf{R}\f$ should be symmetric; if they
 *     are not, they will be replaced with \f$\mathbf{Q}'\f$ and \f$\mathbf{R}'\f$
 *     such that:
 *     \f{align*}
 *       \mathbf{Q}' &= \frac{\mathbf{Q}+\mathbf{Q}^T}{2}, \\
 *       \mathbf{R}' &= \frac{\mathbf{R}+\mathbf{R}^T}{2}.
 *     \f}
 *  -# Such a Discrete-time LQ regulation problem is equivalent to a Discrete-time
 *     LQ state-feedback regulation problem with the following weighting matrices:
 *     \f[
 *       \begin{pmatrix}\mathbf{\bar{Q}} & \mathbf{\bar{N}} \\ \mathbf{\bar{N}}^{T} & \mathbf{\bar{R}}\end{pmatrix} = \begin{pmatrix}\mathbf{C}^{T} & \mathbf{0} \\ \mathbf{D}^{T} & \mathbf{I}\end{pmatrix} \begin{pmatrix}\mathbf{Q} & \mathbf{N} \\ \mathbf{N}^{T} & \mathbf{R}\end{pmatrix} \begin{pmatrix}\mathbf{C} & \mathbf{D} \\ \mathbf{0} & \mathbf{I}\end{pmatrix}
 *     \f]
 *     that is:
 *     \f{align*}
 *       \mathbf{\bar{Q}} &= \mathbf{C}^{T}\mathbf{Q}\mathbf{C}, \\
 *       \mathbf{\bar{R}} &= \mathbf{R}+\mathbf{D}^{T}\mathbf{Q}\mathbf{D}+\mathbf{N}^{T}\mathbf{D}+\mathbf{N}\mathbf{D}^T, \\
 *       \mathbf{\bar{N}} &= \mathbf{C}^{T}\left(\mathbf{Q}\mathbf{D}+\mathbf{N}\right).
 *     \f}
 * .
 *
 * <hr/>
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

#ifndef DCS_CONTROL_DLQRY_HPP
#define DCS_CONTROL_DLQRY_HPP


#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/traits/layout_type.hpp>
#include <complex>
#include <dcs/assert.hpp>
#include <dcs/control/design/dlqr.hpp>
#include <dcs/debug.hpp>
#include <dcs/macro.hpp>
#include <stdexcept>


namespace dcs { namespace control {

namespace detail { namespace /*<unnamed>*/ {

/// Solve the Infinite-horizon, Discrete-time Linear Quadratic Regulator problem.
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
void dlqry2dlqr(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
				::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
				::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
				::boost::numeric::ublas::matrix_expression<QMatrixT>& Q,
				::boost::numeric::ublas::matrix_expression<RMatrixT>& R,
				::boost::numeric::ublas::matrix_expression<NMatrixT>& N)
{
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING(A);
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING(B);

	namespace ublasx = ::boost::numeric::ublasx;
	namespace ublas = ::boost::numeric::ublas;

	typedef typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
					typename ::boost::numeric::ublas::promote_traits<
						typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
						typename ::boost::numeric::ublas::promote_traits<
							typename ::boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
							typename ::boost::numeric::ublas::promote_traits<
								typename ::boost::numeric::ublas::matrix_traits<RMatrixT>::value_type,
								typename ::boost::numeric::ublas::matrix_traits<NMatrixT>::value_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type
		>::promote_type value_type;
	typedef typename ublasx::layout_type<AMatrixT>::type layout_type;
	typedef ublas::matrix<value_type,layout_type> work_matrix_type;


	// Derive parameters of equivalent LQR problem
	// ND = N'*D
	work_matrix_type ND(ublas::prod(ublas::trans(N), D));
	// QQ = C'*Q*C
	work_matrix_type QQ(ublas::prod(ublas::trans(C), Q));
	QQ = ublas::prod(QQ, C);
	// RR = R + D'*Q*D + ND + ND';
	work_matrix_type RR(ublas::prod(ublas::trans(D), Q)); // R == D'*Q
	RR = ublas::prod(RR, D); // RR == D'*Q*D
	RR += R + ND + ublas::trans(ND); // RR == D'*Q*D+R+ND+ND'
	// NN = C'*(Q*D + N);
	work_matrix_type NN(ublas::prod(Q, D) + N); // NN == Q*D+N
	NN = ublas::prod(ublas::trans(C), NN); // NN == C'*(Q*D+N)

	Q() = QQ;
	R() = RR;
	N() = NN;
}

}} // Namespace detail::<unnamed>


/**
 * \brief Infinite-horizon Discrete Linear Quadratic controller for
 *  output-feedback regulation.
 *
 * \tparam RealT The type for real numbers.
 *
 * \author Marco Guazzone, marco.guazzone@gmail.com
 */
template <typename RealT>
class dlqry_controller
{
	public: typedef RealT real_type;
	public: typedef ::std::complex<real_type> complex_type;
	public: typedef ::boost::numeric::ublas::matrix<real_type> matrix_type;
	public: typedef ::boost::numeric::ublas::vector<complex_type> vector_type;


	/// Default constructor
	public: dlqry_controller()
	{
		// empty
	}


	/// A constructor
	public: template <
				typename QMatrixT,
				typename RMatrixT,
				typename NMatrixT
		> dlqry_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						   ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
		: Q_(Q),
		  R_(R),
		  N_(N)
	{
		// Empty
	}


	/// A constructor
	public: template <
				typename QMatrixT,
				typename RMatrixT
		> dlqry_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
		: Q_(Q),
		  R_(R),
		  N_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(Q),::boost::numeric::ublasx::num_rows(R)))
	{
		// Empty
	}


	public: template<typename MatrixT>
		void Q(::boost::numeric::ublas::matrix_expression<MatrixT> const& X)
	{
		Q_ = X;
	}


	public: matrix_type Q() const
	{
		return Q_;
	}


	public: matrix_type& Q()
	{
		return Q_;
	}


	public: template<typename MatrixT>
		void R(::boost::numeric::ublas::matrix_expression<MatrixT> const& X)
	{
		R_ = X;
	}


	public: matrix_type R() const
	{
		return R_;
	}


	public: matrix_type& R()
	{
		return R_;
	}


	public: template<typename MatrixT>
		void N(::boost::numeric::ublas::matrix_expression<MatrixT> const& X)
	{
		N_ = X;
	}


	public: matrix_type N() const
	{
		return N_;
	}


	public: matrix_type& N()
	{
		return N_;
	}


	/**
	 * \brief Solve the Infinite-horizon Discrete-time Linear Quadratic Regulator
	 *  problem.
	 *
	 * \tparam AMatrixT The type of the state matrix for the state-space system
	 *  model.
	 * \tparam BMatrixT The type of the input matrix for the state-space
	 *  system model.
	 * \tparam CMatrixT The type of the output matrix for the state-space
	 *  system model.
	 * \tparam DMatrixT The type of the feedforward matrix for the state-space
	 *  system model.
	 *
	 *  \param A The state matrix for the state-space system model.
	 *  \param B The input matrix for the state-space system model.
	 *  \param C The output matrix for the state-space system model.
	 *  \param D The feedforward matrix for the state-space system model.
	 */
	public: template <typename AMatrixT, typename BMatrixT, typename CMatrixT, typename DMatrixT>
		void solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
		           ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C, ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D)
	{
		// pre: A must be square
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(A) == ::boost::numeric::ublasx::num_columns(A),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Matrix A must be square.")
			);
		// pre: num_rows(A) == num_rows(B)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(A) == ::boost::numeric::ublasx::num_rows(B),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices A and B are not conformant.")
			);
		// pre: num_rows(A) == num_columns(C)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(A) == ::boost::numeric::ublasx::num_columns(C),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices A and C are not conformant.")
			);
		// pre: num_columns(B) == num_columns(D)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_columns(B) == ::boost::numeric::ublasx::num_columns(D),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices B and D are not conformant.")
			);
		// pre: num_rows(C) == num_rows(D)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(C) == ::boost::numeric::ublasx::num_rows(D),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices C and D are not conformant.")
			);
		// pre: Q must be square
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(Q) == ::boost::numeric::ublasx::num_columns(Q),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Matrix Q must be square.")
			);
		// pre: num_rows(C) == num_rows(Q)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(C) == ::boost::numeric::ublasx::num_rows(Q),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices C and Q are not conformant.")
			);
		// pre: R must be square
		DCS_ASSERT(
				::boost::numeric::ublasx::num_rows(R) == ::boost::numeric::ublasx::num_columns(R),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Matrix R must be square.")
			);
		// pre: num_columns(D) == num_rows(R)
		DCS_ASSERT(
				::boost::numeric::ublasx::num_columns(D) == ::boost::numeric::ublasx::num_rows(R),
				throw ::std::invalid_argument("[dcs::control::dlqry_controller] Size of matrices D and R are not conformant.")
			);

		matrix_type QQ(Q_);
		matrix_type RR(R_);
		matrix_type NN(N_);

		detail::dlqry2dlqr(A, B, C, D, QQ, RR, NN);
		dlqr_controller<real_type> lq(QQ, RR, NN);
		lq.solve(A, B);
		K_ = lq.gain();
		S_ = lq.are_solution();
		e_ = lq.eigenvalues();
	}


	/// Return the optimal state-feedback matrix gain.
	public: matrix_type gain() const
	{
		return K_;
	}


	/// Return the associated DARE solution.
	public: matrix_type are_solution() const
	{
		return S_;
	}


	/// Return the closed-loop eigenvalues.
	public: vector_type eigenvalues() const
	{
		return e_;
	}


	public: template <typename VectorExprT>
		vector_type control(::boost::numeric::ublas::vector_expression<VectorExprT> const& x) const
	{
		// pre: size(x) == num_columns(K_)
		DCS_ASSERT(
			::boost::numeric::ublasx::size(x) == ::boost::numeric::ublasx::num_columns(K_),
			throw ::std::invalid_argument("[dcs::control::dlqry_controller::control] Wrong state dimensiion.")
		);

		return -::boost::numeric::ublas::prod(K_, x);
	}


	public: template <typename MatrixExprT>
		matrix_type control(::boost::numeric::ublas::matrix_expression<MatrixExprT> const& X) const
	{
		// pre: num_columns(X) == num_columns(K_)
		DCS_ASSERT(
			::boost::numeric::ublasx::num_columns(X) == ::boost::numeric::ublasx::num_columns(K_),
			throw ::std::invalid_argument("[dcs::control::dlqry_controller::control] Wrong state dimensiion.")
		);

		return -::boost::numeric::ublas::prod(K_, ::boost::numeric::ublas::trans(X)); 
	}


	/// The error weigthed matrix.
	private: matrix_type Q_;
	/// The control weigthed matrix.
	private: matrix_type R_;
	/// The cross-coupling weigthed matrix.
	private: matrix_type N_;
	/// The optimal feedback gain matrix.
	private: matrix_type K_;
	/// The solution to the associated DARE.
	private: matrix_type S_;
	/// The closed-loop eigenvalues.
	private: vector_type e_;
};


/**
 * \brief Infinite-horizon Discrete Linear Quadratic Controller design with
 *  weighting output.
 *
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 * \tparam NMatrixT The type of the cross-coupling weighting matrix.
 *
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \param N The cross-coupling weighting matrix.
 * \return An object representing an Infinite-horizon Discrete-time Linear
 *  Quadratic state-feedback Regulator controller.
 */
template <
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
dlqry_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename boost::numeric::ublas::matrix_traits<RMatrixT>::value_type,
			typename boost::numeric::ublas::matrix_traits<NMatrixT>::value_type
		>::promote_type
	>::promote_type
> dlqry(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
	    ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
	    ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
				typename ::boost::numeric::ublas::promote_traits<
					typename boost::numeric::ublas::matrix_traits<RMatrixT>::value_type,
					typename boost::numeric::ublas::matrix_traits<NMatrixT>::value_type
				>::promote_type
			>::promote_type value_type;

	return dlqry_controller<value_type>(Q, R, N);
}


/**
 * \brief Infinite-horizon Discrete Linear Quadratic Controller design with
 *  weighting output.
 *
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 *
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \return An object representing an Infinite-horizon Discrete-time Linear
 *  Quadratic state-feedback Regulator controller.
 */
template <
	typename QMatrixT,
	typename RMatrixT
>
inline
dlqry_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
		typename boost::numeric::ublas::matrix_traits<RMatrixT>::value_type
	>::promote_type
> dlqry(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
	    ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
				typename boost::numeric::ublas::matrix_traits<RMatrixT>::value_type
			>::promote_type value_type;

	return dlqry_controller<value_type>(Q, R);
}


/**
 * \brief Optimal state-feedback control by Infinite-horizon Discrete Linear
 *  Quadratic Controller design with weighting output.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam CMatrixT The type of the output matrix for the controlled state-space
 *  model.
 * \tparam DMatrixT The type of the feeforward matrix for the controlled
 *  state-space model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 * \tparam NMatrixT The type of the cross-coupling weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param C The output matrix for the controlled state-space model.
 * \param D The feedforward matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \param N The cross-coupling weighting matrix.
 * \return The optimal state-feedback gain matrix \f$\mathbf{K}\f$.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
//::boost::numeric::ublas::matrix<double> dlqry_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
					typename ::boost::numeric::ublas::promote_traits<
						typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
						typename ::boost::numeric::ublas::promote_traits<
							typename ::boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
							typename ::boost::numeric::ublas::promote_traits<
								typename ::boost::numeric::ublas::matrix_traits<RMatrixT>::value_type,
								typename ::boost::numeric::ublas::matrix_traits<NMatrixT>::value_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type
		>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqry_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			  ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
			  ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
			  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			  ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
					typename ::boost::numeric::ublas::promote_traits<
						typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
						typename ::boost::numeric::ublas::promote_traits<
							typename ::boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
							typename ::boost::numeric::ublas::promote_traits<
								typename ::boost::numeric::ublas::matrix_traits<RMatrixT>::value_type,
								typename ::boost::numeric::ublas::matrix_traits<NMatrixT>::value_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type
		>::promote_type value_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type QQ(Q);
	out_matrix_type RR(R);
	out_matrix_type NN(N);

	detail::dlqry2dlqr(A, B, C, D, QQ, RR, NN);

	return dlqr_solve(A, B, QQ, RR, NN);
}


/**
 * \brief Optimal state-feedback control by Infinite-horizon Discrete Linear
 *  Quadratic Controller design with weighting output.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam CMatrixT The type of the output matrix for the controlled state-space
 *  model.
 * \tparam DMatrixT The type of the feedforward matrix for the controlled
 *  state-space  model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param C The output matrix for the controlled state-space model.
 * \param D The feedforward matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \return The optimal state-feedback gain matrix \f$\mathbf{K}\f$.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename QMatrixT,
	typename RMatrixT
>
inline
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
					typename ::boost::numeric::ublas::promote_traits<
						typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
						typename ::boost::numeric::ublas::promote_traits<
							typename ::boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
							typename ::boost::numeric::ublas::matrix_traits<RMatrixT>::value_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type
		>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqry_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			  ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
			  ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
			  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
					typename ::boost::numeric::ublas::promote_traits<
						typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
						typename ::boost::numeric::ublas::promote_traits<
							typename ::boost::numeric::ublas::matrix_traits<QMatrixT>::value_type,
							typename ::boost::numeric::ublas::matrix_traits<RMatrixT>::value_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type
		>::promote_type value_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type QQ(Q);
	out_matrix_type RR(R);
	out_matrix_type NN(
		::boost::numeric::ublas::zero_matrix<value_type>(
			::boost::numeric::ublasx::num_rows(Q),
			::boost::numeric::ublasx::num_rows(R)
		)
	);

	detail::dlqry2dlqr(A, B, C, D, QQ, RR, NN);

	return dlqr_solve(A, B, QQ, RR, NN);
}


}} // Namespace dcs::control


#endif // DCS_CONTROL_DLQRY_HPP
