/**
 * \file dcs/control/design/dlqr.hpp
 *
 * \brief Infinite-horizon Discrete Linear Quadratic state-feedback Regulator.
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
 *   J(\mathbf{u}) &= \sum_{k=1}^{\infty}{\mathbf{x}^T(k)\mathbf{Q}\mathbf{x}(k) + \mathbf{u}^T(k)\mathbf{R}\mathbf{u}(k) + 2\mathbf{x}^T(k)\mathbf{N}\mathbf{u}(k)} \\
 *                 &= \sum_{k=1}^{\infty}{\begin{pmatrix}\mathbf{x}^T(k) & \mathbf{u}^T(k)\end{pmatrix}\begin{pmatrix}\mathbf{Q} & \mathbf{N} \\ \mathbf{N}^T & \mathbf{R}\end{pmatrix}\begin{pmatrix}\mathbf{x}(k) \\ \mathbf{u}(k)\end{pmatrix}}
 * \f}
 * where the <em>error weighting matrix</em>
 *           \f$\mathbf{Q}\in\mathbb{R}^{n \times n}\f$ is a positive
 *           semi-definite real matrix,
 *       the <em>control weighting matrix</em>
 *           \f$\mathbf{R}\in\mathbb{R}^{m \times m}\f$ is a positive
 *           definite real matrix, and
 *       the <em>cross-coupling weighting matrix</em>
 *           \f$\mathbf{N}\in\mathbb{R}^{n \times m}\f$ is a real matrix.
 * Matrices \f$\mathbf{Q}\f$ and \f$\mathbf{R}\f$ should be symmetric; if they
 * are not, they are replaced with \f$\mathbf{Q}'\f$ and \f$\mathbf{R}'\f$ such
 * that:
 * \f{align*}
 *   \mathbf{Q}' &= \frac{\mathbf{Q}+\mathbf{Q}^T}{2}, \\
 *   \mathbf{R}' &= \frac{\mathbf{R}+\mathbf{R}^T}{2}.
 * \f}
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

#ifndef DCS_CONTROL_DLQR_HPP
#define DCS_CONTROL_DLQR_HPP


#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/eigen.hpp>
#include <boost/numeric/ublasx/operation/max.hpp>
#include <boost/numeric/ublasx/operation/min.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/traits/layout_type.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <complex>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <dcs/control/solver/dare.hpp>
#include <stdexcept>


namespace dcs { namespace control {

namespace detail { namespace /*<unnamed>*/ {

/// Solve the Infinite-horizon, Discrete-time Linear Quadratic Regulator problem.
template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT,
	typename KMatrixT,
	typename SMatrixT,
	typename EVectorT
>
void dlqr(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
		  boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
		  boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
		  boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
		  boost::numeric::ublas::matrix_expression<NMatrixT> const& N,
		  bool want_S,
		  bool want_e,
		  boost::numeric::ublas::matrix_container<KMatrixT>& K,
		  boost::numeric::ublas::matrix_container<SMatrixT>& S,
		  boost::numeric::ublas::vector_container<EVectorT>& e)
{
	namespace ublasx = ::boost::numeric::ublasx;
	namespace ublas = ::boost::numeric::ublas;

	typedef AMatrixT A_matrix_type;
	typedef BMatrixT B_matrix_type;
	typedef QMatrixT Q_matrix_type;
	typedef RMatrixT R_matrix_type;
	typedef NMatrixT N_matrix_type;
//	typedef typename ublas::promote_traits<
//				typename ublas::matrix_traits<A_matrix_type>::value_type,
//				typename ublas::promote_traits<
//					typename ublas::matrix_traits<B_matrix_type>::value_type,
//					typename ublas::promote_traits<
//						typename ublas::matrix_traits<Q_matrix_type>::value_type,
//						typename ublas::promote_traits<
//							typename ublas::matrix_traits<R_matrix_type>::value_type,
//							typename ublas::matrix_traits<N_matrix_type>::value_type
//						>::promote_type
//					>::promote_type
//				>::promote_type
//			>::promote_type value_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<A_matrix_type>::value_type,
				typename ublas::matrix_traits<B_matrix_type>::value_type
			>::promote_type value_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<A_matrix_type>::size_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<B_matrix_type>::size_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<Q_matrix_type>::size_type,
						typename ublas::promote_traits<
							typename ublas::matrix_traits<R_matrix_type>::size_type,
							typename ublas::matrix_traits<N_matrix_type>::size_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type size_type;
	typedef typename ublas::type_traits<value_type>::real_type real_type;
	typedef ::std::complex<real_type> complex_type;

	// precondition: Q must be real
	BOOST_STATIC_ASSERT(
		::boost::is_floating_point<typename ublas::matrix_traits<Q_matrix_type>::value_type>::value
		//"dcs_control_lqr_Q_real"
	);
	// precondition: R must be real
	BOOST_STATIC_ASSERT(
		::boost::is_floating_point<typename ublas::matrix_traits<R_matrix_type>::value_type>::value
		//"dcs_control_lqr_R_real"
	);
	// precondition: N must be real
	BOOST_STATIC_ASSERT(
		::boost::is_floating_point<typename ublas::matrix_traits<N_matrix_type>::value_type>::value
		//"dcs_control_lqr_N_real"
	);

	size_type A_nr = ublasx::num_rows(A);
	size_type A_nc = ublasx::num_columns(A);
	size_type B_nr = ublasx::num_rows(B);
	size_type B_nc = ublasx::num_columns(B);
	size_type Q_nr = ublasx::num_rows(Q);
	size_type Q_nc = ublasx::num_columns(Q);
	size_type R_nr = ublasx::num_rows(Q);
	size_type R_nc = ublasx::num_columns(Q);
	size_type N_nr = ublasx::num_rows(N);
	size_type N_nc = ublasx::num_columns(N);

	// precondition: A is square
	DCS_ASSERT(
		A_nr == A_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] State matrix A must be a square matrix.")
	);
	// precondition: num_rows(B) == num_rows(A)
	DCS_ASSERT(
		B_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the input matrix B must be the same of the state matrix A.")
	);
	// precondition: Q is square
	DCS_ASSERT(
		Q_nr == Q_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] Error weighting matrix Q must be a square matrix.")
	);
	// precondition: num_rows(Q) == num_rows(A)
	DCS_ASSERT(
		Q_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the error weighting matrix Q must be the same of the state matrix A.")
	);
	// precondition: R is square
	DCS_ASSERT(
		R_nr == R_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] Control weighting matrix Q must be a square matrix.")
	);
	// precondition: num_rows(R) == num_rows(B)
	DCS_ASSERT(
		R_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the control weighting matrix Q must be the same of the input matrix B.")
	);
	// precondition: num_rows(N) == num_rows(A) && num_columns(N) == num_columns(B)
	DCS_ASSERT(
		N_nr == A_nr && N_nc == B_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] The cross-term weighting matrix N must have the same number of rows of the state matrix A and be the same number of columns of the input matrix B.")
	);


	const value_type eps = ::std::numeric_limits<real_type>::epsilon();

	Q_matrix_type tmp_Q;
	R_matrix_type tmp_R;

	//// Check/enforce symmetry and check positivity
	//if (ublas::norm_1(ublas::trans(Q)-Q) > (100*eps*ublas::norm_1(Q)))
	//{
	//	DCS_DEBUG_TRACE("The error weighting matrix Q is not symmetric and will be replaced by (Q+Q')/2.");
	//	tmp_Q = (Q+ublas::trans(Q))/static_cast<value_type>(2);
	//}
	//else
	//{
	//	tmp_Q = Q;
	//}
	//if (ublas::norm_1(ublas::trans(R)-R) > (100*eps*ublas::norm_1(R)))
	//{
	//	DCS_DEBUG_TRACE("The control weighting matrix R is not symmetric and has been replaced by (R+R')/2.");
	//	tmp_R = (R+ublas::trans(R))/static_cast<value_type>(2);
	//}
	//else
	//{
	//	tmp_R = R;
	//}

	// Enforce symmetry
	tmp_Q = (Q+ublas::trans(Q))/static_cast<real_type>(2);
	tmp_R = (R+ublas::trans(R))/static_cast<real_type>(2);


	ublas::vector<complex_type> v;

	ublas::vector<real_type> v_R;
	ublasx::eigenvalues(tmp_R, v);
	v_R = ublas::real(v);

	//FIXME: why using real(...)?
	//ublas::vector< std::complex<value_type> > v_QNR;
	ublas::vector<real_type> v_QNR;
	ublasx::eigenvalues(
		ublasx::cat_columns(
			ublasx::cat_rows(tmp_Q, N), 
			ublasx::cat_rows(ublas::trans(N), tmp_R) 
		),
		v
	);
	v_QNR = ublas::real(v);

	if (ublasx::min(v_R) <= real_type/*zero*/())
	{
		throw ::std::runtime_error("[dcs::control::detail::dlqr] The control weighting matrix R is not positive definite.");
	}
	else if (ublasx::min(v_QNR) < (-1.0e+2*eps*::std::max(real_type/*zero*/(), ublasx::max(v_QNR))))
	{
		DCS_DEBUG_TRACE("[dcs::control::detail::dlqr] The matrix [Q N;N' R] is not positive semi-definite.");
	}

	ublas::matrix<value_type> E = ublas::identity_matrix<value_type>(A_nr);

	dare_solver<value_type> solver;
	solver.solve(A, B, Q, R, N, E);

	K() = solver.gain();
	if (want_S)
	{
		S() = solver.solution();
	}
	if (want_e)
	{
		e() = solver.eigenvalues();
	}
}

}} // Namespace detail::<unnamed>


/**
 * \brief Infinite-horizon Discrete Linear Quadratic controller for
 *  state-feedback regulation.
 *
 * \tparam RealT The type for real numbers.
 *
 * \author Marco Guazzone, marco.guazzone@gmail.com
 */
template <typename RealT>
class dlqr_controller
{
	public: typedef RealT real_type;
	public: typedef ::std::complex<real_type> complex_type;
	public: typedef ::boost::numeric::ublas::matrix<real_type> matrix_type;
	public: typedef ::boost::numeric::ublas::vector<complex_type> vector_type;


	/// Default constructor
	public: dlqr_controller()
	{
		// empty
	}


	/// A constructor
	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename QMatrixT,
				typename RMatrixT,
				typename NMatrixT
		> dlqr_controller(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						  ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
		: Q_(Q),
		  R_(R),
		  N_(N)
	{
		solve(A, B);
	}


	/// A constructor
	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename QMatrixT,
				typename RMatrixT
		> dlqr_controller(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
		: Q_(Q),
		  R_(R),
		  N_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(B),::boost::numeric::ublasx::num_columns(B)))
	{
		solve(A, B);
	}


	/// A constructor
	public: template <
				typename QMatrixT,
				typename RMatrixT,
				typename NMatrixT
		> dlqr_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						  ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
		: Q_(Q),
		  R_(R),
		  N_(N)
	{
		// Empty
	}


	/**
	 * \brief Solve the Infinite-horizon Discrete-time Linear Quadratic Regulator
	 *  problem.
	 *
	 * \tparam AMatrixT The type of the state matrix for the state-space system
	 *  model.
	 * \tparam BMatrixT The type of the input matrix for the state-space
	 *  system model.
	 *
	 *  \param A The state matrix for the state-space system model.
	 *  \param B The input matrix for the state-space system model.
	 */
	public: template <typename AMatrixT, typename BMatrixT>
		void solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
	{
		detail::dlqr(A, B, Q_, R_, N_, true, true, K_, S_, e_);
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
		// preconditions: size(x) == num_columns(K_)
		DCS_ASSERT(
			::boost::numeric::ublasx::size(x) == ::boost::numeric::ublasx::num_columns(K_),
			throw ::std::invalid_argument("[dcs::control::dlqr_controller::control] Wrong state dimensiion.")
		);

		return -::boost::numeric::ublas::prod(K_, x);
	}


	public: template <typename MatrixExprT>
		matrix_type control(::boost::numeric::ublas::matrix_expression<MatrixExprT> const& X) const
	{
		// preconditions: num_columns(X) == num_columns(K_)
		DCS_ASSERT(
			::boost::numeric::ublasx::num_columns(X) == ::boost::numeric::ublasx::num_columns(K_),
			throw ::std::invalid_argument("[dcs::control::dlqr_controller::control] Wrong state dimensiion.")
		);

		return -::boost::numeric::ublas::prod(K_, ::boost::numeric::ublas::trans(X)); 
	}


	/// The error weigthed matrix.
	private: matrix_type Q_;
	/// The control weigthed matrix.
	private: matrix_type R_;
	/// The cross-coupling weigthed matrix.
	private: matrix_type N_;
	/// The optimal gain matrix.
	private: matrix_type K_;
	/// The solution to the associated DARE.
	private: matrix_type S_;
	/// The closed-loop eigenvalues which gives the closed-loop poles of \f$A-BK\f$.
	private: vector_type e_;
};


/**
 * \brief Infinite-horizon Discrete Linear Quadratic Controller design.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 * \tparam NMatrixT The type of the cross-coupling weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \param N The cross-coupling weighting matrix.
 * \return An object representing an Infinite-horizon Discrete-time Linear
 *  Quadratic state-feedback Regulator controller.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
dlqr_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
	>::promote_type
> dlqr(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
	   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
	   ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
	   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
	   ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
			>::promote_type value_type;

	return dlqr_controller<value_type>(A, B, Q, R, N);
}


/**
 * \brief Infinite-horizon Discrete Linear Quadratic Controller design.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \return An object representing an Infinite-horizon Discrete-time Linear
 *  Quadratic state-feedback Regulator controller.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT
>
inline
dlqr_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
	>::promote_type
> dlqr(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
	   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
	   ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
	   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
			>::promote_type value_type;

	return dlqr_controller<value_type>(A, B, Q, R);
}


/**
 * \brief Optimal state-feedback control by Infinite-horizon Discrete Linear
 *  Quadratic Controller design.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 * \tparam NMatrixT The type of the cross-coupling weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \param N The cross-coupling weighting matrix.
 * \return The optimal state-feedback gain matrix \f$\mathbf{K}\f$.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
		>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			  ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
			>::promote_type value_type;
	typedef ::std::complex<typename ::boost::numeric::ublas::type_traits<value_type>::real_type> complex_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type K;
	out_matrix_type dummy_S;
	::boost::numeric::ublas::vector<complex_type> dummy_e;

	detail::dlqr(
		A,
		B,
		Q,
		R,
		N,
		false,
		false,
		K,
		dummy_S,
		dummy_e
	);

	return K;
}


/**
 * \brief Optimal state-feedback control by Infinite-horizon Discrete Linear
 *  Quadratic Controller design.
 *
 * \tparam AMatrixT The type of the state matrix for the controlled state-space
 *  model.
 * \tparam BMatrixT The type of the input matrix for the controlled state-space
 *  model.
 * \tparam QMatrixT The type of the error weighting matrix.
 * \tparam RMatrixT The type of the control weighting matrix.
 *
 * \param A The state matrix for the controlled state-space model.
 * \param B The input matrix for the controlled state-space model.
 * \param Q The error weighting matrix.
 * \param R The control weighting matrix.
 * \return The optimal state-feedback gain matrix \f$\mathbf{K}\f$.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename KMatrixT
>
inline
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
		>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqr_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			 ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			 ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
{
	typedef typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
				typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type
			>::promote_type value_type;
	typedef ::std::complex<typename ::boost::numeric::ublas::type_traits<value_type>::real_type> complex_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type K;
	out_matrix_type dummy_S;
	::boost::numeric::ublas::vector<complex_type> dummy_e;

	detail::dlqr(
		A,
		B,
		Q,
		R,
		::boost::numeric::ublas::zero_matrix<value_type>(
			::boost::numeric::ublasx::num_rows(B),
			::boost::numeric::ublasx::num_columns(B)
		),
		false,
		false,
		K,
		dummy_S,
		dummy_e
	);

	return K;
}


}} // Namespace dcs::control


#endif // DCS_CONTROL_DLQR_HPP
