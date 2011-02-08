/**
 * \file dcs/control/design/dlqi.hpp
 *
 * \brief Infinite-horizon Discrete Linear Quadratic Integral state-feedback.
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
 * the <em>Infinite-horizon, Discrete Linear-Quadratic (LQ) Integral
 * state-feedback</em> controller problem calculates the optimal \e gain matrix
 * \f$\mathbf{K}\f$ such that the augmented state-feedback law:
 * \f{align}
 *   \mathbf{u}(k) &= -K\mathbf{z}(k) \\
 *                 &= -K(\mathbf{x}(k) \mathbf{x}_i(k))^T
 * \f}
 * minimizes the following quadratic cost function (<em>performance index</em>):
 * \f{align*}
 *   J(\mathbf{u}) &= \sum_{k=1}^{\infty}{\mathbf{z}^T(k)\mathbf{Q}\mathbf{z}(k) + \mathbf{u}^T(k)\mathbf{R}\mathbf{u}(k) + 2\mathbf{z}^T(k)\mathbf{N}\mathbf{u}(k)} \\
 *                 &= \sum_{k=1}^{\infty}{\begin{pmatrix}\mathbf{z}^T(k) & \mathbf{u}^T(k)\end{pmatrix}\begin{pmatrix}\mathbf{Q} & \mathbf{N} \\ \mathbf{N}^T & \mathbf{R}\end{pmatrix}\begin{pmatrix}\mathbf{z}(k) \\ \mathbf{u}(k)\end{pmatrix}}
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
 * Copyright (C) 2009-2011  Distributed Computing System (DCS) Group, Computer
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

#ifndef DCS_CONTROL_DLQI_HPP
#define DCS_CONTROL_DLQI_HPP


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
	typename CMatrixT,
	typename DMatrixT,
	typename RealT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT,
	typename KMatrixT,
	typename SMatrixT,
	typename EVectorT
>
void dlqi(boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
		  boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
		  boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
		  boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
		  RealT ts,
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
	typedef CMatrixT C_matrix_type;
	typedef DMatrixT D_matrix_type;
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
				typename ublas::promote_traits<
					typename ublas::matrix_traits<B_matrix_type>::value_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<C_matrix_type>::value_type,
						typename ublas::matrix_traits<D_matrix_type>::value_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<A_matrix_type>::size_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<B_matrix_type>::size_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<C_matrix_type>::size_type,
						typename ublas::promote_traits<
							typename ublas::matrix_traits<D_matrix_type>::size_type,
							typename ublas::promote_traits<
								typename ublas::matrix_traits<Q_matrix_type>::size_type,
								typename ublas::promote_traits<
									typename ublas::matrix_traits<R_matrix_type>::size_type,
									typename ublas::matrix_traits<N_matrix_type>::size_type
								>::promote_type
							>::promote_type
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
	size_type C_nr = ublasx::num_rows(C);
//	size_type C_nc = ublasx::num_columns(C);
//	size_type D_nr = ublasx::num_rows(D);
//	size_type D_nc = ublasx::num_columns(D);
	size_type Q_nr = ublasx::num_rows(Q);
	size_type Q_nc = ublasx::num_columns(Q);
	size_type R_nr = ublasx::num_rows(R);
	size_type R_nc = ublasx::num_columns(R);
	size_type N_nr = ublasx::num_rows(N);
	size_type N_nc = ublasx::num_columns(N);

	// precondition: A is square
	DCS_ASSERT(
		A_nr == A_nc,
		throw ::std::invalid_argument("[dcs::control::dlqi] State matrix A must be a square matrix.")
	);
	// precondition: num_rows(B) == num_rows(A)
	DCS_ASSERT(
		B_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqi] The number of rows of the input matrix B must be the same of the state matrix A.")
	);
	// precondition: Q is square
	DCS_ASSERT(
		Q_nr == Q_nc,
		throw ::std::invalid_argument("[dcs::control::dlqi] Error weighting matrix Q must be a square matrix.")
	);
	// precondition: num_rows(Q) == num_rows(A)+num_rows(C)
	DCS_ASSERT(
		Q_nr == (A_nr+C_nr),
		throw ::std::invalid_argument("[dcs::control::dlqi] The number of rows of the error weighting matrix Q must be the equal to the sum between the number of rows of the state matrix A and the number of outputs.")
	);
	// precondition: R is square
	DCS_ASSERT(
		R_nr == R_nc,
		throw ::std::invalid_argument("[dcs::control::dlqi] Control weighting matrix Q must be a square matrix.")
	);
	// precondition: num_rows(R) == num_rows(B)
	DCS_ASSERT(
		R_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqi] The number of rows of the control weighting matrix Q must be the same of the input matrix B.")
	);
	// precondition: num_rows(N) == (num_rows(A)+num_rows(C)) && num_columns(N) == num_columns(B)
	DCS_ASSERT(
		N_nr == (A_nr+C_nr) && N_nc == B_nc,
		throw ::std::invalid_argument("[dcs::control::dlqi] The cross-term weighting matrix N must have a number of rows equal to the sum between the number of rows of the state matrix A and the number of outputs, and a number of columns equal to the number of columns of the input matrix B.")
	);


	const value_type eps = ::std::numeric_limits<real_type>::epsilon();

	Q_matrix_type tmp_Q;
	R_matrix_type tmp_R;

	// Enforce symmetry

	tmp_Q = (Q+ublas::trans(Q))/static_cast<real_type>(2);
	tmp_R = (R+ublas::trans(R))/static_cast<real_type>(2);

	// Check positivity

	ublas::vector<complex_type> v_R;
	ublasx::eigenvalues(tmp_R, v_R);

	ublas::vector<complex_type> v_QNR;
	ublasx::eigenvalues(
		// [tmp_Q N; N' tmp_R]
		ublasx::cat_columns(
			ublasx::cat_rows(tmp_Q, N), 
			ublasx::cat_rows(ublas::trans(N), tmp_R) 
		),
		v_QNR
	);

	if (::std::real(ublasx::min(v_R)) <= real_type/*zero*/())
	{
		throw ::std::runtime_error("[dcs::control::detail::dlqi] The control weighting matrix R is not positive definite.");
	}
	else if (::std::real(ublasx::min(v_QNR)) < (-1.0e+2*eps*::std::max(real_type/*zero*/(), ::std::real(ublasx::max(v_QNR)))))
	{
		DCS_DEBUG_TRACE("[dcs::control::detail::dlqi] The matrix [Q N;N' R] is not positive semi-definite.");
		::std::clog << "[Warning] The matrix [Q N;N' R] is not positive semi-definite." << ::std::endl;
	}

	// Form the augmented system with an integrator
	//  xx(k+1) = AAxx(k)+BBu(k)
	// where
	//  AA = [A 0; -C*abs(t_s) I]
	//  BB = [B; -D*abs(t_s)]
	//  xx(k) = [x(k); xi(k)]
	//  t_s = sampling time
	//  xi = integrator output
	// NOTE: integration is based on Forward Euler formula
	//   xi(k+1) = xi(k)+t_s*(r(k)-y(k))
	ts = ::std::abs(ts);
	ublas::matrix<value_type> AA;
	AA = ublasx::cat_columns(
			ublasx::cat_rows(A, ublas::zero_matrix<value_type>(A_nr,C_nr)),
			ublasx::cat_rows(-C*ts, ublas::identity_matrix<value_type>(C_nr,C_nr))
		);
	ublas::matrix<value_type> BB;
	BB = ublasx::cat_columns(B, -D*ts);

	ublas::matrix<value_type> E = ublas::identity_matrix<value_type>(A_nr+C_nr);

	dare_solver<value_type> solver;
	solver.solve(AA, BB, Q, R, N, E);

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
 * \brief Infinite-horizon Discrete Linear Quadratic Integral controller for
 *  state-feedback regulation.
 *
 * \tparam RealT The type for real numbers.
 *
 * \note
 *  Inspired by the \c dlqi MATLAB function.
 *
 * \author Marco Guazzone, marco.guazzone@mfn.unipmn.it
 */
template <typename RealT>
class dlqi_controller
{
	public: typedef RealT real_type;
	public: typedef ::std::complex<real_type> complex_type;
	public: typedef ::boost::numeric::ublas::matrix<real_type> matrix_type;
	public: typedef ::boost::numeric::ublas::vector<complex_type> vector_type;


	/// Default constructor
	public: dlqi_controller()
	{
		// empty
	}


	/// A constructor
	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename CMatrixT,
				typename DMatrixT,
				typename QMatrixT,
				typename RMatrixT,
				typename NMatrixT
		> dlqi_controller(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						  ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
						  ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
						  real_type ts,
						  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						  ::boost::numeric::ublas::matrix_expression<NMatrixT> const& N)
		: Q_(Q),
		  R_(R),
		  N_(N)
	{
		solve(A, B, C, D, ts);
	}


	/// A constructor
	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename CMatrixT,
				typename DMatrixT,
				typename QMatrixT,
				typename RMatrixT
		> dlqi_controller(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						  ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						  ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
						  ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
						  real_type ts,
						  ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
		: Q_(Q),
		  R_(R),
		  N_(::boost::numeric::ublas::zero_matrix<real_type>(
					::boost::numeric::ublasx::num_rows(Q),
					::boost::numeric::ublasx::num_rows(R)
				)
			)
	{
		solve(A, B, C, D, ts);
	}


	/// A constructor
	public: template <
				typename QMatrixT,
				typename RMatrixT,
				typename NMatrixT
		> dlqi_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
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
		> dlqi_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						  ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R)
		: Q_(Q),
		  R_(R),
		  N_(::boost::numeric::ublas::zero_matrix<real_type>(
					::boost::numeric::ublasx::num_rows(Q),
					::boost::numeric::ublasx::num_rows(R)
				)
			)
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
	public: template <typename AMatrixT,
					  typename BMatrixT,
					  typename CMatrixT,
					  typename DMatrixT>
		void solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
				   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
				   ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
					real_type ts)
	{
		detail::dlqi(A, B, C, D, ts, Q_, R_, N_, true, true, K_, S_, e_);
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
			throw ::std::invalid_argument("[dcs::control::dlqi_controller::control] Wrong state dimensiion.")
		);

		return -::boost::numeric::ublas::prod(K_, x);
	}


	public: template <typename MatrixExprT>
		matrix_type control(::boost::numeric::ublas::matrix_expression<MatrixExprT> const& X) const
	{
		// preconditions: num_columns(X) == num_columns(K_)
		DCS_ASSERT(
			::boost::numeric::ublasx::num_columns(X) == ::boost::numeric::ublasx::num_columns(K_),
			throw ::std::invalid_argument("[dcs::control::dlqi_controller::control] Wrong state dimensiion.")
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
 *
 * \note
 *  Inspired by the \c dlqi MATLAB function.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename RealT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
dlqi_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
				typename ::boost::numeric::ublas::promote_traits<
					typename boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
					RealT
				>::promote_type
			>::promote_type
		>::promote_type
	>::promote_type
> dlqi(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
	   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
	   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
	   ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
	   RealT ts,
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
							RealT
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;

	return dlqi_controller<value_type>(A, B, C, D, ts, Q, R, N);
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
 *
 * \note
 *  Inspired by the \c dlqi MATLAB function.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename RealT,
	typename QMatrixT,
	typename RMatrixT
>
inline
dlqi_controller<
	typename ::boost::numeric::ublas::promote_traits<
		typename boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
				typename ::boost::numeric::ublas::promote_traits<
					typename boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
					RealT
				>::promote_type
			>::promote_type
		>::promote_type
	>::promote_type
> dlqi(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
	   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
	   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
	   ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
	   RealT ts,
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
							RealT
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;

	return dlqi_controller<value_type>(A, B, C, D, ts, Q, R);
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
 *
 * \note
 *  Inspired by the \c dlqi MATLAB function.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename RealT,
	typename QMatrixT,
	typename RMatrixT,
	typename NMatrixT
>
inline
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqi_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
				typename ::boost::numeric::ublas::promote_traits<
					typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
					RealT
				>::promote_type
			>::promote_type
		>::promote_type
	>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqi_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			 ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
			 ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
			 RealT ts,
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
							RealT
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;
	typedef ::std::complex<typename ::boost::numeric::ublas::type_traits<value_type>::real_type> complex_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type K;
	out_matrix_type dummy_S;
	::boost::numeric::ublas::vector<complex_type> dummy_e;

	detail::dlqi(
		A,
		B,
		C,
		D,
		ts,
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
 *
 * \note
 *  Inspired by the \c dlqi MATLAB function.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT,
	typename RealT,
	typename QMatrixT,
	typename RMatrixT,
	typename KMatrixT
>
inline
//typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type dlqi_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
::boost::numeric::ublas::matrix<
	typename ::boost::numeric::ublas::promote_traits<
		typename ::boost::numeric::ublas::matrix_traits<AMatrixT>::value_type,
		typename ::boost::numeric::ublas::promote_traits<
			typename ::boost::numeric::ublas::matrix_traits<BMatrixT>::value_type,
			typename ::boost::numeric::ublas::promote_traits<
				typename ::boost::numeric::ublas::matrix_traits<CMatrixT>::value_type,
				typename ::boost::numeric::ublas::promote_traits<
					typename ::boost::numeric::ublas::matrix_traits<DMatrixT>::value_type,
					RealT
				>::promote_type
			>::promote_type
		>::promote_type
	>::promote_type,
	typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type
> dlqi_solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			 ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
			 ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D,
			 RealT ts,
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
							RealT
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;
	typedef ::std::complex<typename ::boost::numeric::ublas::type_traits<value_type>::real_type> complex_type;
	//typedef typename boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type out_matrix_type; //this doesn't work if A and B have different value types
	typedef ::boost::numeric::ublas::matrix<value_type, typename ::boost::numeric::ublasx::layout_type<AMatrixT>::type> out_matrix_type;

	out_matrix_type K;
	out_matrix_type dummy_S;
	::boost::numeric::ublas::vector<complex_type> dummy_e;

	detail::dlqi(
		A,
		B,
		C,
		D,
		ts,
		Q,
		R,
		::boost::numeric::ublas::zero_matrix<value_type>(
			::boost::numeric::ublasx::num_rows(Q),
			::boost::numeric::ublasx::num_columns(R)
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


#endif // DCS_CONTROL_DLQI_HPP
