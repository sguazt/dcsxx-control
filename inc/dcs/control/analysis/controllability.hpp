/**
 * \file dcs/control/analysis/controllability.hpp
 *
 * \brief Controllability for a state-space system.
 *
 * A system with internal state vector \f$\mathbf{x}\f$ is called
 * \e controllable if and only if the system states can be changed by changing
 * the system input.
 * Equivalently, a system is said to be controllable if for any initial state
 * \f$x(k_0)\f$ there exists a control sequence \f$u(k)\f$, with
 * \f$k=k_0,k_0+1,\ldots,k_f-1\f$, such that an arbitrary final state
 * \f$x(k_f)\f$ can be reached in finite time \f$k_f\f$.
 *
 * More specifically, given a continuous-time state-space system model:
 * \f[
 *   \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t)+\mathbf{B}\mathbf{u}(t)
 *   \mathbf{y}(t) = \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f]
 * or a discrete-time state-space system model:
 * \f[
 *   \mathbf{x}(t) = \mathbf{A}\mathbf{x}(t-1)+\mathbf{B}\mathbf{u}(t-1)
 *   \mathbf{y}(t) = \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f]
 * the following definitions apply [1]:
 * - The system is called <em>state controllable</em> at time \f$t_0\f$ if
 *   it there exists an unconstrained input signal \f$\mathbf{u}(t)\f$ that
 *   willl transfer an initial state to any final state in a finite time
 *   interval \f$t_0 \le t \le t_1\f$.
 * - The system is called <em>complete state controllable</em> (or simply
 *   \e controllable) at time \f$t_0\f$ if every state \f$x_0\f$ in the
 *   state-space model is controllable.
 * - The system is called <em>complete output controllable</em> if there exists
 *   an unconstrained input signal \f$\mathbf{u}(t)\f$ that will transfer any
 *   given initial output \f$\mathbf{y}(t_0)\f$ to any final output
 *   \f$\mathbf{y}(t_1)\f$ in a finite time interval \f$t_0 \le t \le t_1\f$.
 * .
 *
 * Formally, system controllability is defined in terms of the
 * <em>controllability matrix</em>.
 *
 * Given a continuous-time state-space system model:
 * \f{align*}
 *   \dot{\mathbf{x}}(t) &= \mathbf{A}\mathbf{x}(t)+\mathbf{B}\mathbf{u}(t),\\
 *   \mathbf{y}(t) &= \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f}
 * or a discrete-time state-space system model:
 * \f{align*}
 *   \mathbf{x}(t) &= \mathbf{A}\mathbf{x}(t-1)+\mathbf{B}\mathbf{u}(t-1),\\
 *   \mathbf{y}(t) &= \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f}
 * where:
 * - \f$\mathbf{A}\f$ is a \f$n \times n}\f$ matrix,
 * - \f$\mathbf{B}\f$ is a \f$n \times m}\f$ matrix,
 * - \f$\mathbf{C}\f$ is a \f$p \times n}\f$ matrix,
 * - \f$\mathbf{D}\f$ is a \f$p \times m}\f$ matrix.
 * .
 * We define:
 * - the <em>(state) controllability matrix</em> as the following
 *   \f$n \times nm\f$ matrix:
 *   \f[
 *     \begin{pmatrix}
 *       \mathbf{B} & \mathbf{A}\mathbf{B} & \cdots & \mathbf{A}^{n-1}\mathbf{B}
 *     \end{pmatrix}
 *   \f]
 * - the <em>output controllability matrix</em> as the following
 *   \f$p \times (n+1)m\f$ matrix:
 *   \f[
 *     \begin{pmatrix}
 *       \mathbf{C}\mathbf{B} & \mathbf{C}\mathbf{A}\mathbf{B} & \cdots & \mathbf{C}\mathbf{A}^{n-1}\mathbf{B} & \mathbf{D}
 *     \end{pmatrix}
 *   \f]
 * .
 * Then, we say that:
 * - A system is <em>complete state controllable</em> (or simply
 *   \e controllable) if and only if the associated controllability matrix is of
 *   rank \f$n\f$.
 * - A system is <em>complete output controllable</em> (or simply
 *   \e output controllable) if and only if the associated output
 *   controllability matrix is of rank \f$p\f$.
 * .
 *
 * \note
 *  Estimating the rank of the controllability matrix is ill-conditioned; that
 *  is, it is very sensitive to roundoff errors and errors in the data.
 *
 * References:
 * -# K. Ogata, <em>Modern Control Engineering</em>, 3rd edition, Prentice Hall, 2002.
 * .
 *
 * \todo
 *  It might be a good idea to add an additional parameter to the function for
 *  checking controllability, representing the tolerance with which estimating
 *  the rank of the controllability matrix. However, this may be a problem if we
 *  change in the future the way we check for controllability. As a matter of
 *  fact, there are several ways for checking the controllability of a system.
 *  For instance, in Octave, controllability is determined by applying Arnoldi
 *  iteration with complete re-orthogonalization to obtain an orthogonal basis
 *  of the Krylov subspace.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2009 Marco Guazzone (marco.guazzone@gmail.com)
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

#ifndef DCS_CONTROL_ANALYSIS_CONTROLLABILITY_HPP
#define DCS_CONTROL_ANALYSIS_CONTROLLABILITY_HPP


#include <algorithm>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/eps.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/rank.hpp>
#include <boost/numeric/ublasx/operation/rot90.hpp>
#include <boost/numeric/ublasx/operation/sum.hpp>
#include <boost/numeric/ublasx/operation/svd.hpp>
#include <limits>


namespace dcs { namespace control {

/**
 * \brief Type traits for building the type of the controllability matrix given
 *  the type of the state and the input matrices.
 *
 * \tparam AMatrixT The type of the state matrix for a state-space system model.
 * \tparam BMatrixT The type of the input matrix for a state-space system model.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
struct controllability_matrix_traits
{
	/// The type of the controllability matrix.
	typedef typename ::boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type type;
};


/**
 * \brief Type traits for building the type of the output controllability matrix
 *  given the type of the state, the input, the output, and the direct
 *  transmission matrices.
 *
 * \tparam AMatrixT The type of the state matrix for a state-space system model.
 * \tparam BMatrixT The type of the input matrix for a state-space system model.
 * \tparam CMatrixT The type of the control matrix for a state-space system
 *  model.
 * \tparam DMatrixT The type of the direct transmission matrix for a state-space
 *  system model.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
struct output_controllability_matrix_traits
{
	/// The type of the output controllability matrix.
	typedef typename ::boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type type;
};


/**
 * \brief Compute the (state) controllability matrix.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return The state controllability matrix.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
typename controllability_matrix_traits<
	AMatrixT,
	BMatrixT
>::type make_controllability_matrix(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
									::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename controllability_matrix_traits<AMatrixT,BMatrixT>::type matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::value_type value_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	size_type n = ublasx::num_rows(B);
	size_type m = ublasx::num_columns(B);

	// [B AB A^2*B ... A^{n-1}*B]
	matrix_type C = ublas::zero_matrix<value_type>(n, n*m);
	ublas::subrange(C, 0, n, 0, m) = B;
	for (size_type k = 1; k < n; ++k)
	{
		size_type km = k*m;

		// C[:,k*m+1:(k+1)*m] = A*C[:,(k-1)*m+1:k*m]
		ublas::subrange(C, 0, n, km, km+m) = ublas::prod(A, ublas::subrange(C, 0, n, km-m, km));
	}

	return C;
}


/**
 * \brief Compute the output controllability matrix.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return The state controllability matrix.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
typename output_controllability_matrix_traits<
	AMatrixT,
	BMatrixT,
	CMatrixT,
	DMatrixT
>::type make_output_controllability_matrix(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
										   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
										   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
										   ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename output_controllability_matrix_traits<AMatrixT,BMatrixT,CMatrixT,DMatrixT>::type matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	size_type n = ublasx::num_rows(B);
	size_type m = ublasx::num_columns(B);
	size_type p = ublasx::num_rows(C);
	size_type nr = p;
	size_type nc = (n+1)*m;

	// [C*B C*A*B C*A^2*B ... C*A^{n-1}*B D]
	// Compute: [B A*B A^2*B ... A^{n-1}*B]
	matrix_type P = make_controllability_matrix(A, B);
	size_type xnr = ublasx::num_rows(P);
	size_type xnc = ublasx::num_columns(P);

	// Make enough space for the output controllability matrix
	P.resize(::std::max(nr, xnr), nc, true);

	// Compute: C*[B A*B A^2*B ... A^{n-1}*B] == [C*B C*A*B C*A^{2}*B C*A^{n-1}*B]
	ublas::subrange(P, 0, nr, 0, xnc) = ublas::prod(C, ublas::subrange(P, 0, xnr, 0, xnc));

	// Add D to the last m columns
	ublas::subrange(P, 0, nr, xnc, nc) = D;

	// Make P of the right size
	P.resize(nr, nc, true);

	return P;
}


/**
 * \brief Test if a system is complete state controllable.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return \c true if the state-space system \f$(\mathbf{A},\mathbf{B})\f$ is
 *  complete state controllable; \c false otherwise.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
bool is_controllable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
{
	return	::boost::numeric::ublasx::rank(
				make_controllability_matrix(A, B)
			)
			==
			::boost::numeric::ublasx::num_rows(B);
}


/**
 * \brief Test if a system is complete output controllable.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \tparam CMatrixT The type of the output matrix for the state-space system
 *  model.
 * \tparam DMatrixT The type of the direct transmission matrix for the
 *  state-space system  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \param C The output matrix for the state-space system model.
 * \param D The direct transmission matrix for the state-space system model.
 * \return \c true if the state-space system \f$(\mathbf{A},\mathbf{B},
 *  \mathbf{C},\mathbf{D})\f$ is complete output controllable; \c false
 *  otherwise.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
bool is_output_controllable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
							::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
							::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
							::boost::numeric::ublas::matrix_expression<DMatrixT> const& D)
{
	return	::boost::numeric::ublasx::rank(
				make_output_controllability_matrix(A, B, C, D)
			)
			==
			::boost::numeric::ublasx::num_rows(C);
}


template <typename ValueT>
class controllable_decomposition
{
	public: typedef ValueT value_type;
	public: typedef ::std::size_t size_type;
	public: typedef ::boost::numeric::ublas::matrix<value_type> matrix_type;
	public: typedef ::boost::numeric::ublas::vector<size_type> vector_type;


	public: controllable_decomposition()
	{
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT
			>
		controllable_decomposition(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
								   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
	{
		decompose(A, B);
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename RealT
			>
		controllable_decomposition(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
								   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
								   RealT tol)
	{
		decompose(A, B, tol);
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT
			>
		void decompose(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		//decompose(A, B, ublasx::num_rows(A)*ublas::norm_1(A)*::std::numeric_limits<value_type>::epsilon());, 
		//decompose(A, B, ::std::max(ublasx::num_rows(A), ublasx::num_columns(A))*ublasx::eps(ublas::norm_2(A)));//FIXME: this is right
		decompose(A, B, ::std::max(ublasx::num_rows(A), ublasx::num_columns(A))*ublasx::eps(value_type(ublas::norm_1(A))));//FIXME: this is wrong
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename RealT
			>
		void decompose(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					   RealT tol)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		typedef typename ublas::promote_traits<
					typename ublas::promote_traits<
						typename ublas::matrix_traits<AMatrixT>::size_type,
						typename ublas::matrix_traits<BMatrixT>::size_type
					>::promote_type,
					size_type
				>::promote_type work_size_type;
//		typedef typename ublas::promote_traits<
//					typename ublas::matrix_traits<AMatrixT>::value_type,
//					typename ublas::promote_traits<
//						typename ublas::matrix_traits<BMatrixT>::value_type,
//						typename ublas::promote_traits<
//							RealT,
//							value_type
//						>::promote_type
//					>::promote_type
//				>::promote_type work_value_type;

		size_type A_nr(ublasx::num_rows(A));
//		size_type A_nc(ublasx::num_columns(A));
//		size_type B_nr(ublasx::num_rows(B));
//		size_type B_nc(ublasx::num_columns(B));

		// Initialization
		k_ = ublas::zero_vector<size_type>(A_nr);
		matrix_type PTj(ublas::identity_matrix<value_type>(A_nr, A_nr));
		matrix_type Aj(A);
		matrix_type Bj(B);
//		size_type rojn1 = B_nc;
		size_type deltaj = 0;
//		size_type sigmaj(A_nr);

		// Body loop
		for (work_size_type j = 0; j < A_nr; ++j)
		{
			ublasx::svd_decomposition<value_type> svd;
			svd.decompose(Bj);
			matrix_type Uj(svd.U());
			matrix_type Sj(svd.S());
//			matrix_type Vj(svd.V());
			size_type Sj_nr(ublasx::num_rows(Sj));
//			size_type Sj_nc(ublasx::num_columns(Sj));

			matrix_type P(ublasx::rot90(ublas::identity_matrix<value_type>(Sj_nr, Sj_nr), 1));
			//ublasx::rot90(ublas::identity_matrix<value_type>(Sj_nr, Sj_nr));
			//matrix_type P;
			Uj = ublas::prod(Uj, P);

			matrix_type BB(ublas::prod(ublas::trans(Uj), Bj));
			size_type BB_nr(ublasx::num_rows(BB));
//			size_type BB_nc(ublasx::num_columns(BB));

			size_type roj(ublasx::rank(BB, tol));
			size_type sigmaj(BB_nr-roj);
			k_(j) = roj;
			if (roj == 0 || sigmaj == 0)
			{
				break;
			}
			//matrix_type UAU = ublas::prod(ublas::prod(ublas::trans(Uj), Aj), Uj);
			matrix_type UAU(ublas::prod(ublas::trans(Uj), Aj));
			UAU = ublas::prod(UAU, Uj);
			Aj = ublas::subrange(UAU, 0, sigmaj, 0, sigmaj);
			Bj = ublas::subrange(UAU, 0, sigmaj, sigmaj, sigmaj+roj);
			size_type Uj_nr(ublasx::num_rows(Uj));
			size_type Uj_nc(ublasx::num_columns(Uj));
			PTj = ublas::prod(
					PTj,
					ublasx::cat_columns(
						ublasx::cat_rows(Uj, ublas::zero_matrix<value_type>(Uj_nr, deltaj)),
						ublasx::cat_rows(ublas::zero_matrix<value_type>(deltaj, Uj_nc), ublas::identity_matrix<value_type>(deltaj, deltaj))
					)
				);
			deltaj += roj;
		}

		// Finalization
		T_ = ublas::trans(PTj);
		//Ab_ = ublas::prod(ublas::prod(T_, A), PTj);
		Ab_ = ublas::prod(T_, A);
		Ab_ = ublas::prod(Ab_, PTj);
		Bb_ = ublas::prod(T_, B);
	}


	/// Return the similarity transformation matrix T such that A_bar=TAT' and B_bar=TB.
	public: matrix_type const& T() const
	{
		return T_;
	}


	/// Return the vector containing the number of controllable states at each iteration
	public: vector_type const& k() const
	{
		return k_;
	}


	public: size_type num_controllable_states() const
	{
		return ::boost::numeric::ublasx::sum(k_);
	}


	/// Return the transformed A matrix.
	public: matrix_type const& A_bar() const
	{
		return Ab_;
	}


	/// Return the transformed B matrix.
	public: matrix_type const& B_bar() const
	{
		return Bb_;
	}


	public: template <typename MatrixT>
		matrix_type C_bar(::boost::numeric::ublas::matrix_expression<MatrixT> const& C) const
	{
		namespace ublas = ::boost::numeric::ublas;

		return ublas::prod(C, ublas::trans(T_));
	}


	private: matrix_type T_;
	private: vector_type k_;
	private: matrix_type Ab_;
	private: matrix_type Bb_;
};


/*
template <
	typename AMatrixT,
	typename BMatrixT
>
controllable_decomposition(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<AMatrixT>::value_type,
				typename ublas::matrix_traits<BMatrixT>::value_type
			>::promote_type value_type;

	return controllable_decomposition(A, B, ublasx::num_rows(A)*ublas::norm_1(A)*::std::numeric_limits<value_type>::epsilon());, 
}


template <
	typename AMatrixT,
	typename BMatrixT,
	typename RealT
>
controllable_decomposition(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
						   RealT tol)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<AMatrixT>::size_type,
				typename ublas::matrix_traits<BMatrixT>::size_type
			>::promote_type size_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<AMatrixT>::value_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<BMatrixT>::value_type,
					RealT
				>::promote_type
			>::promote_type value_type;
	typedef ublas::matrix<value_type> work_matrix_type;
	typedef ublas::vector<value_type> work_vector_type;

	size_type A_nr(ublasx::num_rows(A));
	size_type A_nc(ublasx::num_columns(A));
	size_type B_nr(ublasx::num_rows(B));
//	size_type B_nc(ublasx::num_columns(B));

	// Initialization
	work_matrix_type PTj(ublas::identity_matrix<value_type>(A_nr, A_nr));
	work_matrix_type Aj(A);
	work_matrix_type Bj(B);
//	size_type rojn1 = B_nc;
	size_type deltaj = 0;
//	size_type sigmaj(A_nr);
	work_vector_type k(ublas::zero_vector<value_type>(A_nr));

	// Body loop
	for (size_type j = 0; j < A_nr; ++j)
	{
		ublasx::svd_decomposition<value_type> svd;
		svd.decompose(Bj);
		work_matrix_type Uj(svd.U());
		work_matrix_type Sj(svd.S());
		work_matrix_type Vj(svd.V());
		size_type S_nr(ublasx::num_rows(S));
		size_type S_nc(ublasx::num_columns(S));

		work_matrix_type P(ublasx::rot90(ublas::identity_matrix<value_type>(S_nr, S_nr)));
		Uj = ublas::prod(Uj, P);

		work_matrix_type BB(ublas::prod(ublas::trans(Uj), Bj));
		BB_nr = ublasx::num_rows(BB);
		BB_nc = ublasx::num_columns(BB);

		size_type roj(ublasx::rank(BB, tol));
		size_type sigmaj(BB_nr-roj);
		k(jj) = roj;
		if (roj == 0 || sigmaj == 0)
		{
			break;
		}
		abxy = ublas::prod(ublas::prod(ublas::trans(Uj), Aj), Uj);
		Aj = ublas::subrange(abxy, 0, sigmaj, 0, sigmaj);
		Bj = ublas::subrange(abxy, 0, sigmaj, sigmaj, sigmaj+roj);
		size_type Uj_nr(ublasx::num_rows(Uj));
		size_type Uj_nc(ublasx::num_columns(Uj));
		PTj = ublas::prod(
				PTj,
				ublasx::cat_columns(
					ublasx::cat_rows(U, ublas::zero_matrix<value_type>(Uj_nr, deltaj)),
					ublasx::cat_rows(ublas::zero_matrix<value_type>(deltaj, Uj_nc, ublas::identity_matrix<value_type>(deltaj, deltaj)))
				)
			);
		deltaj += roj;
	}

	// Finalization
	T = ublas::trans(PTj);
	//A_bar = ublas::prod(ublas::prod(T, a), ublas::trans(T));
	A_bar = ublas::prod(ublas::prod(t, a), PTj);
	B_bar = ublas::prod(T, b);
	//cbar = ublas::prod(c, ublas::trans(T));
	C_bar = ublas::prod(c, PTj);
}
*/

}} // Namespace dcs::control


#endif // DCS_CONTROL_ANALYSIS_CONTROLLABILITY_HPP
