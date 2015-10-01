/**
 * \file dcs/control/solver/qp.hpp
 *
 * \brief Quadratic Programming (QP) problem solver
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

#ifndef DCS_CONTROL_SOLVER_QP_HPP
#define DCS_CONTROL_SOLVER_QP_HPP

#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <cstddef>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <dcs/exception.hpp>
#include <dcs/logging.hpp>
#include <sstream>
#include <stdexcept>

#if defined(DCS_CONTROL_QP_USE_CPLEX)
# include <ilconcert/iloalg.h>
# include <ilconcert/iloenv.h>
# include <ilconcert/iloexpression.h>
# include <ilconcert/ilomodel.h>
# include <ilcplex/ilocplex.h>
#elif defined(DCS_CONTROL_QP_USE_GUROBI)
# include <gurobi_c++.h>
# include <vector>
#elif defined(DCS_CONTROL_QP_USE_QPOASES)
# include <boost/smart_ptr.hpp>
# include <qpOASES.hpp>
#else
# error Unable to find a suitable QP solver
#endif // DCS_CONTROL_QP_USE_CPLEX


namespace dcs { namespace control {

namespace detail {

#ifdef DCS_CONTROL_QP_USE_CPLEX
template <typename RealT,
		  typename QMatrixT,
		  typename CVectorT,
		  typename AMatrixT,
		  typename BVectorT>//,
		  //typename RealT>
boost::numeric::ublas::vector<RealT> qp_solve_by_cplex(const boost::numeric::ublas::matrix_expression<QMatrixT>& Q,
													   const boost::numeric::ublas::vector_expression<CVectorT>& c,
													   const boost::numeric::ublas::matrix_expression<AMatrixT>& A,
													   const boost::numeric::ublas::vector_expression<BVectorT>& b)
{
	namespace ublas = boost::numeric::ublas;
	namespace ublasx = boost::numeric::ublasx;

	ublas::vector<RealT> sol;

	DCS_ASSERT(ublasx::num_rows(Q) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q is not a square matrix"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(b) == ublasx::num_rows(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector b have incompatible dimensions"));

	const std::size_t nb = ublasx::size(b);
	const std::size_t nx = ublasx::size(c);

	try
	{
		// Initialize the Concert Technology app
		IloEnv env;

		IloModel model(env);
		model.setName("MPC - QP solver");

		// Decision variables

		// Variable x_{i}
		IloNumVarArray x(env, nx);
		for (std::size_t i = 0; i < nx; ++i)
		{
			std::ostringstream oss;
			oss << "x[" << i << "]";
			x[i] = IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT, oss.str().c_str());
			model.add(x[i]);
		}

		// Constraints

		std::size_t cc = 0; // Constraint counter
		std::size_t csc = 0; // Constraint subcounter

		// C1: Ax \le b
		++cc;
		for (std::size_t i = 0; i < nb; ++i)
		{
/*
			++csc;

			std::ostringstream oss;
			oss << "C" << cc << "_{" << csc << "}";

			IloExpr lhs(env);
			for (std::size_t j = 0; j < nx; ++j)
			{
				lhs += A()(i,j)*x[j];
			}

			IloConstraint cons;
			if (std::isfinite(b()(i)))
			{
				cons = (lhs <= b()(i));
			}
			else
			{
				if (b()(i) > 0)
				{
					cons = (lhs <= IloInfinity);
				}
				else
				{
					cons = (lhs <= -IloInfinity);
				}
			}
			cons.setName(oss.str().c_str());
			model.add(cons);
*/
			if (std::isfinite(b()(i)))
			{
				++csc;

				std::ostringstream oss;
				oss << "C" << cc << "_{" << csc << "}";

				IloExpr lhs(env);
				for (std::size_t j = 0; j < nx; ++j)
				{
					lhs += A()(i,j)*x[j];
				}

				IloConstraint cons(lhs <= b()(i));
				model.add(cons);
			}
		}

		// Objective
		IloObjective z;
		{
			IloExpr expr(env);
			for (std::size_t i = 0; i < nx; ++i)
			{
				// Quadratic term
				for (std::size_t j = 0; j < nx; ++j)
				{
					expr += 0.5*x[i]*Q()(i,j)*x[j];
				}
				// Linear term
				expr += c()(i)*x[i];
			}
			z = IloMinimize(env, expr);
		}
		model.add(z);

		// Solve the model.
		// NOTE: The default solver (i.e., IloCplex::RootAlg parameter set to
		//  IloCplex::Auto) seems to have troubles maybe due to numerical
		//  problems.
		//  The primal simplex algorithm (i.e., IloCplex::IloPrimal) seems to be
		//  more robust, instead.

		IloCplex solver(model);
		//solver.setParam(IloCplex::RootAlg, IloCplex::Primal);
#ifndef DCS_DEBUG
		solver.setOut(env.getNullStream());
		solver.setWarning(env.getNullStream());
#else // DCS_DEBUG
		solver.exportModel("cplex-mpc_qp.sav");
#endif // DCS_DEBUG

		const bool solved = solver.solve();

		IloAlgorithm::Status status = solver.getStatus();
		switch (status)
		{
			case IloAlgorithm::Optimal: // The algorithm found an optimal solution.
			case IloAlgorithm::Feasible: // The algorithm found a feasible solution, though it may not necessarily be optimal.
				break;
			case IloAlgorithm::Infeasible: // The algorithm proved the model infeasible (i.e., it is not possible to find an assignment of values to variables satisfying all the constraints in the model).
			case IloAlgorithm::Unbounded: // The algorithm proved the model unbounded.
			case IloAlgorithm::InfeasibleOrUnbounded: // The model is infeasible or unbounded.
			case IloAlgorithm::Error: // An error occurred and, on platforms that support exceptions, that an exception has been thrown.
			case IloAlgorithm::Unknown: // The algorithm has no information about the solution of the model.
			{
				::std::ostringstream oss;
				oss << "Optimization was stopped with status = " << status << " (CPLEX status = " << solver.getCplexStatus() << ", sub-status = " << solver.getCplexSubStatus() << ")";
				dcs::log_warn(DCS_LOGGING_AT, oss.str());
			}
		}

		if (solved)
		{
#ifdef DCS_DEBUG
			DCS_DEBUG_TRACE( "-------------------------------------------------------------------------------[" );
			DCS_DEBUG_TRACE( "- Objective value: " << static_cast<RealT>(solver.getObjValue()) );

			DCS_DEBUG_TRACE( "- Decision variables: " );
			// Output x_{i}
			for (std::size_t i = 0; i < nx; ++i)
			{
				DCS_DEBUG_STREAM << x[i].getName() << " = " << solver.getValue(x[i]) << ::std::endl;
			}

			DCS_DEBUG_TRACE( "]-------------------------------------------------------------------------------" );
#endif // DCS_DEBUG

			sol.resize(nx);
			for (std::size_t i = 0; i < nx; ++i)
			{
				sol[i] = static_cast<RealT>(solver.getValue(x[i]));
			}
		}

		z.end();
		x.end();

		// Close the Concerty Technology app
		env.end();
	}
	catch (IloException const& e)
	{
		std::ostringstream oss;
		oss << "Got exception from CPLEX: " << e.getMessage();
		DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
	}
	catch (...)
	{
		DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected error during the optimization");
	}

	return sol;
}
#endif // DCS_CONTROL_QP_USE_CPLEX

#ifdef DCS_CONTROL_QP_USE_GUROBI
template <typename RealT,
		  typename QMatrixT,
		  typename CVectorT,
		  typename AMatrixT,
		  typename BVectorT>//,
		  //typename RealT>
boost::numeric::ublas::vector<RealT> qp_solve_by_gurobi(const boost::numeric::ublas::matrix_expression<QMatrixT>& Q,
														const boost::numeric::ublas::vector_expression<CVectorT>& c,
														const boost::numeric::ublas::matrix_expression<AMatrixT>& A,
														const boost::numeric::ublas::vector_expression<BVectorT>& b)
{
	namespace ublas = boost::numeric::ublas;
	namespace ublasx = boost::numeric::ublasx;

	ublas::vector<RealT> sol;

	DCS_ASSERT(ublasx::num_rows(Q) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q is not a square matrix"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(b) == ublasx::num_rows(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector b have incompatible dimensions"));

	const std::size_t nb = ublasx::size(b);
	const std::size_t nx = ublasx::size(c);

	try
	{
		GRBEnv env = GRBEnv();

		GRBModel model = GRBModel(env);

		model.set(GRB_StringAttr_ModelName, "MPC - QP solver");

#ifdef DCS_DEBUG
		env.set(GRB_IntParam_OutputFlag, 1);
		env.set(GRB_IntParam_LogToConsole, 1);
#else
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_LogToConsole, 0);
#endif // DCS_DEBUG

		// Decision variables

		// Variable x_{i}
		std::vector<GRBVar> x(nx);
		for (std::size_t i = 0; i < nx; ++i)
		{
			std::ostringstream oss;
			oss << "x[" << i << "]";
			x[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, oss.str());
		}

		// Integrate variables into the model
		model.update();

		// Constraints

		std::size_t cc = 0; // Constraint counter
		std::size_t csc = 0; // Constraint subcounter

		// C1: Ax \le b
		++cc;
		for (std::size_t i = 0; i < nb; ++i)
		{
			++csc;

			std::ostringstream oss;
			oss << "C" << cc << "_{" << csc << "}";

			GRBLinExpr lhs;
			for (std::size_t j = 0; j < nx; ++j)
			{
				lhs += A()(i,j)*x[j];
			}

			if (std::isfinite(b()(i)))
			{
				model.addConstr(lhs, GRB_LESS_EQUAL, b()(i), oss.str());
			}
			else
			{
				if (b()(i) > 0)
				{
					model.addConstr(lhs, GRB_LESS_EQUAL, GRB_INFINITY, oss.str());
				}
				else
				{
					model.addConstr(lhs, GRB_LESS_EQUAL, -GRB_INFINITY, oss.str());
				}
			}
		}

		// Objective
		GRBQuadExpr obj;
		for (std::size_t i = 0; i < nx; ++i)
		{
			// Quadratic term
			for (std::size_t j = 0; j < nx; ++j)
			{
				obj += 0.5*x[i]*Q()(i,j)*x[j];
			}
			// Linear term
			obj += c()(i)*x[i];
		}
		model.setObjective(obj, GRB_MINIMIZE);

		// Integrate variables into the model
		model.update();

#ifdef DCS_DEBUG
		model.write("gurobi-mpc_qp.lp");
#endif // DCS_DEBUG

		model.optimize();

		bool solved = false;

		int status = model.get(GRB_IntAttr_Status);
		switch (status)
		{
			case GRB_OPTIMAL: // The algorithm found an optimal solution.
			case GRB_SUBOPTIMAL: // The algorithm found a feasible solution, though it may not necessarily be optimal.
				solved = true;
				break;
			default:
				{
					::std::ostringstream oss;
					oss << "Optimization was stopped with status = " << status;
					dcs::log_warn(DCS_LOGGING_AT, oss.str());
				}
		}

		if (solved)
		{
#ifdef DCS_DEBUG
			DCS_DEBUG_TRACE( "-------------------------------------------------------------------------------[" );
			DCS_DEBUG_TRACE( "- Objective value: " << static_cast<RealT>(model.get(GRB_DoubleAttr_ObjVal)) );

			DCS_DEBUG_TRACE( "- Decision variables: " );
			// Output x_{i}
			for (std::size_t i = 0; i < nx; ++i)
			{
				DCS_DEBUG_STREAM << x[i].get(GRB_StringAttr_VarName) << " = " << x[i].get(GRB_DoubleAttr_X) << " (" << static_cast<RealT>(x[i].get(GRB_DoubleAttr_X)) << ")" << ::std::endl;
			}

			DCS_DEBUG_TRACE( "]-------------------------------------------------------------------------------" );
#endif // DCS_DEBUG

			sol.resize(nx);
			for (std::size_t i = 0; i < nx; ++i)
			{
				sol[i] = static_cast<RealT>(x[i].get(GRB_DoubleAttr_X));
			}
		}
	}
	catch (GRBException const& e)
	{
		std::ostringstream oss;
		oss << "Got exception from GUROBI: " << e.getMessage() << " (Error code: " << e.getErrorCode() << ")";
		DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
	}
	catch (...)
	{
		DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected error during the optimization");
	}

	return sol;
}
#endif // DCS_CONTROL_QP_USE_GUROBI

#ifdef DCS_CONTROL_QP_USE_QPOASES
template <typename RealT,
		  typename QMatrixT,
		  typename CVectorT,
		  typename AMatrixT,
		  typename BVectorT>//,
		  //typename RealT>
boost::numeric::ublas::vector<RealT> qp_solve_by_qpoases(const boost::numeric::ublas::matrix_expression<QMatrixT>& Q,
														const boost::numeric::ublas::vector_expression<CVectorT>& c,
														const boost::numeric::ublas::matrix_expression<AMatrixT>& A,
														const boost::numeric::ublas::vector_expression<BVectorT>& b)
{
	namespace ublas = boost::numeric::ublas;
	namespace ublasx = boost::numeric::ublasx;

	ublas::vector<RealT> sol;

	DCS_ASSERT(ublasx::num_rows(Q) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q is not a square matrix"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(Q),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix Q and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(c) == ublasx::num_columns(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector c have incompatible dimensions"));
	DCS_ASSERT(ublasx::size(b) == ublasx::num_rows(A),
			   DCS_EXCEPTION_THROW(std::invalid_argument, "Matrix A and vector b have incompatible dimensions"));

	const std::size_t nb = ublasx::size(b);
	const std::size_t nx = ublasx::size(c);

	try
	{
		qpOASES::QProblem problem(nx, nb);

		qpOASES::Options options;
#ifdef DCS_DEBUG
		options.printLevel = qpOASES::PL_DEBUG_ITER;
#else
		options.printLevel = qpOASES::PL_NONE;
#endif // DCS_DEBUG
		problem.setOptions(options);

		boost::scoped_array<qpOASES::real_t> tmp_Q(new qpOASES::real_t[nx*nx]);
		for (std::size_t i = 0; i < nx; ++i)
		{
			const std::size_t offs = i*nx;
			for (std::size_t j = 0; j < nx; ++j)
			{
				tmp_Q[offs+j] = Q()(i,j);
			}
		}
		boost::scoped_array<qpOASES::real_t> tmp_c(new qpOASES::real_t[nx]);
		for (std::size_t i = 0; i < nx; ++i)
		{
			tmp_c[i] = c()(i);
		}
		boost::scoped_array<qpOASES::real_t> tmp_A(new qpOASES::real_t[nb*nx]);
		for (std::size_t i = 0; i < nb; ++i)
		{
			const std::size_t offs = i*nx;
			for (std::size_t j = 0; j < nx; ++j)
			{
				tmp_A[offs+j] = A()(i,j);
			}
		}
		boost::scoped_array<qpOASES::real_t> tmp_b(new qpOASES::real_t[nb]);
		for (std::size_t i = 0; i < nb; ++i)
		{
			tmp_b[i] = b()(i);
		}
		boost::scoped_array<qpOASES::real_t> tmp_lbA(new qpOASES::real_t[nb]);
		for (std::size_t i = 0; i < nb; ++i)
		{
			tmp_lbA[i] = -std::numeric_limits<qpOASES::real_t>::infinity();
		}

		int nWSR = 10;

		qpOASES::returnValue ret;
		ret = problem.init(tmp_Q.get(),
						   tmp_c.get(),
						   tmp_A.get(),
						   0,
						   0,
						   tmp_lbA.get(),
						   tmp_b.get(),
						   nWSR);

		if (ret != qpOASES::SUCCESSFUL_RETURN)
		{
			std::ostringstream oss;
			oss << "Unable to initialize the qpOASES problem (status: " << ret << ")";
			DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
		}

		boost::scoped_array<qpOASES::real_t> x(new qpOASES::real_t[nx]);
		problem.getPrimalSolution(x.get());

#ifdef DCS_DEBUG
		DCS_DEBUG_TRACE( "-------------------------------------------------------------------------------[" );
		DCS_DEBUG_TRACE( "- Objective value: " << static_cast<RealT>(problem.getObjVal()) );

		DCS_DEBUG_TRACE( "- Decision variables: " );
		// Output x_{i}
		for (std::size_t i = 0; i < nx; ++i)
		{
			DCS_DEBUG_STREAM << "x[" << i << "]" << " = " << x[i] << " (" << static_cast<RealT>(x[i]) << ")" << ::std::endl;
		}

		DCS_DEBUG_TRACE( "]-------------------------------------------------------------------------------" );
#endif // DCS_DEBUG

		sol.resize(nx);
		for (std::size_t i = 0; i < nx; ++i)
		{
			sol[i] = static_cast<RealT>(x[i]);
		}
	}
	catch (std::exception const& e)
	{
		std::ostringstream oss;
		oss << "Got exception from qpOASES: " << e.what();
		DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
	}
	catch (...)
	{
		DCS_EXCEPTION_THROW(std::runtime_error, "Unexpected error during the optimization");
	}

	return sol;
}
#endif // DCS_CONTROL_QP_USE_QPOASES

} // Namespace detail


/**
 * Solve a quadratic programming (QP) problem.
 *
 * Solve the following quadratic programming problem:
 * \f{align}
 *  \text{minimize} &\quad \frac{1}{2} x^T Q x + c^T x,
 *  \text{subject to} &\quad Ax \le b.
 * \f}
 *
 * \param Q A real symmetric matrix representing the quadratic term in the objective expression
 * \param c A real vector representing the linear term in the objective expression
 * \param A A real matrix representing the linear coefficients in the constraints \f$Ax \le b\f$
 * \param b A real vector representing the constant vector in the constraints \f$Ax \le b\f$
 * \return The vector \f$x\f$ that solves the quadratic programming problem.
 *
 * \note Soft constraints are not supported yet.
 */
template <typename RealT,
		  typename QMatrixT,
		  typename CVectorT,
		  typename AMatrixT,
		  typename BVectorT>//,
		  //typename RealT>
boost::numeric::ublas::vector<RealT> qp_solve(const boost::numeric::ublas::matrix_expression<QMatrixT>& Q,
											  const boost::numeric::ublas::vector_expression<CVectorT>& c,
											  const boost::numeric::ublas::matrix_expression<AMatrixT>& A,
											  const boost::numeric::ublas::vector_expression<BVectorT>& b)
{
#if defined(DCS_CONTROL_QP_USE_CPLEX)
	return detail::qp_solve_by_cplex<RealT>(Q, c, A, b);
#elif defined(DCS_CONTROL_QP_USE_GUROBI)
	return detail::qp_solve_by_gurobi<RealT>(Q, c, A, b);
#elif defined(DCS_CONTROL_QP_USE_QPOASES)
	return detail::qp_solve_by_qpoases<RealT>(Q, c, A, b);
#else
# error Unable to find a suitable QP solver
#endif // DCS_CONTROL_QP_USE_...
}

}} // Namespace dcs::control

#endif // DCS_CONTROL_SOLVER_QP_HPP
