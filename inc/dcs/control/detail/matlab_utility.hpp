/**
 * \file dcs/control/detail/matlab_utility.hpp
 *
 * \brief Utilities to interface with the MATLAB environment.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * \copyright
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

#ifndef DCS_CONTROL_DETAIL_MATLAB_UTILITY_HPP
#define DCS_CONTROL_DETAIL_MATLAB_UTILITY_HPP


#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <cctype>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <dcs/exception.hpp>
#include <dcs/system/posix_process.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace dcs { namespace control { namespace detail {

template <typename ValueT>
std::string to_matlab_str(const std::vector<ValueT>& v, bool column = true)
{
	typedef std::vector<ValueT> vector_type;
	typedef typename vector_type::size_type size_type;

	std::ostringstream oss;
	oss << "[";
	for (size_type i = 0,
				   n = v.size();
		 i < n;
		 ++i)
	{
		if (i > 0)
		{
			if (column)
			{
				oss << ";";
			}
			else
			{
				oss << " ";
			}
		}
		oss << v[i];
	}
	oss << "]";

	return oss.str();
}

template <typename VectorExprT>
std::string to_matlab_str(const boost::numeric::ublas::vector_expression<VectorExprT>& v, bool column = true)
{
	namespace ublas = boost::numeric::ublas;
	namespace ublasx = boost::numeric::ublasx;

	typedef VectorExprT vector_type;
	typedef typename ublas::vector_traits<vector_type>::size_type size_type;

	std::ostringstream oss;
	oss << "[";
	for (size_type i = 0,
				   n = ublasx::size(v);
		 i < n;
		 ++i)
	{
		if (i > 0)
		{
			if (column)
			{
				oss << ";";
			}
			else
			{
				oss << " ";
			}
		}
		oss << v()(i);
	}
	oss << "]";

	return oss.str();
}

template <typename MatrixExprT>
std::string to_matlab_str(const boost::numeric::ublas::matrix_expression<MatrixExprT>& A)
{
	namespace ublas = boost::numeric::ublas;
	namespace ublasx = boost::numeric::ublasx;

	typedef MatrixExprT matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	const size_type nr = ublasx::num_rows(A);
	const size_type nc = ublasx::num_columns(A);
	std::ostringstream oss;
	oss << "[";
	for (size_type r = 0; r < nr; ++r)
	{
		if (r > 0)
		{
			oss << ";";
		}
		for (size_type c = 0; c < nc; ++c)
		{
			if (c > 0)
			{
				oss << " ";
			}
			oss << A()(r,c);
		}
	}
	oss << "]";

	return oss.str();
}

template <typename T>
void parse_matlab_str(const std::string& text, T& x)
{
	std::istringstream iss(text);
	while (iss.good())
	{
		bool ok = true;
		char ch = iss.peek();

		if (std::isspace(ch))
		{
			// Skip space
			iss.get();
		}
		else if (std::isdigit(ch) || ch == '+' || ch == '-' || ch == '.')
		{
			// Found the beginning of a number

			iss >> x;
			if (iss.fail())
			{
				ok = false;
			}
		}
		else
		{
			ok = false;
		}

		if (!ok)
		{
			// Look for special numbers like NaN and +/-inf
			bool neg = (ch == '-');
			if (ch == '+' || ch == '-')
			{
				iss.get();
			}
			std::string s(3, 0);
			iss.readsome(&*s.begin(), s.length());
			if (s == "NaN")
			{
				x = std::numeric_limits<T>::quiet_NaN();
				ok = true;
			}
			else if (s == "inf")
			{
				if (neg)
				{
					x = -std::numeric_limits<T>::infinity();
				}
				else
				{
					x = +std::numeric_limits<T>::infinity();
				}
				ok = true;
			}
		}
		if (!ok)
		{
			DCS_EXCEPTION_THROW(std::runtime_error, "Unable to parse a MATLAB number");
		}
	}
}

template <typename T>
void parse_matlab_str(const std::string& text, boost::numeric::ublas::vector<T>& v)
{
	typename boost::numeric::ublas::vector<T>::size_type n(0);

	std::istringstream iss(text);
	bool inside = false;
	bool done = false;
	while (iss.good() && !done)
	{
		char ch = iss.peek();
		bool ko = false;

		if (inside)
		{
			if (std::isspace(ch) || ch == ';')
			{
				// Found an element separator
				iss.get();
//					while (iss.good() && (ch = iss.peek()) && std::isspace(ch))
//					{
//						iss.get();
//					}
			}
			else if (ch == ']')
			{
				// Found the end of the vector
//					iss.get();
//					inside = false;
				done = true;
			}
			else if (std::isdigit(ch) || ch == '+' || ch == '-' || ch == '.')
			{
				// Found the beginning of a number
				T x;
				iss >> x;
				if (iss.fail())
				{
					ko = true;
				}
				else
				{
					v.resize(n+1, true);
					v(n) = x;
					++n;
				}
			}
			else
			{
				ko = true;
			}
			if (ko)
			{
				// Look for special numbers like NaN and +/-inf
				bool neg = (ch == '-');
				if (ch == '+' || ch == '-')
				{
					iss.get();
				}
				std::string s(3, 0);
				T x = 0;
				iss.readsome(&*s.begin(), s.length());
				if (s == "NaN")
				{
					x = std::numeric_limits<T>::quiet_NaN();
					ko = false;
				}
				else if (s == "inf")
				{
					if (neg)
					{
						x = -std::numeric_limits<T>::infinity();
					}
					else
					{
						x = +std::numeric_limits<T>::infinity();
					}
					ko = false;
				}
				if (!ko)
				{
					v.resize(n+1, true);
					v(n) = x;
					++n;
				}
			}
		}
		else
		{
			if (ch == '[')
			{
				iss.get();
				v.resize(0, false);
				inside = true;
			}
			else if (std::isspace(ch))
			{
				iss.get();
			}
			else
			{
				ko = true;
			}
		}

		if (ko)
		{
			DCS_EXCEPTION_THROW(std::runtime_error, "Unable to parse a MATLAB vector");
		}
	}
}

template <typename T>
void parse_matlab_str(const std::string& text, boost::numeric::ublas::matrix<T>& A)
{
	typename boost::numeric::ublas::matrix<T>::size_type r(0);
	typename boost::numeric::ublas::matrix<T>::size_type c(0);
	typename boost::numeric::ublas::matrix<T>::size_type nc(0);

	std::istringstream iss(text);
	bool inside = false;
	bool done = false;
	while (iss.good() && !done)
	{
		char ch = iss.peek();
		bool ko = false;

		if (inside)
		{
			if (std::isspace(ch))
			{
				// Found a column separator
				iss.get();
//					while (iss.good() && (ch = iss.peek()) && std::isspace(ch))
//					{
//						iss.get();
//					}
			}
			else if (ch == ';')
			{
				// Found a row separator
				iss.get();
				++r;
				c = 0;
			}
			else if (ch == ']')
			{
				// Found the end of the matrix
//					iss.get();
//					inside = false;
				done = true;
			}
			else if (std::isdigit(ch) || ch == '+' || ch == '-' || ch == '.')
			{
				// Found the beginning of a number
				T x;
				iss >> x;
				if (iss.fail())
				{
					ko = true;
				}
				else
				{
					if (nc <= c)
					{
						nc = c+1;
					}
					A.resize(r+1, nc, true);
					A(r,c) = x;
					++c;
				}
			}
			else
			{
				ko = true;
			}
			if (ko)
			{
				// Look for special numbers like NaN and +/-inf
				bool neg = (ch == '-');
				if (ch == '+' || ch == '-')
				{
					iss.get();
				}
				std::string s(3, 0);
				T x = 0;
				iss.readsome(&*s.begin(), s.length());
				if (s == "NaN")
				{
					x = std::numeric_limits<T>::quiet_NaN();
					ko = false;
				}
				else if (s == "inf")
				{
					if (neg)
					{
						x = -std::numeric_limits<T>::infinity();
					}
					else
					{
						x = +std::numeric_limits<T>::infinity();
					}
					ko = false;
				}
				if (!ko)
				{
					if (nc <= c)
					{
						nc = c+1;
					}
					A.resize(r+1, nc, true);
					A(r,c) = x;
					++c;
				}
			}
		}
		else
		{
			// Note: outside of a matrix, only two types of character are
			// allowed: spaces and '['

			if (ch == '[')
			{
				iss.get();
				A.resize(0, 0, false);
				inside = true;
			}
			else if (std::isspace(ch))
			{
				iss.get();
			}
			else
			{
				ko = true;
			}
		}

		if (ko)
		{
			DCS_EXCEPTION_THROW(std::runtime_error, "Unable to parse a MATLAB matrix");
		}
	}
}

template <typename ArgsIterT, typename ConsumerT>
void run_matlab(std::string const& cmd,
				ArgsIterT args_first, ArgsIterT args_last,
				ConsumerT& consumer)
{
	dcs::system::posix_process process(cmd);
	process.asynch(true);
	process.run(args_first, args_last, false, true);
	if (process.status() != dcs::system::running_process_status)
	{
		std::ostringstream oss;
		oss << "Unable to start MATLAB: " << strerror(errno);
		DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
	}

	consumer(process);

	process.wait();

	if (process.status() != dcs::system::terminated_process_status)
	{
		std::ostringstream oss;
		oss << "MATLAB terminated unexpectedly: " << strerror(errno);

		DCS_EXCEPTION_THROW(std::runtime_error, oss.str());
	}
}

}}} // Namespace dcs::control::detail


#endif // DCS_CONTROL_DETAIL_MATLAB_UTILITY_HPP
