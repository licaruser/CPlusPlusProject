#pragma once
/*
* Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 2 or any later version.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
* more details. A copy of the GNU General Public License is available at:
* http://www.fsf.org/licensing/licenses
*/


/*****************************************************************************
*                                    random.h
*
* Random Number Generator.
*
* The class Random can generate an "long int" random number. The default
* initial state is set to 1, of cause it can be set to any integer through
* routine "seed(int)".
*
* Based on the "Random" class, we can generate some usually used distributed
* pseudorandom number or pseudorandom sequence. These distributions include:
*      Uniform     Normal      Exponential     Rayleigh        Poisson
*      Bernoulli
*
* Zhang Ming, 2010-10, Xi'an Jiaotong University.
*****************************************************************************/

#ifndef RANDOMGEN_H
#define RANDOMGEN_H


#include <vector.h>


namespace splab
{

	class Random
	{

	public:

		explicit Random(long int initValue = 1)
		{
			if (initValue < 0)
				initValue += M;

			state = initValue;
			if (state == 0)
				state = 1;
		};
		~Random() {};

		void seed(long int value);
		long int random();
		long int getM() const;

	private:

		static const long int A = 48271;
		static const long int M = 2147483647;
		static const long int Q = 44488;
		static const long int R = 3399;

		long int state;

	};

	template<typename Type> Type randu(int, const Type&, const Type&);
	template<typename Type> Type randn(int, const Type&, const Type&);
	template<typename Type> Type rande(int, const Type&);
	template<typename Type> Type randr(int, const Type&);
	template<typename Type> int randp(int, const Type&);
	template<typename Type> int randb(int, const Type&);

	template<typename Type> Vector<Type> randu(int, const Type&,
		const Type&, int);
	template<typename Type> Vector<Type> randn(int, const Type&,
		const Type&, int);
	template<typename Type> Vector<Type> rande(int, const Type&, int);
	template<typename Type> Vector<Type> randr(int, const Type&, int);
	template<typename Type> Vector<int>  randp(int, const Type&, int);
	template<typename Type> Vector<int>  randb(int, const Type&, int);

	/**
	* Set the internal state.
	*/
	inline void Random::seed(long int value)
	{
		state = value;
	}


	/**
	* Return a pseudorandom int, and then change the internal state.
	*/
	inline long int Random::random()
	{
		long int tmpState = A*(state%Q) - R*(state / Q);

		if (tmpState >= 0)
			state = tmpState;
		else
			state = tmpState + M;

		return state;
	}


	/**
	* Return M.
	*/
	inline long int Random::getM() const
	{
		return M;
	}


	/**
	* Return an Uniform[low, high] distributed pseudorandom number.
	*/
	template <typename Type>
	Type randu(int seed, const Type &low, const Type &high)
	{
		static Random rg(seed);
		long double u_0_1 = rg.random() / (long double)rg.getM();

		return low + Type(u_0_1*(high - low));
	}


	/**
	* Return a Normal~(mu, sigma) distributed pseudorandom number.
	*/
	template <typename Type>
	Type randn(int seed, const Type &mu, const Type &sigma)
	{
		static Random rg(seed);
		long double u_0_1 = 0.0;

		for (int i = 0; i<12; ++i)
			u_0_1 += rg.random() / (long double)rg.getM();
		u_0_1 -= 6.0;

		return Type(sigma*u_0_1 + mu);
	}


	/**
	* Return a Exponential~(beta) distributed pseudorandom number.
	*/
	template <typename Type>
	Type rande(int seed, const Type &beta)
	{
		static Random rg(seed);
		long double u_0_1 = rg.random() / (long double)rg.getM();

		return Type(-beta*log(u_0_1));
	}


	/**
	* Return a Rayleigh~(sigma) distributed pseudorandom number.
	*/
	template <typename Type>
	Type randr(int seed, const Type &sigma)
	{
		static Random rg(seed);
		long double u_0_1 = rg.random() / (long double)rg.getM();

		return Type(sigma * sqrt(-2.0*log(u_0_1)));
	}


	/**
	* Return a Poisson~(lambda) distributed pseudorandom number.
	*/
	template <typename Type>
	int randp(int seed, const Type &lambda)
	{
		static Random rg(seed);
		long double u_0_1 = 0.0,
			q = exp(-lambda),
			p = 1.0;
		int k = 0;

		while (p >= q)
		{
			u_0_1 = rg.random() / (long double)rg.getM();
			p *= u_0_1;
			k++;
		}
		return k - 1;
	}


	/**
	* Return a Bernoulli~(p) distributed pseudorandom number.
	*/
	template <typename Type>
	int randb(int seed, const Type &p)
	{
		static Random rg(seed);
		long double u_0_1 = rg.random() / (long double)rg.getM();

		return (Type(u_0_1) <= p) ? 1 : 0;
	}


	/**
	* Return an Uniform[low, high] distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<Type> randu(int seed, const Type &low, const Type &high, int N)
	{
		long double u_0_1 = 0.0;
		Random rg(seed);
		Vector<Type> rs(N);

		for (int i = 0; i<N; ++i)
		{
			u_0_1 = rg.random() / (long double)rg.getM();
			rs[i] = low + Type(u_0_1*(high - low));
		}

		return rs;
	}


	/**
	* Return an Normal~(mu, sigma) distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<Type> randn(int seed, const Type &mu, const Type &sigma, int N)
	{
		long double u_0_1;
		Random rg(seed);
		Vector<Type> rs(N);

		for (int i = 0; i<N; ++i)
		{
			u_0_1 = 0.0;
			for (int j = 0; j<12; ++j)
				u_0_1 += rg.random() / (long double)rg.getM();
			u_0_1 -= 6.0;

			rs[i] = Type(sigma*u_0_1) + mu;
		}

		return rs;
	}


	/**
	* Return an Exponential~(beta) distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<Type> rande(int seed, const Type &beta, int N)
	{
		long double u_0_1 = 0.0;
		Random rg(seed);
		Vector<Type> rs(N);

		for (int i = 0; i<N; ++i)
		{
			u_0_1 = rg.random() / (long double)rg.getM();
			rs[i] = Type(-beta*log(u_0_1));
		}

		return rs;
	}


	/**
	* Return an Rayleigh~(sigma) distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<Type> randr(int seed, const Type &sigma, int N)
	{
		long double u_0_1 = 0.0;
		Random rg(seed);
		Vector<Type> rs(N);

		for (int i = 0; i<N; ++i)
		{
			u_0_1 = rg.random() / (long double)rg.getM();
			rs[i] = Type(sigma * sqrt(-2.0*log(u_0_1)));
		}

		return rs;
	}


	/**
	* Return a Poisson~(lambda) distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<int> randp(int seed, const Type &lambda, int N)
	{
		Random rg(seed);
		Vector<int> rs(N);

		long double u_0_1 = 0.0,
			q = exp(-lambda),
			p = 1.0;
		int k = 0;

		for (int i = 0; i<N; ++i)
		{
			k = 0;
			p = 1.0;
			while (p >= q)
			{
				u_0_1 = rg.random() / (long double)rg.getM();
				p *= u_0_1;
				k++;
			}
			rs[i] = k - 1;
		}

		return rs;
	}


	/**
	* Return a Bernoulli~(p) distributed pseudorandom sequence of size N.
	*/
	template <typename Type>
	Vector<int> randb(int seed, const Type &p, int N)
	{
		long double u_0_1 = 0.0;
		Random rg(seed);
		Vector<int> rs(N);

		for (int i = 0; i<N; ++i)
		{
			u_0_1 = rg.random() / (long double)rg.getM();
			rs[i] = (Type(u_0_1) <= p) ? 1 : 0;
		}

		return rs;
	}

}
// namespace splab


#endif //RANDOMGEN_H
