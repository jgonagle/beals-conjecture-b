/*
 * BealsConjectureII.c
 *
 *  Created on: Jul 28, 2011
 *      Author: John McGonagle
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include "beals_conjecture_ii.h"

#define PRINT_DEBUG

int main()
{
	do_precompute();

	return 0;
}

void do_precompute()
{
	num_primes = sizeof(primes) / sizeof(primes[0]);
	num_positive_be_quad_combos_low = 1;
	num_positive_be_quad_combos_high = 1;
	num_negative_be_quad_combos_low = 1;
	num_negative_be_quad_combos_high = 1;
	max_congruency_class_size = 0;
	power_two_num_primes = power_two_exp(num_primes);
	power_two_biggest_prime = power_two_exp(primes[num_primes - 1]);
	power_two_biggest_prime_squared = power_two_exp(pow(primes[num_primes - 1], 2));
	power_two_max_congruency_class_size_squared = 0;

	//table precomputing
	partial_primorial_precompute();
	crt_coefficient_precompute();
	be_pair_precompute();
	positive_be_quad_precompute();
	negative_be_quad_precompute();

	//free memory allocated for precomputed tables
	free_table_memory();
}

void free_table_memory()
{
	free(partial_primorial_table);
	free(crt_coefficient_table);
	free(be_pair_table);
	free(positive_be_quad_table);
	free(negative_be_quad_table);
	free(congruency_class_size_table);
	free(positive_be_quad_prime_size_table);
	free(negative_be_quad_prime_size_table);
}

//calculates the modulus of A, B, C in the solution
void partial_primorial_precompute()
{
	partial_primorial_table = (unsigned long long *) calloc((1 << num_primes),
															sizeof(unsigned long long));

	partial_primorial_table[0] = 1;
	partial_primorial_table[1] = 2;

	unsigned int i, j,
			 bit_mask;

	for (i = 1; i <  num_primes; i++)
	{
		bit_mask = (1 << i) - 1;

		for (j = (1 << i); j < (1 << (i + 1)); j++)
		{
			partial_primorial_table[j] = primes[i] * partial_primorial_table[j & bit_mask];

			#ifdef PRINT_DEBUG
				printf("%I64d\n", partial_primorial_table[j]);
			#endif
		}
	}
}

void crt_coefficient_precompute()
{
	crt_coefficient_table = (unsigned long long *) calloc((1 << (num_primes +
																 power_two_num_primes)),
														  sizeof(unsigned long long));

	unsigned int i, j;
	unsigned long long coefficient;

	for (i = 0; i < (1 << num_primes); i++)
	{
		for(j = 0; j < num_primes; j++)
		{
			if ((i >> j) & 1)
			{
				coefficient = ext_euclidean_coefficient((partial_primorial_table[i] / primes[j]),
													   ((unsigned long long) primes[j]));
				coefficient *= (partial_primorial_table[i] / primes[j]);

				ulonglong_3D_array_set(coefficient,
									   i, j, 0,
									   crt_coefficient_table,
									   power_two_num_primes,
									   0);

				#ifdef PRINT_DEBUG
					printf("Coefficient for bitstring %u and prime %u is %I64d\n", i, primes[j], coefficient);
				#endif
			}
		}
	}
}

void be_pair_precompute()
{
	unsigned int n, base, exponent,
				 congruency_class, congruency_class_size;
	be_pair cur_pair;

	congruency_class_size_table = (unsigned int *) calloc((1 << (power_two_num_primes +
																 power_two_biggest_prime)),
														  sizeof(unsigned int));

	be_pair_table = (be_pair *) calloc((1 << (power_two_num_primes +
											  power_two_biggest_prime +
											  power_two_biggest_prime_squared)),
									   sizeof(be_pair));

	for (n = 0; n < num_primes; n++)
	{
		for(base = 1; base < primes[n]; base++)
		{
			for(exponent = 0; exponent < (primes[n] - 1); exponent++)
			{
				congruency_class = modular_exponentiation(base, exponent, primes[n]);

				cur_pair.base = base;
				cur_pair.exponent = exponent;
				cur_pair.congruence = congruency_class;
				cur_pair.modulus = primes[n];

				congruency_class_size = uint_3D_array_get(n, cur_pair.congruence, 0,
														  congruency_class_size_table,
														  power_two_biggest_prime,
														  0);

				be_pair_3D_array_set(cur_pair,
									 n, cur_pair.congruence, congruency_class_size,
									 be_pair_table,
									 power_two_biggest_prime,
									 power_two_biggest_prime_squared);

				uint_3D_array_set((congruency_class_size + 1),
								  n, cur_pair.congruence, 0,
								  congruency_class_size_table,
								  power_two_biggest_prime,
								  0);

				if (++congruency_class_size > max_congruency_class_size)
				{
					max_congruency_class_size = congruency_class_size;
				}

				#ifdef PRINT_DEBUG
					printf("(%u, %u, %u) = %u\n", base, exponent, primes[n], congruency_class);
				#endif
			}
		}
	}
}

void positive_be_quad_precompute()
{
	unsigned int n, congruency_class, num_first_be_pair, num_second_be_pair,
				 prime_be_quad_matches, congruency_class_size;
	be_quad cur_quad;

	power_two_max_congruency_class_size_squared = 2 * power_two_exp(max_congruency_class_size);

	positive_be_quad_table = (be_quad *) calloc((1 << (power_two_num_primes +
											  	  	   power_two_biggest_prime +
											  	  	   power_two_max_congruency_class_size_squared)),
											  	sizeof(be_quad));

	positive_be_quad_prime_size_table = (unsigned int *) calloc(num_primes,
												  	  	   	    sizeof(unsigned int));

	for (n = 0; n < num_primes; n++)
	{
		prime_be_quad_matches = 0;

		for(congruency_class = 1; congruency_class < primes[n]; congruency_class++)
		{
			congruency_class_size = uint_3D_array_get(n, congruency_class, 0,
													  congruency_class_size_table,
													  power_two_biggest_prime,
													  0);

			for (num_first_be_pair = 0; num_first_be_pair < congruency_class_size; num_first_be_pair++)
			{
				for (num_second_be_pair = 0; num_second_be_pair < congruency_class_size; num_second_be_pair++)
				{
					cur_quad.first_be_pair = be_pair_3D_array_get(n, congruency_class, num_first_be_pair,
															   	  be_pair_table,
															   	  power_two_biggest_prime,
															   	  power_two_biggest_prime_squared);

					cur_quad.second_be_pair = be_pair_3D_array_get(n, congruency_class, num_second_be_pair,
															   	   be_pair_table,
															   	   power_two_biggest_prime,
															   	   power_two_biggest_prime_squared);

					be_quad_3D_array_set(cur_quad,
										 n, prime_be_quad_matches++, 0,
										 positive_be_quad_table,
										 (power_two_biggest_prime + power_two_max_congruency_class_size_squared),
										 0);

					#ifdef PRINT_DEBUG
						printf("(%u, %u, %u) = (%u, %u, %u) = %u\n", cur_quad.first_be_pair.base,
																	 cur_quad.first_be_pair.exponent,
																	 cur_quad.first_be_pair.modulus,
																	 cur_quad.second_be_pair.base,
																	 cur_quad.second_be_pair.exponent,
																	 cur_quad.second_be_pair.modulus,
																	 cur_quad.first_be_pair.congruence);
					#endif
				}
			}

			positive_be_quad_prime_size_table[n] = prime_be_quad_matches;
		}

		if (n % 2)
		{
			num_positive_be_quad_combos_high *= prime_be_quad_matches;
		}
		else
		{
			num_positive_be_quad_combos_low *= prime_be_quad_matches;
		}

		#ifdef PRINT_DEBUG
			printf("Number of positive base-power pair matches for %u: %u\n", primes[n], prime_be_quad_matches);
		#endif
	}

	#ifdef PRINT_DEBUG
	printf("Total number of positive base-power pair combos for %u primes: %I64d * %I64d\n",
		   num_primes, num_positive_be_quad_combos_low, num_positive_be_quad_combos_high);
	#endif
}

void negative_be_quad_precompute()
{
	unsigned int n, congruency_class, num_first_be_pair, num_second_be_pair,
				 prime_be_quad_matches, first_congruency_class_size, second_congruency_class_size;
	be_quad cur_quad;

	power_two_max_congruency_class_size_squared = 2 * power_two_exp(max_congruency_class_size);

	negative_be_quad_table = (be_quad *) calloc((1 << (power_two_num_primes +
													   power_two_biggest_prime +
													   power_two_max_congruency_class_size_squared)),
												sizeof(be_quad));

	negative_be_quad_prime_size_table = (unsigned int *) calloc(num_primes,
													  	   	    sizeof(unsigned int));

	for (n = 0; n < num_primes; n++)
	{
		prime_be_quad_matches = 0;

		for(congruency_class = 1; congruency_class < primes[n]; congruency_class++)
		{
			first_congruency_class_size = uint_3D_array_get(n, congruency_class, 0,
															congruency_class_size_table,
															power_two_biggest_prime,
															0);

			second_congruency_class_size = uint_3D_array_get(n, (primes[n] - congruency_class), 0,
															 congruency_class_size_table,
															 power_two_biggest_prime,
															 0);

			for (num_first_be_pair = 0; num_first_be_pair < first_congruency_class_size; num_first_be_pair++)
			{
				for (num_second_be_pair = 0; num_second_be_pair < second_congruency_class_size; num_second_be_pair++)
				{
					cur_quad.first_be_pair = be_pair_3D_array_get(n, congruency_class, num_first_be_pair,
															   	  be_pair_table,
															   	  power_two_biggest_prime,
															   	  power_two_biggest_prime_squared);

					cur_quad.second_be_pair = be_pair_3D_array_get(n, (primes[n] - congruency_class), num_second_be_pair,
															   	   be_pair_table,
															   	   power_two_biggest_prime,
															   	   power_two_biggest_prime_squared);

					be_quad_3D_array_set(cur_quad,
										 n, prime_be_quad_matches++, 0,
										 negative_be_quad_table,
										 (power_two_biggest_prime + power_two_max_congruency_class_size_squared),
										 0);

					#ifdef PRINT_DEBUG
						printf("(%u, %u, %u) = -(%u, %u, %u) = %u\n", cur_quad.first_be_pair.base,
																  	  cur_quad.first_be_pair.exponent,
																  	  cur_quad.first_be_pair.modulus,
																  	  cur_quad.second_be_pair.base,
																  	  cur_quad.second_be_pair.exponent,
																  	  cur_quad.second_be_pair.modulus,
																  	  cur_quad.first_be_pair.congruence);
					#endif
				}
			}

			negative_be_quad_prime_size_table[n] = prime_be_quad_matches;
		}

		if (n % 2)
		{
			num_negative_be_quad_combos_high *= prime_be_quad_matches;
		}
		else
		{
			num_negative_be_quad_combos_low *= prime_be_quad_matches;
		}

		#ifdef PRINT_DEBUG
			printf("Number of negative base-power pair matches for %u: %u\n", primes[n], prime_be_quad_matches);
		#endif
	}

	#ifdef PRINT_DEBUG
		printf("Total number of negative base-power pair combos for %u primes: %I64d * %I64d\n",
			   num_primes, num_negative_be_quad_combos_low, num_negative_be_quad_combos_high);
	#endif
}

inline int moduli_set_solver(be_quad *congruent_quad_combo,
			   	   	  	  	 unsigned int bit_mask_top, unsigned int bit_mask_bottom,
			   	   	  	  	 full_crt_solution *solution)
{
	unsigned int n, variable_code, primorial_mask;
	be_quad cur_quad;
	unsigned long long	a_congruence = 0, b_congruence = 0, c_congruence = 0,
						*crt_coefficients;
	partial_crt_solution x_solution, y_solution, z_solution;

	primorial_mask = bit_mask_top | bit_mask_bottom;

	//coefficient array for this prime mask
	crt_coefficients = &(ulonglong_3D_array_get(primorial_mask, 0, 0,
											  	crt_coefficient_table,
											  	power_two_num_primes, 0));

	x_solution.congruence = 0;
	x_solution.modulus = 1;
	y_solution.congruence = 0;
	y_solution.modulus = 1;
	z_solution.congruence = 0;
	z_solution.modulus = 1;

	for (n = 0; n < num_primes; n++)
	{
		variable_code = ((bit_mask_top & 1) << 1) + (bit_mask_bottom & 1);
		cur_quad = congruent_quad_combo[n];

		switch(variable_code)
		{
			case 1:
				b_congruence += cur_quad[n].first_be_pair.base * crt_coefficients[n];
				c_congruence += cur_quad[n].second_be_pair.base * crt_coefficients[n];

				if ((succ_sub_solver(&y_solution, cur_quad.first_be_pair.exponent, (cur_quad.first_be_pair.modulus - 1)) < 0) ||
					(succ_sub_solver(&z_solution, cur_quad.second_be_pair.exponent, (cur_quad.second_be_pair.modulus - 1)) < 0))
				{
					return -1;
				}

				break;
			case 2:
				a_congruence += congruent_quad_combo[n].first_be_pair.base * crt_coefficients[n];
				c_congruence += congruent_quad_combo[n].second_be_pair.base * crt_coefficients[n];

				if ((succ_sub_solver(&x_solution, cur_quad.first_be_pair.exponent, (cur_quad.first_be_pair.modulus - 1)) < 0) ||
					(succ_sub_solver(&z_solution, cur_quad.second_be_pair.exponent, (cur_quad.second_be_pair.modulus - 1)) < 0))
				{
					return -1;
				}

				break;
			case 3:
				a_congruence += congruent_quad_combo[n].first_be_pair.base * crt_coefficients[n];
				b_congruence += congruent_quad_combo[n].second_be_pair.base * crt_coefficients[n];

				if ((succ_sub_solver(&x_solution, cur_quad.first_be_pair.exponent, (cur_quad.first_be_pair.modulus - 1)) < 0) ||
					(succ_sub_solver(&y_solution, cur_quad.second_be_pair.exponent, (cur_quad.second_be_pair.modulus - 1)) < 0))
				{
					return -1;
				}

				break;
			//cover the case of zero as well
			default :
				break;
		}

		bit_mask_top = (bit_mask_top >> 1);
		bit_mask_bottom = (bit_mask_bottom >> 1);
	}

	solution->a_congruence = a_congruence;
	solution->b_congruence = b_congruence;
	solution->c_congruence = c_congruence;
	solution->x_congruence = x_solution.congruence;
	solution->y_congruence = y_solution.congruence;
	solution->z_congruence = z_solution.congruence;
	solution->abc_modulus = partial_primorial_table[primorial_mask];
	solution->x_modulus = x_solution.modulus;
	solution->y_modulus = y_solution.modulus;
	solution->z_modulus = z_solution.modulus;

	return 1;
}


inline int succ_sub_solver(partial_crt_solution *solution_m,
						   unsigned long long solution_a_congruence,
						   unsigned long long solution_a_modulus)
{
	unsigned long long gcd, multip_inverse;

	if (gcd_of_three(solution_m->modulus,
					 solution_a_congruence - solution_m->congruence,
					 solution_a_modulus,
					 &gcd) < 0)
	{
		return -1;
	}
	else
	{
		if (ge_multip_inverse(((solution_m-> modulus) / gcd),
							  (solution_a_modulus / gcd),
							  &multip_inverse) < 0)
		{
			return -1;
		}
		else
		{
			solution_m->congruence = (multip_inverse * (((solution_a_congruence - solution_m->congruence) / gcd) *
									  	  	  	  	  	solution_m->modulus)) +
									 solution_m->congruence;
			solution_m->modulus = solution_m->modulus * (solution_a_modulus / gcd);

			return 1;
		}
	}
}

inline int gcd_of_three(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long *gcd)
{
	return 0;
}

inline int get_multip_inverse(unsigned long long congruence, unsigned long long moduli, unsigned long long *multip_inverse)
{
	return 0;
}

unsigned int power_two_exp(unsigned int number)
{
	return (((unsigned int) (log(number) / log(2))) + 1);
}

inline unsigned int modular_exponentiation(unsigned int base, unsigned int exponent, unsigned int modulus)
{
	unsigned int congruency_class = 1;

	while (exponent > 0)
	{
		if ((exponent & 1) == 1)
		{
			congruency_class = (congruency_class * base) % modulus;
		}

		exponent = exponent >> 1;
		base = (base * base) % modulus;
	}

	return congruency_class;
}

inline unsigned long long ext_euclidean_coefficient(unsigned long long x, unsigned long long y)
{
	unsigned long long a_old = 1, b_old = 0, remainder_old = x,
					   a_new = 0, b_new = 1, remainder_new = y, multiplier,
					   a_temp, b_temp, remainder_temp;
	unsigned int x_bigger = 1;

	if (x < y)
	{
		x_bigger = 0;

		remainder_old = y;
		remainder_new = x;
	}

	while (remainder_new != 0)
	{
		multiplier = remainder_old / remainder_new;

		a_temp = a_new;
		b_temp = b_new;
		remainder_temp = remainder_new;

		a_new = a_old - (multiplier * a_new);
		b_new = b_old - (multiplier * b_new);
		remainder_new = remainder_old % remainder_new;

		a_old = a_temp;
		b_old = b_temp;
		remainder_old = remainder_temp;
	}

	if (x_bigger)
	{
		return a_old;
	}
	else
	{
		return b_old;
	}
}

inline void uint_3D_array_set(unsigned int value,
						  	  unsigned int x, unsigned int y, unsigned int z,
						  	  unsigned int *table,
						  	  unsigned int table_dimension_y,
						  	  unsigned int table_dimension_z)
{
	table[(((x << table_dimension_y) + y) << table_dimension_z) + z] = value;
}

inline unsigned int uint_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									  unsigned int *table,
									  unsigned int table_dimension_y,
									  unsigned int table_dimension_z)
{
	return table[(((x << table_dimension_y) + y) << table_dimension_z) + z];
}

inline void ulonglong_3D_array_set(unsigned long long value,
						  	  	   unsigned int x, unsigned int y, unsigned int z,
						  	  	   unsigned long long *table,
						  	  	   unsigned int table_dimension_y,
						  	  	   unsigned int table_dimension_z)
{
	table[(((x << table_dimension_y) + y) << table_dimension_z) + z] = value;
}

inline unsigned long long ulonglong_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									  	  	     unsigned long long *table,
									  	  	     unsigned int table_dimension_y,
									  	  	     unsigned int table_dimension_z)
{
	return table[(((x << table_dimension_y) + y) << table_dimension_z) + z];
}

inline void be_pair_3D_array_set(be_pair value,
						  	  	 unsigned int x, unsigned int y, unsigned int z,
						  	  	 be_pair *table,
						  	  	 unsigned int table_dimension_y,
						  	  	 unsigned int table_dimension_z)
{
	table[(((x << table_dimension_y) + y) << table_dimension_z) + z] = value;
}

inline be_pair be_pair_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									be_pair *table,
									unsigned int table_dimension_y,
									unsigned int table_dimension_z)
{
	return table[(((x << table_dimension_y) + y) << table_dimension_z) + z];
}

inline void be_quad_3D_array_set(be_quad value,
						  	  	 unsigned int x, unsigned int y, unsigned int z,
						  	  	 be_quad *table,
						  	  	 unsigned int table_dimension_y,
						  	  	 unsigned int table_dimension_z)
{
	table[(((x << table_dimension_y) + y) << table_dimension_z) + z] = value;
}

inline be_quad be_quad_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									be_quad *table,
									unsigned int table_dimension_y,
									unsigned int table_dimension_z)
{
	return table[(((x << table_dimension_y) + y) << table_dimension_z) + z];
}

/*
void timer_start()
{
	QueryPerformanceCounter(&start);
}

void timer_finish(char *segment_name)
{
	QueryPerformanceCounter(&finish);

	LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);

    printf("\n##############################################\n");
	printf("%s took %8.8f microseconds", segment_name,
		   (((float) (finish.QuadPart - start.QuadPart)) / frequency.QuadPart));
	printf("\n##############################################\n");
}
*/
