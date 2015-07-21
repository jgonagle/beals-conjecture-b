/*
 * BealsConjectureII.h
 *
 *  Created on: Jul 30, 2011
 *      Author: John McGonagle
 */

typedef struct pair_be {
	unsigned int base;
	unsigned int exponent;
	unsigned int congruence;
	unsigned int modulus;
} be_pair;

typedef struct quad_be {
	be_pair first_be_pair;
	be_pair second_be_pair;
} be_quad;

typedef struct solution_crt_partial {
	unsigned long long congruence;
	unsigned long long modulus;
} partial_crt_solution;

typedef struct solution_crt_full {
	unsigned long long a_congruence;
	unsigned long long b_congruence;
	unsigned long long c_congruence;
	unsigned long long x_congruence;
	unsigned long long y_congruence;
	unsigned long long z_congruence;
	unsigned long long abc_modulus;
	unsigned long long x_modulus;
	unsigned long long y_modulus;
	unsigned long long z_modulus;
} full_crt_solution;

//primes and exponents related to storage size of tables
unsigned int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
unsigned int num_primes,
			 max_congruency_class_size,
			 power_two_num_primes,
			 power_two_biggest_prime,
			 power_two_biggest_prime_squared,
			 power_two_max_congruency_class_size_squared;

unsigned long long num_positive_be_quad_combos_low,
				   num_positive_be_quad_combos_high,
				   num_negative_be_quad_combos_low,
				   num_negative_be_quad_combos_high;

//tables used for precomputing
unsigned long long *partial_primorial_table;
unsigned long long *crt_coefficient_table;
be_pair *be_pair_table;
be_quad *positive_be_quad_table;
be_quad *negative_be_quad_table;
unsigned int *congruency_class_size_table;
unsigned int *positive_be_quad_prime_size_table;
unsigned int *negative_be_quad_prime_size_table;

//for timer
LARGE_INTEGER start, finish;

void timer_start();
void timer_finish();

//does all precomputing and deallocation of dynamic memory
void do_precompute();

//handles dynamic memory allocation for precomputing tables
void free_table_memory();

//functions to separate precomputing
void partial_primorial_precompute();
void crt_coefficient_precompute();
void be_pair_precompute();
void positive_be_quad_precompute();
void negative_be_quad_precompute();

//functions to solve moduli sets
inline int moduli_set_solver(be_quad *congruent_quad_combo,
			   	   	  	  	 unsigned int bit_mask_top, unsigned int bit_mask_bottom,
			   	   	  	  	 full_crt_solution *solution);
inline int succ_sub_solver(partial_crt_solution *compound_solution,
						   unsigned long long simple_solution_congruence,
						   unsigned long long simple_solution_modulus);
inline int gcd_of_three(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long *gcd);
inline int get_multip_inverse(unsigned long long congruence, unsigned long long moduli, unsigned long long *multip_inverse);

//helper functions
unsigned int power_two_exp(unsigned int number);
inline unsigned int modular_exponentiation(unsigned int base,
										   unsigned int exponent,
										   unsigned int modulus);
inline unsigned long long ext_euclidean_coefficient(unsigned long long x,
													unsigned long long y);

//set and get for multidimensional arrays
inline void uint_3D_array_set(unsigned int value,
							  unsigned int x, unsigned int y, unsigned int z,
							  unsigned int *table,
							  unsigned int table_dimension_y,
							  unsigned int table_dimension_z);

inline unsigned int uint_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									  unsigned int *table,
									  unsigned int table_dimension_y,
									  unsigned int table_dimension_z);

inline void ulonglong_3D_array_set(unsigned long long value,
						  	  	   unsigned int x, unsigned int y, unsigned int z,
						  	  	   unsigned long long *table,
						  	  	   unsigned int table_dimension_y,
						  	  	   unsigned int table_dimension_z);

inline unsigned long long ulonglong_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									  	  	     unsigned long long *table,
									  	  	     unsigned int table_dimension_y,
									  	  	     unsigned int table_dimension_z);

inline void be_pair_3D_array_set(be_pair value,
						  	  	 unsigned int x, unsigned int y, unsigned int z,
						  	  	 be_pair *table,
						  	  	 unsigned int table_dimension_y,
						  	  	 unsigned int table_dimension_z);

inline be_pair be_pair_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									be_pair *table,
									unsigned int table_dimension_y,
									unsigned int table_dimension_z);

inline void be_quad_3D_array_set(be_quad value,
						  	  	 unsigned int x, unsigned int y, unsigned int z,
						  	  	 be_quad *table,
						  	  	 unsigned int table_dimension_y,
						  	  	 unsigned int table_dimension_z);

inline be_quad be_quad_3D_array_get(unsigned int x, unsigned int y, unsigned int z,
									be_quad *table,
									unsigned int table_dimension_y,
									unsigned int table_dimension_z);
