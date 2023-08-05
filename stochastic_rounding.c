#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <stdint.h>
#include <stdbool.h>


union DoubleToInt {
  double dVal;
  uint64_t iVal;
};

/*
  x function that rounds a binary64 value to a binary32 value
  stochastically. Implemented by treating FP number representations
  as integer values.
*/
float SR(double x) {
  union DoubleToInt temp;
  temp.dVal = x;
  uint32_t r = rand() & 0x1FFFFFFF;
  temp.iVal += r;
  temp.iVal = temp.iVal & 0xFFFFFFFFE0000000;

  return (float)temp.dVal;
}

const bool plotting = true;
const bool test_mode = false;
const bool verbose = false;

/* --------------------------------- */
/*              PART 1               */
/* --------------------------------- */
// Generate a random probability drawn from a uniform distribution over the interval [0, 1)
float random_generator(int v){
	if(v == 1){
		return arc4random_uniform((double)RAND_MAX)/((float)RAND_MAX);
	}
	else{
		return (float) rand() / ((double) (RAND_MAX) + 1);
	}
}

/*
 * Calculates the normalized distance between a given value "x" and a line segment
 * defined by two endpoint values "val1" and "val2". Returned value is a single-precision floating-point number.
 */
float normalized_distance(double x, float val1, float val2){
	return (float) fabs(x - val1) / fabs(val1 - val2);
}

// Implement RZ(x) - rounding toward zero
float round_to_zero(double x) {
	return nextafterf((float) x, 0);
}

// Implement RA(x) - rounding away from zero
float round_away_from_zero(double x) {
	float x_rounded = (float) x;
	if(x > 0){
		if(x_rounded > x){
			return x_rounded;
		}
		return nextafterf(x_rounded, INFINITY);
	}
	if(x_rounded > x){
		return nextafterf(x_rounded, -INFINITY);
	}
	return x_rounded;
}

/*
 * Implement SR_alternative according to the Eqn 1.
 * "Stochastic rounding avoids stagnation and that the computed result has expected value equal to the exact sum."
 * - Stochastic Rounding and Its Probabilistic Backward Error Analysis by
 * Connolly, Michael P. and Higham, Nicholas J. and Mary, Theo
*/
float SR_alternative(double x) {
	float RZ = round_to_zero(x);
	float RA = round_away_from_zero(x);
	float probability = random_generator(1);
	float norm_distance = normalized_distance(x, RZ, RA);
	if(verbose) {
		printf("RZ(x):              %.60f \n", RZ);
		printf("RA(x):              %.60f \n", RA);
		printf("P:                  %.60f \n", probability);
		printf("p:                  %.60f \n", norm_distance);
	}
	if(probability < norm_distance){
		return RA;
	}
	return RZ;
}

void plot_absolute_error_part_1(double *provided_error, double *alt_error, int n, int step) {
	char * commandsForGnuplot[] = {
			"set title \"Absolute Error in Stochastic Rounding\"",
			"set ylabel \"Absolute Rounding Error\"",
			"set xlabel \"Number of Roundings\"",
			"plot 'part1.temp' using 1:2 with lines dt 2 lt 9 lc 7 title 'Absolute Error - Provided SR', 'part1.temp' using 1:3 with lines lc 6 title 'Absolute Error - Alternative SR'"
	};
	FILE * temp = fopen("part1.temp", "w");
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	for(int i = 0; i < n; i++){
		fprintf(temp, "%d %.60f %.60f\n", i * step, provided_error[i], alt_error[i]);
	}
	for (int i=0; i < 4; i++){
		fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
	}
}

/*
 * This method computes the sum of two double-precision floating-point numbers `a` and `b` using
 * the fast 2-sum algorithm, and stores the result and error in the variables pointed to by s
 * and err, respectively. The result is stored as a single-precision floating-point number.
*/
void fast2sum(double a, double b, float *s, float *err){
	*s = a + b;
	*err = b - (*s - a);
}

void plot_absolute_error_part_2(double *provided_error, double *alt_error, double *comp_error, int n, int step) {
	char * commandsForGnuplot[] = {
			"set title \"Absolute Error in Approximating Taylor Series\"",
			"set ylabel \"Absolute Rounding Error\"",
			"set xlabel \"Number of Roundings\"",
			"set logscale x 10",
			"set logscale y 10",
			"plot 'part2.temp' using 1:2 with lines dt 2 lt 9 lc 7 title 'binary32 RN', 'part2.temp' using 1:3 with lines title 'binary32 SR' lt rgb 'blue', 'part2.temp' using 1:4 with lines dt 2 lt 2 title 'binary32 COMP'"
	};
	FILE * temp = fopen("part2.temp", "w");
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	for(int i = 0; i < n; i++){
		fprintf(temp, "%d %.60f %.60f %.60f\n", i * step, provided_error[i], alt_error[i], comp_error[i]);
	}
	for (int i=0; i <= 5; i++){
		fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
	}
}

double mean_average(double *arr, int size){
	double mean = 0.0;
	for (int i = 0; i < size; i++) {
		mean += arr[i];
	}
	return mean / size;
}

/*
 * This method checks whether two floating-point values are almost equal.
 * It returns true if the absolute difference between x and y is less than or equal to
 * the maximum relative difference allowed, which is FLT_EPSILON times the larger of the
 * absolute values of x and y. Otherwise, it returns false.
*/
bool are_almost_equal(float x, float y){
	float maxRelDiff = FLT_EPSILON;  // machine precision
	float diff = fabs(x - y);
	x = fabs(x);
	y = fabs(y);
	float largest;
	if(y > x){
		largest = y;
	}
	else{
		largest = x;
	}
	/*
	 * Use a relative instead of an absolute tolerance to make the comparison less sensitive to `x` and `y`'s magnitudes
	 */
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}

double random_generator_double(){
	return (double) arc4random_uniform(RAND_MAX)/((double)RAND_MAX);
}

void plot_absolute_error_part_3(double *provided_error, double *alt_error, int n, int step) {
	char * commandsForGnuplot[] = {
			"set title \"Absolute Error in calculating Vector Inner Product\"",
			"set ylabel \"Absolute Rounding Error\"",
			"set xlabel \"Number of iterations\"",
			"set logscale y 10",
			"plot 'part3.temp' using 1:2 with lines dt 2 lt 9 lc 7 title 'binary32 RN', 'part3.temp' using 1:3 with lines lc 6 title 'binary32 SR'"
	};
	FILE * temp = fopen("part3.temp", "w");
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	for(int i = 0; i < n; i++){
		fprintf(temp, "%d %.60f %.60f \n", i * step, provided_error[i], alt_error[i]);
	}
	for (int i=0; i <= 4; i++){
		fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
	}
}

const long int K = 5000000;

int main() {
	// Seed random generator
	if(test_mode){
		srand(1);
	}
	else {
		srand(time(NULL));
	}

	// An arbitrary value for alternative_tmp.
	double sample = M_PI;
	double avg = 0;
	double avg_alternative = 0;

	// Calculate the neighbouring binary32 values.
	float closest = (float)sample;
	float down, up;
	if (closest > sample) {
		up = closest;
		down = nextafterf(closest, -INFINITY);
	}
	else {
		up = nextafterf(closest, INFINITY);
		down = closest;
	}

	// Round many times, and calculate the average values as well as count
	// the numbers of times alternative_tmp was up/down.
	int count_up_avg = 0;
	int count_up_avg_alt = 0;
	double avg_tmp;
	double alternative_tmp;

	int max_step = 1000;
	int step = 0;
	int error_index = 0;
	int n = K / max_step;
	double absolute_error[n];
	double absolute_error_alt[n];
	for (int i = 1; i <= K; i++) {
		// Counting of provided SR method
		avg_tmp = SR(sample);
		if(avg_tmp > sample){
			count_up_avg += 1;
		}
		avg += avg_tmp;

		// Counting of implemented SR method
		alternative_tmp = SR_alternative(sample);
		if(alternative_tmp > sample){
			count_up_avg_alt += 1;
		}
		avg_alternative += alternative_tmp;

		if(step == max_step){
			error_index = (int)i / max_step;
			absolute_error[error_index] = fabs((double) avg / i - sample);
			absolute_error_alt[error_index] = fabs((double) avg_alternative / i - sample);
			step = 0;
		}
		step += 1;
	}
	avg /= K;
	avg_alternative /= K;

	printf("------------------- PART I -------------------\n");
	printf("Value being rounded:           %.60f \n", sample);
	printf("SR average value:              %.60f \n", avg);
	printf("SR_alternative average value:  %.60f \n", avg_alternative);
	printf("Binary32 value before:         %.60f \n", down);
	printf("Binary32 value after:          %.60f \n", up);
	printf("Closest binary32:              %.60f \n", closest);

	// Print out the average of all rounded values
	if(plotting) {
		plot_absolute_error_part_1(absolute_error, absolute_error_alt, n, max_step);
	}

	// Check that SR_alternative function is correct by comparing the probabilities
	// of alternative_tmp up/down, and the expected probability. Print them out
	// below.
	float prob_RU, prob_RD;
	prob_RU = (float) count_up_avg_alt/K;
	prob_RD = 1 - prob_RU;
	printf("SR Alternative Approximate Probability RU:   %.6f \n", prob_RU);
	printf("SR Alternative Approximate Probability RD:   %.6f \n", prob_RD);

	prob_RD = normalized_distance(sample, up, down);
	prob_RU = 1 - prob_RD;
	printf("SR Expected Probability RU:          %.6f \n", prob_RU);
	printf("SR Expected Probability RD:          %.6f \n", prob_RD);

	prob_RU = (float) count_up_avg/K;
	prob_RD = 1 - prob_RU;
	printf("SR Provided Probability RU:          %.6f \n", prob_RU);
	printf("SR Provided Probability RD:          %.6f \n", prob_RD);

	/* --------------------------------- */
	/*              PART 2               */
	/* --------------------------------- */
	printf("------------------- PART II -------------------\n");
	long int N = 500000000;
	float fharmonic = 0;
	float fharmonic_next = 0;
	float fharmonic_sr = 0;
	float fharmonic_comp = 0;
	double dharmonic = 0;

	// Error term in the compensated summation.
	float t = 0;
	max_step = 1000000;
	step = 0;

	n = N / max_step;
	bool stagnation = false;
	double abs_error_32[n];
	double abs_error_32_alt[n];
	double abs_error_32_comp[n];
	for (int i = 1; i <= N; i++) {
		// Recursive sum, binary32 RN
		fharmonic += (float)1/(float)i;
		fharmonic_next = fharmonic + (float)1/(float)(i+1);
		// Recursive sum, binary32 SR Alternative
		fharmonic_sr = SR_alternative((double) fharmonic_sr + (double)1/(double)i);
		// Compensated sum, binary32 RN
		float addend = (float)1/(float)i + t;
		fast2sum(fharmonic_comp, addend, &fharmonic_comp, &t);
		// Recursive sum, binary64 RN
		dharmonic += (double)1/(double)i;
		if(!stagnation && are_almost_equal(fharmonic, fharmonic_next)){
			printf("Sum is stagnating after %d iterations\n", i);
			stagnation = true;
		}
		if(step == max_step){
			error_index = (int)i / max_step;
			abs_error_32[error_index] = fabs((double) fharmonic - dharmonic);
			abs_error_32_alt[error_index] = fabs((double) fharmonic_sr - dharmonic);
			abs_error_32_comp[error_index] = fabs((double) fharmonic_comp - dharmonic);
			step = 0;
		}
		step += 1;
	}
	double abs_error_32_final = fabs((double) fharmonic - dharmonic);
	double abs_error_32_alt_final = fabs((double) fharmonic_sr - dharmonic);
	double abs_error_32_comp_final = fabs((double) fharmonic_comp - dharmonic);

	printf("Values of the harmonic series after %ld iterations \n", N);
	printf("Recursive summation, binary32:          %.30f \n", fharmonic);
	printf("Recursive summation with SR, binary32:  %.30f \n", fharmonic_sr);
	printf("Compensated summation, binary32:        %.30f \n", fharmonic_comp);
	printf("Recursive summation, binary64:          %.30f \n", dharmonic);

	printf("Absolute errors of the harmonic series  \n");
	printf("Recursive summation error, binary32:          %.30f \n", abs_error_32_final);
	printf("Recursive summation with SR error, binary32:  %.30f \n", abs_error_32_alt_final);
	printf("Compensated summation error, binary32:        %.30f \n", abs_error_32_comp_final);

	if(abs_error_32_final > abs_error_32_alt_final && abs_error_32_final > abs_error_32_comp_final){
		printf("Sorted absolute errors of the harmonic series summation:\n");
		printf("Recursive summation error, binary32 has the highest error\n");
		if(abs_error_32_comp_final < abs_error_32_final && abs_error_32_comp_final < abs_error_32_alt_final){
			printf("Compensated summation error, binary32 has the lowest error\n");
		}
	}

	// Print out the absolute errors with respect to the double precision Taylor series
	if(plotting){
		plot_absolute_error_part_2(abs_error_32, abs_error_32_alt, abs_error_32_comp, n, max_step);
	}

	/* --------------------------------- */
	/*              PART 3               */
	/* --------------------------------- */
	printf("------------------- PART III -------------------\n");
	float fvec_res = 0;
	float fvec_res_sr = 0;
	double dvec_res = 0;
	int VEC_SIZE = K;
	// Calculate vector inner product directly instead of storing large sequences of memory
	double vec1_tmp, vec2_tmp;
	max_step = 10000;
	step = 0;
	n = VEC_SIZE / max_step;
	double abs_error_32_vec[n];
	double abs_error_32_alt_vec[n];

	// Simulate calculating vector Inner Product with different Rounding modes
	for (int i = 0; i < VEC_SIZE; i++) {
		// Generate two vector elements and simulate dot product
		vec1_tmp = random_generator_double();
		vec2_tmp = random_generator_double();
		fvec_res += (float) (vec1_tmp * vec2_tmp);
		fvec_res_sr = SR_alternative((double) fvec_res_sr + (double) (vec1_tmp * vec2_tmp));
		dvec_res += (double) (vec1_tmp * vec2_tmp);
		if(step == max_step){
			error_index = (int)i / max_step;
			abs_error_32_vec[error_index] = fabs((double) fvec_res - dvec_res);
			abs_error_32_alt_vec[error_index] = fabs((double) fvec_res_sr - dvec_res);
			step = 0;
		}
		step += 1;
	}
	printf("Absolute error in vector inner product calculation:\n");
	printf("in binary32 with RN:  %.5f \n", fabs((double) fvec_res - dvec_res));
	printf("in binary32 with SR:  %.5f \n", fabs((double) fvec_res_sr - dvec_res));

	if(plotting){
		plot_absolute_error_part_3(abs_error_32_vec, abs_error_32_alt_vec, n, max_step);
	}
	return 0;
}
