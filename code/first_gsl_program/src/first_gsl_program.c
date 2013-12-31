/*
 ============================================================================
 Name        : first_gsl_program.c
 Author      : Mark Greenaway
 Version     :
 Copyright   : Your copyright notice
 Description : Generate one hundred normally distirbuted numbers, then show a histogram
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

const int SAMPLES = 1e5;
const int COLUMNS = 80;

void generate_samples(gsl_rng* rnd, double *norm)
{
	int i;
	/* Generate some samples from N(100, 1) */
	for (i = 0; i < SAMPLES; i++) {
		norm[i] = 100.0 + gsl_ran_gaussian(rnd, 1.0);
	}
}

void print_samples(const double norm[])
{
	int i;

	for (i = 0; i < SAMPLES; i++) {
		if (i % 5 == 0 && i > 0) {
			printf("\n");
		}
		printf("%d: %f\t", i, norm[i]);
	}
	printf("\n");
}

int main(void) {
	int i;
	double norm[SAMPLES];

	gsl_rng *rnd = gsl_rng_alloc(gsl_rng_taus);
	gsl_histogram *norm_hist = gsl_histogram_alloc(COLUMNS);

	/* Generate some samples from N(100, 1) */
	generate_samples(rnd, norm);
	print_samples(norm);

	/* Construct a histogram of the samples */
	gsl_histogram_set_ranges_uniform(norm_hist, 98.0, 102.0);

	for (i = 0; i < SAMPLES; i++) {
		gsl_histogram_increment(norm_hist, norm[i]);
	}

	/* Display the histogram */
	for (i = 0; i < COLUMNS; i++) {
		printf("Bin %d: %f\n", i, gsl_histogram_get(norm_hist, i));
	}

	gsl_rng_free(rnd);
	gsl_histogram_free(norm_hist);

	return 0;
}
