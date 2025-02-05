#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
gsl_complex Zta(const double s, const long long max) {
	gsl_complex sum = gsl_complex_rect(0.0, 0.0);
	for (int n = 1; n <= max; ++n) {
		const double ex = -s * log(n);
		gsl_complex term = gsl_complex_polar(1.0 / pow(n, s), ex);
		sum = gsl_complex_add(sum, term);
	}
	return sum;
}
double Theta (const double s) {
	gsl_complex z = gsl_complex_rect(0.25, s / 2.0);
	gsl_sf_result lnr, lni;
	gsl_sf_lngamma_complex_e(GSL_REAL(z), GSL_IMAG(z), &lnr, &lni);
	const double G = atan2(lni.val, lnr.val);
	double theta_t = G - (s / 2.0) * log(M_PI);
	return theta_t;
}
double Zeta(const int s, const long long max) {
	if (s > 1) {
		return gsl_sf_zeta(s);
	}
	else if (s == 1) {
		double res = 0.0;
		res = 0.0;
		for (long long n = 1; n <= max; ++n) {
			res += 1.0/n;
		}
		return res;
	}
	gsl_complex zeta = Zta(s, max);
	const double TT = Theta (s);
	double Z = GSL_REAL(zeta) * cos(TT) + GSL_IMAG(zeta) * sin(TT);
	return Z; 
}
int main (int argc, char **argv) {
	if (!(argc == 4 || argc == 3)) {
		fprintf(stderr, "Error: to many or few arguments. Usage: rmn [S] [MAX] *[PREC]\n");
		return 1;
	}
	const int S = atoi(argv[1]);
	const long long MAX = atol(argv[2]);
	const double R = Zeta(S, MAX);
	if (argc == 4) {
		const int PREC = atoi(argv[3]);
		printf("%.*f\n", PREC, R);
		return 0;
	}
	printf("%f\n", R);
	return 0;
}

