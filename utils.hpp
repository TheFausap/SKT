#ifndef utils_h__
#define utils_h__

#include <iostream>
#include <complex>
#include <armadillo>
#include <assert.h>
#include <algorithm>
#include <numeric>


using namespace arma;
using namespace std;

namespace utils {

	cx_mat matrixify(complex<double> d) {
		cx_mat a;
		a << d << endr;
		return(a);
	};

//	cx_mat matrixify(std::list<complex<double>> d) {
//		cx_mat a = d;
//		a.reshape(d.size() / 2, d.size() / 2);
//		return(a);
//	};

	std::vector<int> range(int d1, int d2) {
		std::vector <int> tt(d2 - d1);

		if (d1 != d2) {
			std::iota(tt.begin(), tt.end(), d1);
			return(tt);
		}
		else {
			return(tt);
		}
	}

	double md(cx_mat A, cx_mat B)
	{
		complex<double> z11,z22;

		z11=conj(A(0,0))*B(0,0)+conj(A(1,0))*B(1,0);
		z22=conj(A(0,1))*B(0,1)+conj(A(1,1))*B(1,1);
	
		return( pow(abs(z11+z22),2) );
	};

	double md_tri(cx_mat A, cx_mat B)
	{
		return( sqrt(2-sqrt(md(A,B)))/2);
	};

	double operator_norm(cx_mat A)
	{
		vec eigvals(A.n_rows);

		eigvals = eig_sym(A);

		return( max(max(abs(eigvals),0)) );
	};

	double trace_norm(cx_mat A)
	{
		complex<double> tr;

		tr = trace( (A * A.t()) ); // A is complex matrix so .t() is Hermitian transpose.

		return( real(sqrt(tr)) );
	};

	double trace_distance(cx_mat A, cx_mat B)
	{
		int sz = A.n_cols;
		double pp;
		cx_mat mdiff(sz,sz), mdiag(sz,sz), product(sz,sz);
		vec tvals(A.n_rows);

		tvals.fill(0.0);

		mdiff = A - B;
		mdiag = mdiff.t();          // The same as before.
		product = mdiff*mdiag;
		tvals = eig_sym(product);
		pp = norm(tvals,"fro");
		return(pp);
	};

	double fowler_distance(cx_mat A, cx_mat B)
	{
		int sz = A.n_cols;
		complex<double> tr;
		double frac;
	
		cx_mat prod(sz,sz);
		assert(A.n_cols==B.n_cols);
		prod = A.t() * B;
		tr = trace(prod);
		frac = (1.0*(sz - abs(tr))) / sz;
		return(sqrt(abs(frac)));
	};

	cx_mat matrix_direct_sum(cx_mat A, cx_mat B)
	{
		int sz = A.n_cols+B.n_cols;
		cx_mat C(sz,sz);
	
		C.fill(0.0);
			C.submat(0, 0, A.n_rows - 1, A.n_cols - 1) = A;
			C.submat(A.n_rows, A.n_cols, sz - 1, sz - 1) = B;
		return(C);
	};

	string dagger_and_simplify(string name) {
		size_t name_len = name.length();

		if (name_len < 1) return(name);
		if (name[name_len-1] == 'd') return(name.substr(0,name_len-1));
		string new_name = name + "d";
		return(new_name);
	};


}; // end of namespace utils
#endif // utils_h__