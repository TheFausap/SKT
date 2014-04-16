// Module for constructing complete, standard, orthogonal unitary bases
// for arbitrary SU(d)
// (Generalized Pauli matrices)

#ifndef basis1_h__
#define basis1_h__

#include <map>
#include <cmath>
#include <iterator>
#include <deque>
#include <iostream>
#include <vector>
#include <complex>
#include <armadillo>
#include <string>
#include <numeric>


using namespace arma;
using namespace std;

typedef Oper basOper;
typedef std::pair<int, int> dims;
typedef std::tuple<std::string, int, int> basIdx;
typedef std::map<basIdx, basOper> elemBasis;  // like python dictionaries
typedef std::map<dims, basOper> unitaryBasis; // like python dictionaries


class Basis {
public:
	int d;
	elemBasis elB;
	unitaryBasis uB;
	dims unitary_idKey;
	std::string keys_minus_identity;
	basIdx identity_key;
	basOper identity;

	Basis(int, elemBasis, basIdx);
	Basis(int, unitaryBasis, dims);
	Basis();
	basOper get(basIdx);
	basOper get(dims);
	std::deque<elemBasis> items();
	void print_string();
};

Basis::Basis(int d1, elemBasis eBB, basIdx idKey) {
	d = d1;
	elB = eBB;
	identity_key = idKey;
	identity = eBB[idKey];
};

Basis::Basis(int d1, unitaryBasis l_ub, dims dd) {
	d = d1;
	uB = l_ub;
	unitary_idKey = dd;
	identity = l_ub[dd];
}

Basis::Basis() {
};

basOper Basis::get(basIdx iid) {
	return(elB[iid]);
};

basOper Basis::get(dims iid) {
	return(uB[iid]);
};

void Basis::print_string() {
	cout << "SU(" + std::to_string(d) + ") Basis" << endl;
	for (elemBasis::iterator it = elB.begin(); it != elB.end(); ++it) {
		it->second.print();
		it->second.print_matrix();
	};
};

//////////////////////////////////////////////////////////////////////////////
// Get the standard vector basis for R^n, that is, n n - dimensional vectors
// { v[i] } where v[i] has a 1 in the ith component and 0 everywhere else

std::deque<cx_rowvec> get_standard_vector_basis(int d) {
	cx_rowvec new_vec(d);
	std::deque<cx_rowvec> vv{};
	std::vector<int> rr(d);

	rr = utils::range(0, d);
	for (auto i : rr) {
		new_vec.zeros();
		new_vec(i) = 1;
		vv.push_back(new_vec);
	};
	return(vv);
};


Basis get_hermitian_basis(int d) {
	// Returns a dictionary of generalized hermitian
	// Pauli matrices as a basis for SU(d), i.e.Gell - Mann matrices
	// i.e. d x d orthonormal matrices that are a complete basis for C^{ dxd }
	// From http ://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices#Construction

	elemBasis B;
	Basis B1;
	dims d1;
	basIdx bId;
	cx_mat op_matrix, Id1, E_kj, E_jk;
	basOper h_k_d1;
	std::string op_name;
	std::deque<cx_rowvec> vecbasis{};
	cx_colvec cv1;
	std::vector<int> range_d, range_j, range_k;

#ifdef _DEBUG
	cout << "Getting Hermitian basis for SU(" + std::to_string(d) + ")" << endl;
#endif

	// Base case of recursion, just pass back h_{ 1, 1 }
	if (d == 1) {
		bId = std::make_tuple("h", 1, 1);
		B[bId] = get_identity(1);
		Basis bas1(d, B, bId);
		return(bas1);
	};

	B1 = get_hermitian_basis(d - 1);

	// Generate the h functions
	// # h_{1,d} is always I_d
	op_matrix.eye(d, d);
	op_name = "h_{" + std::to_string(d) + "}";
	bId = std::make_tuple("h", 1, d);
	B[bId] = Oper(op_name, op_matrix);

	// h_{k,d} for 1 < k < d
	range_d = utils::range(2, d);
	for (auto k : range_d) {
		bId = std::make_tuple("h", k, d - 1);
		h_k_d1 = B1.get(bId);
		op_matrix = utils::matrix_direct_sum(h_k_d1.matrix, utils::matrixify(0));
		op_name = "h_{" + std::to_string(k) + "," + std::to_string(d) + "}";
		bId = std::make_tuple("h", k, d);
		B[bId] = Oper(op_name, op_matrix);
	};

	// h_{d,d}
	double scale = sqrt(2.0 / (d*(d - 1)));

#ifdef _DEBUG
	cout << "scale= " + std::to_string(scale) << endl;
#endif

	Id1.eye(d - 1, d - 1);
	op_matrix = scale*utils::matrix_direct_sum(Id1, utils::matrixify(1 - d));
	op_name = "h_{" + std::to_string(d) + "," + std::to_string(d) + "}";
	bId = std::make_tuple("h", d, d);
	B[bId] = Oper(op_name, op_matrix);

#ifdef _DEBUG
	cout << op_name << "= " << op_matrix << endl;
#endif

	// It is easier to generate E_jk with outer product of standard vector basis
	vecbasis = get_standard_vector_basis(d);

	// Generate the f functions
	// First, for 1 < k < j <= d
	// subtract one from vector indices since math indices begin at 1
	// but computer indices begin at 0
	range_j = utils::range(1, d + 1);
	for (auto j : range_j) {
		range_k = utils::range(1, j);
		for (auto k : range_k) {
			cv1 = vecbasis[k - 1].st();  // colvec = transpose of rowvec
			E_kj = cv1*vecbasis[j - 1];  // outer product (ARMADILLO) is done multiplying a colvec * rowvec
			cv1 = vecbasis[j - 1].st();
			E_jk = cv1*vecbasis[k - 1];
			op_matrix = E_kj + E_jk;
			op_name = "f_{" + std::to_string(k) + "," + std::to_string(j) + "}";
			bId = std::make_tuple("f", k, j);
			B[bId] = Oper(op_name, op_matrix);
		}
	};

	// Next, for 1 < k < j <= d
	// subtract one from vector indices since math indices begin at 1
	// but computer indices begin at 0
	range_k = utils::range(1, d + 1);
	for (auto k : range_k) {
		range_j = utils::range(1, k);
		for (auto j : range_j) {
			cv1 = vecbasis[k - 1].st();  // colvec = transpose of rowvec
			E_kj = cv1*vecbasis[j - 1]; // outer product (ARMADILLO) is done multiplying a colvec * rowvec
			cv1 = vecbasis[j - 1].st();
			E_jk = cv1*vecbasis[k - 1];
			op_matrix = complex<double>(0, -1) * (E_kj - E_jk);
			op_name = "f_{" + std::to_string(k) + "," + std::to_string(j) + "}";
			bId = std::make_tuple("f", k, j);
			B[bId] = Oper(op_name, op_matrix);
		}
	};

	bId = std::make_tuple("h", 1, d);
	Basis hermBas(d, B, bId);

	return(hermBas);
};

////////////////////////////////////////////////////////////////////
// Returns a dictionary of generalized unitary
// Pauli matrices as a basis for SU(d)
// i.e. d x d orthonormal matrices that are a complete basis for C^{dxd}
Basis get_unitary_basis(int d) {
	std::vector<int> range_d;
	std::deque<cx_rowvec> v{};
	cx_colvec cvv;
	cx_rowvec rvv;
	cx_mat outer_prod, jkm_term;
	string name;
	unitaryBasis S;
	dims d1;

	range_d = utils::range(0, d);
	v = get_standard_vector_basis(d);

	complex<double> zeta = exp(complex<double>(0, 2) * complex<double>(datum::pi, 0) / complex<double>(d, 0));

#ifdef _DEBUG
	cout << "zeta= " << zeta << endl;
#endif

	for (auto j : range_d) {
		for (auto k : range_d) {
			cx_mat sum(d, d);
			sum.fill(0);
			for (auto m : range_d) {
				cvv = v[m].t();
				outer_prod = cvv * v[(m + k) % d];
#ifdef _DEBUG
				cout << "v[" << std::to_string(m) + "]v[" + std::to_string((m + k) % d) + "]= " << outer_prod << endl;
#endif
				jkm_term = std::pow(zeta, (j*m)) * outer_prod;
#ifdef _DEBUG
				cout << "(" + to_string(j) + "," + to_string(k) + "," + to_string(m) + ")= " << jkm_term << endl;
#endif
				sum = sum + jkm_term;
			}
			name = "S_{" + to_string(j) + "," + to_string(k) + "}";
			d1 = std::make_pair(j, k);
			S[d1] = Oper(name, sum);
#ifdef _DEBUG
			Oper new_op(name, sum);
			new_op.print();
			new_op.print_matrix();
#endif
		}
	}

	d1 = std::make_pair(0, 0);
	Basis uniBas(d, S, d1);
	return(uniBas);
};
#endif // basis1_h__