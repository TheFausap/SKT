// Oper class 
#ifndef operator_h__
#define operator_h__

#include <iostream>
#include <map>
#include <armadillo>


using namespace arma;
using namespace std;

class Oper {
public:
	cx_mat matrix;
	string name;
	string ancestors;

	Oper (string,cx_mat,string);
	Oper (string,string);
	Oper (string,cx_mat);
	Oper ();
	void print_matrix() { cout << "Matrix: " << matrix << endl; }
	void print() { cout << "Operator: " << name << endl << "  Ancestors: " << ancestors << endl; }
	Oper add_ancestor(Oper, string);
	void matrix_from_ancestors(Oper[], Oper, float);
	Oper multiply(Oper,string);
	Oper dagger();
	Oper scale(int,string);
	Oper scale(int);
	bool operator==(const Oper&);
};

typedef Oper tOper;

tOper T, I2, H, SZ, SY, SX, T_inv;
cx_mat H_matrix(2, 2), SX_matrix(2, 2), SY_matrix(2, 2), SZ_matrix(2, 2);
cx_mat T_matrix;

Oper::Oper (string n, cx_mat a, string anc) {
	name = n;
	matrix = a;
	ancestors = anc;
};

Oper::Oper (string n, string anc) {
	name = n;
	matrix = zeros<cx_mat>(0,0);
	ancestors = anc;
};

Oper::Oper (string n, cx_mat a) {
	name = n;
	matrix = a;
	ancestors = n;
};

Oper::Oper () {
};

Oper Oper::add_ancestor(Oper other, string new_name) {
	string new_ancestor = ancestors + other.ancestors;
	Oper new_op (new_name,new_ancestor);
	return(new_op);
};

void Oper::matrix_from_ancestors(Oper isect_dict[], Oper identity, float tolerance) {
	double dist;
	
	for ( string::iterator it=ancestors.begin(); it!=ancestors.end(); ++it) {
		matrix = matrix * isect_dict[*it].matrix;
	};
	matrix = identity.matrix;
	dist = utils::fowler_distance(matrix,identity.matrix);
};

bool Oper::operator==(const Oper &b) {
	return ( (name == b.name) && (ancestors == b.ancestors) );
};

Oper Oper::multiply(Oper other, string new_name) {
	cx_mat new_matrix = matrix * other.matrix;
	string new_ancestors = ancestors + other.ancestors;
	Oper new_op ("",new_matrix,new_ancestors);
	return(new_op);
};

Oper Oper::dagger() {
	cx_mat new_matrix = matrix.t();
	string rev_ancestors = ancestors; 
	string new_ancestors = "";
	string::iterator it;

	reverse(rev_ancestors.begin(),rev_ancestors.end());
	for (it=rev_ancestors.begin(); it!=rev_ancestors.end(); ++it) {
		new_ancestors.append(utils::dagger_and_simplify(ancestors));
	};
	string new_name = utils::dagger_and_simplify(name);
	Oper new_op(new_name,new_matrix,new_ancestors);
	return(new_op);
};

Oper Oper::scale(int scalar, string new_name) {
	cx_mat new_matrix = matrix * scalar;
	string new_ancestors = ancestors + to_string(scalar);
	Oper new_op(new_name,new_matrix,new_ancestors);
	return(new_op);
};

Oper Oper::scale(int scalar) {
	cx_mat new_matrix = matrix * scalar;
	string new_ancestors = ancestors + to_string(scalar);
	Oper new_op("",new_matrix,new_ancestors);
	return(new_op);
};

Oper get_identity(int d) {
	Oper new_op("I",eye<cx_mat>(d,d),"I");
	return(new_op);
};

/////////////////////////////////////////////////////////
// SU(2) constants

void initOperConstants() {

	complex<double> J = complex<double>(0, 1);
	double Hcoeff = 1 / sqrt(2.0);
	complex<double> Tcoeff = exp(J*(datum::pi / 4.0));

	// 2x2 identity matrix
	Oper locI2 = get_identity(2);
	I2 = locI2;

	// Pauli X matrix
	SX_matrix << 0 << 1 << endr
		<< 1 << 0 << endr;
	Oper locSX("SX", SX_matrix);
	SX = locSX;

	// Pauli Y matrix
	SY_matrix << 0 << -J << endr
		<< J << 0 << endr;
	Oper locSY("SY", SY_matrix);
	SY = locSY;

	// Pauli Z matrix
	SZ_matrix << 1 << 0 << endr
		<< 0 << -1 << endr;
	Oper locSZ("SZ", SZ_matrix);
	SZ = locSZ;

	// Hadamard gate
	H_matrix << 1 << 1 << endr
		<< 1 << -1 << endr;
	Oper locH("H", Hcoeff*H_matrix);
	H = locH;

	// pi/8 gate
	T_matrix << 1 << 0 << endr
		<< 0 << Tcoeff << endr;
	Oper locT("T", T_matrix);
	T = locT;

	// Inverse pi/8 gate
	Oper locT_inv("Td", T_matrix.i());
	T_inv = locT_inv;
};

#endif // operator_h__
