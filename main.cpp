// TEST FOR SIMPLIFY ENGINE
//
#include <string>
#include <armadillo>
#include "config.hpp"

void main() {

	ruleSet rSet;
	//tArrayOp iset2 = { "H","T","T_inv" };

	Basis H2;

	initOperConstants();

	tArrayOp iset2 = { H, T, T_inv };
	tArrayOp t8 = { T, T, T, T, T, T, T, T };
	tArrayOp Td8 = { T_inv, T_inv, T_inv, T_inv, T_inv, T_inv, T_inv, T_inv };

	// Simplifying rules
	ProductFactory pRule;
	SimplifyRule *identity_rule = pRule.Make(0, {});
	SimplifyRule *double_H_rule = pRule.Make(1, { H });
	SimplifyRule *adjoint_rule = pRule.Make(2, {});
	SimplifyRule *T8_rule = pRule.Make(3, t8);
	SimplifyRule *Td8_rule = pRule.Make(3, Td8);

	rSet = { identity_rule, double_H_rule, adjoint_rule, T8_rule, Td8_rule };
	
	//Td8_rule->print();
	
	//H2 = get_hermitian_basis(2);

	//T_inv.print();
	//T_inv.print_matrix();

	SimplifyEngine sse(rSet);
	tSimplified pp = sse.simplify(t8);
	cout << pp.first << endl;

};