// TEST FOR SIMPLIFY ENGINE
//
#include <string>
#include <armadillo>
#include "config.hpp"

void main() {

	ruleSet rSet;

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

	BasicApproxSettings settings;

	settings.set_iset(iset2);
	settings.init_simplify_engine(rSet);
	settings.set_identity(I2);
	//settings.basis = H2;

	//tSimplified pp = sse.simplify(t8);
	//cout << pp.first << endl;

};