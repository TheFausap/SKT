// Matrix approximation using a fixed set of basic matrices
//

#include <armadillo>
#include <stdexcept>
#include "simplify.hpp"

class BasicApproxSettings {
public:
	Oper identity;
	elemBasis iset_dict;
	tArrayOp iset;
	ruleSet rSet;
	SimplifyEngine sse;

	BasicApproxSettings::BasicApproxSettings() {
		cx_mat matrix;
		matrix = zeros<cx_mat>(0, 0);
		Oper identity ("N",matrix);
		iset = {};
		rSet = {};
	};

	void check_iset(tArrayOp is) {
		size_t m = is.size();
		cout << to_string(m) + " instructions found" << endl;
		tOper first_op = is[0];
		int first_cols = first_op.matrix.n_cols;
		int first_rows = first_op.matrix.n_rows;
		if ((first_cols == 0) || (first_rows == 0)) domain_error("First operator is not a matrix!");
		if (first_rows != first_cols) domain_error("First operator is not a square matrix");
		for (size_t i = 0; i < m; i++) {
			int i_cols = is[i].matrix.n_cols;
			if (i_cols != first_cols) domain_error("Operator" + to_string(i) + "'s shape does not match first shape!");
		};
	};

	SimplifyEngine init_simplify_engine(ruleSet rSet) {
		SimplifyEngine sse(rSet);
		return(sse);
	}

	tSimplified BasicApproxSettings::simplify(tArrayOp seq) { 
		return(sse.simplify(seq));
	};
};