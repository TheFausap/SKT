// Matrix approximation using a fixed set of basic matrices
//

#include <armadillo>
#include <stdexcept>
#include "simplify.hpp"
#include <ctime>
#include "config.hpp"

typedef std::map<std::string, tOper> tIsetDict;

class BasicApproxSettings {
public:
	Oper identity;
	tIsetDict iset_dict;
	tArrayOp iset;
	ruleSet rSet;
	SimplifyEngine sse;

	BasicApproxSettings::BasicApproxSettings() {
		cx_mat matrix;
		matrix = zeros<cx_mat>(0, 0);
		Oper identity ("N",matrix);
		iset = {};
		iset_dict.clear();
		rSet = {};
	};

	void set_iset(tArrayOp new_iset) {
		check_iset(new_iset);
		iset = new_iset;
		for (auto i : new_iset) iset_dict[i.name] = i;
	};

	void set_identity(Oper new_identity) {
		identity = new_identity;
		iset_dict[new_identity.name] = new_identity;
	};

	void print_iset() {
		cout << "INSTRUCTION SET" << endl;
		for (auto i : iset) i.print();
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

	bool simplify_new(BasicApproxSettings &ss1, tOper &new_op) {
		// Have to convert a string of ancestors in array of ancestors operators
		std::string ancs = new_op.ancestors;
		tOper new_ancestor;
		tArrayOp arrayAncestors;
		std::string simplified_ancestors = "";

		// Array creation from ancestors label
		size_t l_ancs = ancs.length();
		for (int i = 0; i < l_ancs; i++) {
			for (auto j : baseOpSet) {
				if (ancs[i] == j.name) arrayAncestors.push_back(j) ;
			};
		};

		// calling simplify method
		tSimplified ancsimp = ss1.simplify(arrayAncestors);
		
	};

	void gen_basic_approx_generation(BasicApproxSettings &ss1, tArrayOp &s1) {
		//reset_global_sequences();
		//reset_generation_stats();
		extern bool simplify_new;

		tArrayOp simplified_ancestors = {};
		for (auto i : s1) {
			for (auto insn : iset) {
				tOper new_op = i.add_ancestor(insn,"");
				bool already_done = simplify_new(ss1, new_op);
			};
		};
	};

	// Generate table of basic approximations as preprocessing
	// ll_0 - fixed length of sequences to generate for preprocessing table
	void basic_approxes(int &ll0, BasicApproxSettings &sett) {
		BasicApproxSettings settings = sett;

		//reset_global_stats();
		//set_filename_suffix("g1");
		gen_basic_approx_generation(sett, { settings.identity });
		//print_generation_stats(1);
	};

	void generate_approxes(int l0, BasicApproxSettings setts) {
		time_t begin_time, end_time;
		double seconds;

		// set_filename_suffix("iset");
		setts.print_iset();
		// dump_to_file(settings.iset);

		// Start generation timer
		time(&begin_time);
		basic_approxes(l0, setts);
		time(&end_time);
		seconds = difftime(end_time, begin_time);
		cout << "Generation time: " + to_string(seconds / 60) + " minutes." << endl;
	};
};

