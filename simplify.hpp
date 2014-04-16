// Simplification engine for matrix operators
//
#ifndef simplify_h__
#define simplify_h__

#include <string>
#include <iostream>
#include <deque>
#include <vector>
#include "config.hpp"

using namespace std;
typedef std::vector<tOper> tArrayOp;
typedef std::pair<bool, tArrayOp> tRuleOut;
typedef std::pair<size_t, tArrayOp> tSimplified;

class SimplifyRule {
public:
	std::string slogan, id_sym;
	size_t arg_count;

	SimplifyRule::SimplifyRule(std::string slo, int arg_c, std::string ii) {
		slogan = slo;
		arg_count = arg_c;
		id_sym = ii;
	};

	SimplifyRule::SimplifyRule(tArrayOp aOp, std::string new_sym) {
		slogan = "";
		id_sym = new_sym;
		arg_count = aOp.end() - aOp.begin();
		for (auto arg : aOp) slogan += arg.name;
		slogan += " = " + id_sym;
	};

	SimplifyRule::SimplifyRule(tArrayOp aOp) {
		slogan = "";
		id_sym = "I";
		arg_count = aOp.end() - aOp.begin();
		for (auto arg : aOp) slogan += arg.name;
		slogan += " = " + id_sym;
	};

	// this doesn't work in c++, cause is statically typed language.
	// in python works because is dynamically typed language.
	// tRuleOut simplify() { <subclass method> _simplify__(); }
	// FACTORY!

	virtual tRuleOut simplify(tArrayOp pp) {
		return(make_pair(false, pp ));
	};

	void print() {
		cout << "Name : " + slogan + " #### " << endl;
	};
};

typedef std::deque<SimplifyRule *> ruleSet;
typedef SimplifyRule tSimpRule;

class IdentityRule: public SimplifyRule {
public:
	std::string id_sym;

	IdentityRule::IdentityRule(int d) : SimplifyRule("I*Q = Q", 2,"I") {
		id_sym = "I";
	};

	tRuleOut IdentityRule::simplify(tArrayOp);
};

tRuleOut IdentityRule::simplify(tArrayOp OpArr) {
	bool activated = false;
	tOper A, B, C;
	A = OpArr[0];
	B = OpArr[1];
	//C = "";
	if (A.name == id_sym) {
		activated = true;
		C = B;
	}
	else if (B.name == id_sym) {
		activated = true;
		C = A;
	}

	// This code is the same for each simplification rule
	if (activated) {
		for (int i = 0; i < arg_count; i++) OpArr.erase(OpArr.begin());
		OpArr.insert(OpArr.begin(), C);
	}
	return(make_pair(activated,OpArr));
};

class DoubleIdentityRule: public SimplifyRule {
public:
	std::string id_sym;
	tOper symbol;

	DoubleIdentityRule::DoubleIdentityRule(tOper symb) : SimplifyRule("Q*Q = I", 2, "I") {
		symbol = symb;
		id_sym = "I";
	};

	tRuleOut DoubleIdentityRule::simplify(tArrayOp);
};

tRuleOut DoubleIdentityRule::simplify(tArrayOp OpArr) {
	bool activated = false;
	tOper A, B, C;
	A = OpArr[0];
	B = OpArr[1];
	//C = "";
	if ((A.name == symbol.name) && (B.name == symbol.name)){
		activated = true;
		// creates an ad-hoc identity matrix with the right number
		// of rows and cols (square matrix by definition rows=cols).
		C = Oper (id_sym,eye<cx_mat>(A.matrix.n_rows,A.matrix.n_rows));
	}

	// This code is the same for each simplification rule
	if (activated) {
		for (int i = 0; i < arg_count; i++) OpArr.erase(OpArr.begin());
		OpArr.insert(OpArr.begin(), C);
	}
	return(make_pair(activated, OpArr));
};

class AdjointRule: public SimplifyRule {
public:
	std::string id_sym;

	AdjointRule::AdjointRule(int d) : SimplifyRule("Q*Q\\dagger = I", 2, "I") {
		id_sym = "I";
	};

	tRuleOut AdjointRule::simplify(tArrayOp);
};

tRuleOut AdjointRule::simplify(tArrayOp OpArr) {
	bool activated = false;
	tOper A, B, C;
	A = OpArr[0];
	B = OpArr[1];
	// C = "";
	if ((A.name	== B.name.substr(B.name.size()-1)) && (B.name.back() == 'd')) {
		activated = true;
		C = Oper(id_sym, eye<cx_mat>(A.matrix.n_rows, A.matrix.n_rows));
	}
	else if ((B.name == A.name.substr(A.name.size() - 1)) && (A.name.back() == 'd')) {
		activated = true;
		C = Oper(id_sym, eye<cx_mat>(A.matrix.n_rows, A.matrix.n_rows));
	}

	// This code is the same for each simplification rule
	if (activated) {
		for (int i = 0; i < arg_count; i++) OpArr.erase(OpArr.begin());
		OpArr.insert(OpArr.begin(), C);
	}
	return(make_pair(activated, OpArr));
};

class GeneralRule: public SimplifyRule {
public:
	std::string id_sym;
	tArrayOp sequence;

	GeneralRule::GeneralRule(tArrayOp seqs, std::string new_sym) : SimplifyRule(seqs, new_sym) {
		id_sym = new_sym;
		sequence = seqs;
	};

	GeneralRule::GeneralRule(tArrayOp seqs) : SimplifyRule(seqs) {
		id_sym = "I";
		sequence = seqs;
	};

	tRuleOut GeneralRule::simplify(tArrayOp);
};

tRuleOut GeneralRule::simplify(tArrayOp OpArr) {
	bool activated = false;
	tArrayOp nullOpArr;
	tOper C;

	nullOpArr.clear();

	for (int i = 0; i < arg_count; i++) {
		if (sequence[i].name != OpArr[i].name) return(make_pair(false,nullOpArr));
	};

	activated = true;
	C = Oper(id_sym, eye<cx_mat>(OpArr[0].matrix.n_rows, OpArr[0].matrix.n_rows));

#ifdef _DEBUG
	cout << "GeneralRule._simplify__: "; 
	for (auto i: OpArr) cout << i.name;
	cout << " -> " + id_sym << endl;
#endif

	// This code is the same for each simplification rule
	if (activated) {
#ifdef _DEBUG
		cout << slogan + "(" + to_string(arg_count) + ")" + " OBTAINS!" << endl;
#endif
		for (int i = 0; i < arg_count; i++) OpArr.erase(OpArr.begin());
		OpArr.insert(OpArr.begin(), C);
#ifdef _DEBUG
		for (auto i : OpArr) cout << i.name;
		cout << endl;
#endif
	};
	return(make_pair(activated, OpArr));
};

class SimplifyEngine {
public:
	ruleSet rs;
	size_t max_arg_count;

	SimplifyEngine::SimplifyEngine(ruleSet rs_1) {
		rs = rs_1;
		max_arg_count = 0;
		for (auto rule : rs_1) {
			if (rule->arg_count > max_arg_count) max_arg_count = rule->arg_count;
		};
	};

	SimplifyEngine::SimplifyEngine() {
		max_arg_count = 0;
	}

	std::pair<bool, tArrayOp> transfer_to_scratch(tArrayOp &sequence, tArrayOp &scratch) {
		size_t sequence_len = sequence.size();
		if (sequence_len > 0) {
			tOper new_op = sequence.back();
			sequence.pop_back();
			scratch.insert(scratch.begin(), new_op);
			return(make_pair(true, sequence));
		}
		else {
			return(make_pair(false, sequence));
		};
#ifdef _DEBUG
		cout << "transfer to scratch ";
		for (auto i : scratch) cout << i.name;
		cout << endl;
#endif
	};

	tArrayOp fill_scratch_sequence(tArrayOp &sequence, tArrayOp &scratch) {
		bool long_enough = true;
		std::pair<bool, tArrayOp> oo;

		while (long_enough && (scratch.size() < max_arg_count)) {
			oo = transfer_to_scratch(sequence, scratch);
		};
		return(sequence);
	};

	tSimplified simplify(tArrayOp &sequence) {
		size_t simplify_length = sequence.size();
		tArrayOp scratch_sequence {};
		tRuleOut resRule;

		bool global_obtains = true;

		while (global_obtains) {
			global_obtains = false;
			sequence = fill_scratch_sequence(sequence, scratch_sequence);
#ifdef _DEBUG
			cout << "sequence= ";
			for (auto i : sequence) cout << i.name;
			cout << endl << "fill_scratch= ";
			for (auto i : scratch_sequence) cout << i.name;
			cout << endl;
#endif
			for (auto rule : rs) {
				if (scratch_sequence.size() >= rule->arg_count) {
					size_t split = scratch_sequence.size() - rule->arg_count;
					tArrayOp scratch_excess{}, scratch_subset{};
					for (size_t i = 0; i < split;i++) scratch_excess.push_back(scratch_sequence[i]);
					for (size_t i = split; i < scratch_sequence.size(); i++) scratch_subset.push_back(scratch_sequence[i]);
					resRule = rule->simplify(scratch_subset);
					scratch_subset = resRule.second;
					for (size_t i = 0; i < scratch_excess.size(); i++) scratch_sequence[i] = scratch_excess[i];
					for (size_t i = 0; i < scratch_subset.size(); i++) scratch_sequence[i + scratch_excess.size()] = scratch_subset[i];
					
					if (resRule.first) {
#ifdef _DEBUG
						cout << "*** ";
						for (auto i : scratch_sequence) cout << i.name;
						cout << endl;
#endif
						global_obtains = true;
					};
				};
			};
#ifdef _DEBUG
			cout << "globals= " + to_string(resRule.first) << endl;
#endif
		};
		for (size_t i = 0; i < scratch_sequence.size(); i++) sequence.push_back(scratch_sequence[i]);
		simplify_length -= sequence.size();

		return(make_pair(simplify_length, sequence));
	};
};

class ProductFactory
{
public:
	virtual SimplifyRule *Make(int type, tArrayOp OP)
	{
		switch (type)
		{
		case 0:
			return new IdentityRule(0);
		case 1:
			return new DoubleIdentityRule(OP.front());
		case 2:
			return new AdjointRule(0);
		case 3:
			return new GeneralRule(OP);
		default:
			return new IdentityRule(0);
		};
	};
};

#endif // simplify_h__