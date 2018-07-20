#pragma once

#include <exception>
#include <vector>
#include <stack>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <memory>

#include "TOKEN.h"

class INTERPRETER {
	std::vector<std::shared_ptr<TOKEN>> tokens;
	std::vector<std::shared_ptr<TOKEN>> rpn; //reverse polish notation
	std::string expression;

public:
	INTERPRETER(std::string expression) : expression(expression) {
		interpret();
	}
	void interpret() {
		//std::cout << "INPUT: " << expression << std::endl;

		removeBlanks();
		interpretExpression();
		cleanUpForShuntingYard();

		std::cout << "TOKENS: \t";
		for (unsigned int i = 0; i < tokens.size(); i++) {
			std::cout << *(tokens[i]) << " ";
		}
		std::cout << std::endl;

		shuntingYard();

		/*std::cout << "RPN: ";
		for (unsigned int i = 0; i < rpn.size(); i++) {
			std::cout << *(rpn[i]) << " ";
		}
		std::cout << std::endl;*/
	}

private:
	void error(std::string errorstring, unsigned int pos) {
		std::stringstream tmp;
		tmp << errorstring << " at position: " << std::to_string(pos + 1);
		throw tmp.str();
	}

	void error(std::string errorstring) {
		throw errorstring;
	}

	void removeBlanks() {
		expression.assign(expression.begin(), remove_if(expression.begin(), expression.end(), &isspace));
	}

	bool isOperator(char token) {
		return token == '+' || token == '-' || token == '*' || token == '/' || token == '^';
	}

	bool isNumber(char token) {
		return token == '0' || token == '1' || token == '2' || token == '3' || token == '4' || token == '5' || token == '6' || token == '7' || token == '8' || token == '9';
	}

	bool isLetter(char token) {
		return	token == 'a' || token == 'b' || token == 'c' || token == 'd' || token == 'e' || token == 'f' || token == 'g' || token == 'h' || token == 'i' || token == 'j' || token == 'k' ||
			token == 'l' || token == 'm' || token == 'n' || token == 'o' || token == 'p' || token == 'q' || token == 'r' || token == 's' || token == 't' || token == 'u' || token == 'v' ||
			token == 'w' || token == 'x' || token == 'y' || token == 'z' || token == 'A' || token == 'B' || token == 'C' || token == 'D' || token == 'E' || token == 'F' || token == 'G' ||
			token == 'H' || token == 'I' || token == 'J' || token == 'K' || token == 'L' || token == 'M' || token == 'N' || token == 'O' || token == 'P' || token == 'Q' || token == 'R' ||
			token == 'S' || token == 'T' || token == 'U' || token == 'V' || token == 'W' || token == 'X' || token == 'Y' || token == 'Z';
	}

	bool isFunction(std::string test) {
		return test == "exp" || test == "sin" || test == "cos" || test == "tan" || test == "arcsin" || test == "arccos" || test == "arctan" || test == "abs" || test == "sinh" || test == "cosh" || test == "tanh" || test == "arsinh" || test == "arcossh" || test == "artanh" || test == "min" || test == "max" || test == "ln" || test == "log10" || test == "log2";
	}

	bool interpretExpression() {
		char token;
		for (unsigned int i = 0; i < expression.size(); i++) {
			token = expression[i];

			//OPERATOR:
			if (isOperator(token)) {
				tokens.push_back(std::shared_ptr<OPERATOR>(new OPERATOR(token)));
				continue;
			}

			//SCALAR:
			if (isNumber(token) || token == '.') {
				std::string numberstring;
				bool decimal_separator = false;
				do {
					if (token == '.') {
						if (!decimal_separator) {
							decimal_separator = true;
							if (numberstring.empty()) {
								numberstring.push_back('0');
							}
						}
						else {
							error("SCALAR HAS MORE THAN ONE DECIMAL SEPERATOR!", tokens.size());
							return false;
						}
					}
					numberstring.push_back(token);
					if (++i >= expression.size()) {
						break;
					}
					token = expression[i];
				} while (isNumber(token) || token == '.');
				if (token == '.') {
					error("SCALAR ENDS WITH A DECIMAL SEPERATOR!", tokens.size());
					return false;
				}

				tokens.push_back(std::shared_ptr<SCALAR>(new SCALAR(std::stod(numberstring))));

				i--;
				continue;
			}

			//FUNCTION or VARIABLE:
			if (isLetter(token)) {
				std::string letterstring;
				do {
					letterstring.push_back(token);
					if (++i >= expression.size()) {
						break;
					}
					token = expression[i];
				} while (isLetter(token));

				//FUNCTION:
				if (isFunction(letterstring)) {
					tokens.push_back(std::shared_ptr<FUNCTION>(new FUNCTION(letterstring)));
				}

				//VARIABLE:
				else {
					tokens.push_back(std::shared_ptr<VARIABLE>(new VARIABLE(letterstring)));
				}

				i--;
				continue;
			}

			//ARGUMENT_SEPARATOR:
			if (token == ',') {
				tokens.push_back(std::shared_ptr<ARGUMENT_SEPARATOR>(new ARGUMENT_SEPARATOR));
				continue;
			}

			//BRACKETS:
			if (token == '(') {
				tokens.push_back(std::shared_ptr<OPENING_BRACKET>(new OPENING_BRACKET));
				continue;
			}
			if (token == ')') {
				tokens.push_back(std::shared_ptr<CLOSING_BRACKET>(new CLOSING_BRACKET));
				continue;
			}

			std::stringstream tmp;
			tmp << "TOKEN '" << token << "' DOES NOT BELONG TO ANY CATEGORY!";
			error(tmp.str(), tokens.size());
			return false;
		}
		return true;
	}

	bool cleanUpForShuntingYard() {
		if (tokens.empty()) {
			tokens.push_back(std::shared_ptr<SCALAR>(new SCALAR(0)));
		}

		//TEST ON SIGN OPERATOR AT POSITION 0: like -2..., -x..., -(..., -sin..., --...
		bool rerun = true;
		while (rerun) {
			rerun = false;
			if (tokens[0]->getType() == 0) {
				std::shared_ptr<OPERATOR> a = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[0]);
				if (a->getOPERATOR() == '-' || a->getOPERATOR() == '+') {
					if (tokens.size() > 1) {
						rerun = true;
						if (a->getOPERATOR() == '-') {
							if (tokens[1]->getType() == 1) {
								std::shared_ptr<SCALAR> b = std::static_pointer_cast<SCALAR, TOKEN>(tokens[1]);
								tokens[1] = std::shared_ptr<SCALAR>(new SCALAR(-b->getSCALAR()));
								b.reset();
							}
							else if (tokens[1]->getType() == 2 || tokens[1]->getType() == 3 || tokens[1]->getType() == 10) {
								tokens.insert(tokens.begin() + 1, std::shared_ptr<SCALAR>(new SCALAR(-1)));
								tokens.insert(tokens.begin() + 2, std::shared_ptr<OPERATOR>(new OPERATOR('*')));
							}
							else if (tokens[1]->getType() == 0) {
								std::shared_ptr<OPERATOR> b = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[1]);
								if (b->getOPERATOR() == '-') {
									tokens.insert(tokens.begin() + 1, std::shared_ptr<OPERATOR>(new OPERATOR('+')));
									b.reset();
									tokens.erase(tokens.begin() + 2);
								}
								else if (b->getOPERATOR() == '+') {
									tokens.insert(tokens.begin() + 1, std::shared_ptr<OPERATOR>(new OPERATOR('-')));
									b.reset();
									tokens.erase(tokens.begin() + 2);
								}
								else {
									std::stringstream tmp;
									tmp << "OPERATOR '" << b->getOPERATOR() << "' IS NOT A SIGN OPERATOR / TWO NOT COMPATIBLE TOKENS '" << a->getOPERATOR() << "' AND '" << b->getOPERATOR() << "' MEET!";
									error(tmp.str(), 0);
									return false;
								}
							}
							else if (tokens[1]->getType() == 5 || tokens[1]->getType() == -10) {
								std::stringstream tmp;
								tmp << "TWO NOT COMPATIBLE TOKENS '" << a->getOPERATOR() << "' AND '" << (tokens[1]->getType() == 5 ? ',' : ')') << "' MEET!";
								error(tmp.str(), 0);
								return false;
							}
							else {
								error("NOT REACHABLE!");
								return false;
							}
						}
						a.reset();
						tokens.erase(tokens.begin() + 0);
					}
					else {
						std::stringstream tmp;
						tmp << "EXPRESSION ENDS WITH A SIGN OPERATOR: '" << a->getOPERATOR() << "'!";
						error(tmp.str(), 0);
						return false;
					}
				}
				else {
					std::stringstream tmp;
					tmp << "OPERATOR '" << a->getOPERATOR() << "' IS NOT A SIGN OPERATOR!";
					error(tmp.str(), 0);
					return false;
				}
			}
		}

		for (unsigned int i = 0; i < tokens.size() - 1; i++) {

			//TEST ON SIGN OPERATOR AFTER OTHER OPERATORS OR AFTER OPENING BRACKETS OR AFTER ARGUMENT SEPERATOR: like ...*-2..., ...^-x..., .../-(..., ...+-sin... or like ...(-2..., ...(-x..., ...(-(..., ...(-sin... or like ...,-(...
			rerun = true;
			while (rerun) {
				rerun = false;
				if ((tokens[i]->getType() == 0 || tokens[i]->getType() == 10 || tokens[i]->getType() == 5) && tokens[i + 1]->getType() == 0) {
					std::shared_ptr<OPERATOR> b = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i + 1]);
					if (b->getOPERATOR() == '-' || b->getOPERATOR() == '+') {
						if (tokens.size() >(i + 2)) {
							rerun = true;
							if (b->getOPERATOR() == '-') {
								if (tokens[i + 2]->getType() == 1) {
									std::shared_ptr<SCALAR> c = std::static_pointer_cast<SCALAR, TOKEN>(tokens[i + 2]);
									tokens[i + 2] = std::shared_ptr<SCALAR>(new SCALAR(-c->getSCALAR()));
									c.reset();
								}
								else if (tokens[i + 2]->getType() == 2 || tokens[i + 2]->getType() == 3 || tokens[i + 2]->getType() == 10) {
									tokens.insert(tokens.begin() + i + 2, std::shared_ptr<SCALAR>(new SCALAR(-1)));
									tokens.insert(tokens.begin() + i + 3, std::shared_ptr<OPERATOR>(new OPERATOR('*')));
								}
								else if (tokens[i + 2]->getType() == 0) {
									std::shared_ptr<OPERATOR> c = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i + 2]);
									if (c->getOPERATOR() == '-') {
										tokens.insert(tokens.begin() + i + 2, std::shared_ptr<OPERATOR>(new OPERATOR('+')));
										c.reset();
										tokens.erase(tokens.begin() + i + 3);
									}
									else if (c->getOPERATOR() == '+') {
										tokens.insert(tokens.begin() + i + 2, std::shared_ptr<OPERATOR>(new OPERATOR('-')));
										c.reset();
										tokens.erase(tokens.begin() + i + 3);
									}
									else {
										std::stringstream tmp;
										tmp << "OPERATOR '" << c->getOPERATOR() << "' IS NOT A SIGN OPERATOR / TWO NOT COMPATIBLE TOKENS '" << b->getOPERATOR() << "' AND '" << c->getOPERATOR() << "' MEET!";
										error(tmp.str(), i + 1);
										return false;
									}
								}
								else if (tokens[i + 2]->getType() == 5 || tokens[i + 2]->getType() == -10) {
									std::stringstream tmp;
									tmp << "TWO NOT COMPATIBLE TOKENS '" << b->getOPERATOR() << "' AND '" << (tokens[i + 2]->getType() == 5 ? ',' : ')') << "' MEET!";
									error(tmp.str(), i + 1);
									return false;
								}
								else {
									error("NOT REACHABLE!");
									return false;
								}
							}
							b.reset();
							tokens.erase(tokens.begin() + i + 1);
						}
						else {
							std::stringstream tmp;
							tmp << "EXPRESSION ENDS WITH A SIGN OPERATOR: '" << b->getOPERATOR() << "'!";
							error(tmp.str(), 0);
							return false;
						}
					}
					else {
						std::stringstream tmp;
						tmp << "OPERATOR '" << b->getOPERATOR() << "' IS NOT A SIGN OPERATOR / TWO NOT COMPATIBLE TOKENS '" << (tokens[i]->getType() == 0 ? std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i])->getOPERATOR() : tokens[i]->getType() == 10 ? '(' : ',') << "' AND '" << b->getOPERATOR() << "' MEET!";
						error(tmp.str(), i);
						return false;
					}
				}
			}

			//ADD MULTIPLICATION SIGNS: 
			if ((tokens[i]->getType() == 1 || tokens[i]->getType() == 2 || tokens[i]->getType() == -10) && (tokens[i + 1]->getType() == 1 || tokens[i + 1]->getType() == 2 || tokens[i + 1]->getType() == 3 || tokens[i + 1]->getType() == 10)) {
				tokens.insert(tokens.begin() + i + 1, std::shared_ptr<OPERATOR>(new OPERATOR('*')));
			}

			//IF YOU USE FUNCTIONS WITHOUT BRACKETS: like ...sin2..., ...cos-2...
			rerun = true;
			while (rerun) {
				rerun = false;
				if (tokens[i]->getType() == 3) {
					if (tokens[i + 1]->getType() == 1 || tokens[i + 1]->getType() == 2) {
						tokens.insert(tokens.begin() + i + 1, std::shared_ptr<OPENING_BRACKET>(new OPENING_BRACKET));
						tokens.insert(tokens.begin() + i + 3, std::shared_ptr<CLOSING_BRACKET>(new CLOSING_BRACKET));
					}
					else if (tokens[i + 1]->getType() == 0) {
						std::shared_ptr<OPERATOR> b = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i + 1]);
						if (b->getOPERATOR() == '-' || b->getOPERATOR() == '+') {
							if (tokens.size() > (i + 2)) {
								rerun = true;
								if (b->getOPERATOR() == '-') {
									if (tokens[i + 2]->getType() == 1) {
										std::shared_ptr<SCALAR> c = std::static_pointer_cast<SCALAR, TOKEN>(tokens[i + 2]);
										tokens[i + 2] = std::shared_ptr<SCALAR>(new SCALAR(-c->getSCALAR()));
										c.reset();
									}
									else if (tokens[i + 2]->getType() == 2 || tokens[i + 2]->getType() == 3 || tokens[i + 2]->getType() == 10) {
										tokens.insert(tokens.begin() + i + 2, std::shared_ptr<OPENING_BRACKET>(new OPENING_BRACKET));
										tokens.insert(tokens.begin() + i + 3, std::shared_ptr<SCALAR>(new SCALAR(-1)));
										tokens.insert(tokens.begin() + i + 4, std::shared_ptr<OPERATOR>(new OPERATOR('*')));
										tokens.insert(tokens.begin() + i + 6, std::shared_ptr<CLOSING_BRACKET>(new CLOSING_BRACKET));
									}
									else if (tokens[i + 2]->getType() == 0) {
										std::shared_ptr<OPERATOR> c = std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i + 2]);
										if (c->getOPERATOR() == '-') {
											tokens.insert(tokens.begin() + i + 2, std::shared_ptr<OPERATOR>(new OPERATOR('+')));
											c.reset();
											tokens.erase(tokens.begin() + i + 3);
										}
										else if (c->getOPERATOR() == '+') {
											tokens.insert(tokens.begin() + i + 2, std::shared_ptr<OPERATOR>(new OPERATOR('-')));
											c.reset();
											tokens.erase(tokens.begin() + i + 3);
										}
										else {
											std::stringstream tmp;
											tmp << "OPERATOR '" << c->getOPERATOR() << "' IS NOT A SIGN OPERATOR / TWO NOT COMPATIBLE TOKENS '" << b->getOPERATOR() << "' AND '" << c->getOPERATOR() << "' MEET!";
											error(tmp.str(), i + 1);
											return false;
										}
									}
									else if (tokens[i + 2]->getType() == 5 || tokens[i + 2]->getType() == -10) {
										std::stringstream tmp;
										tmp << "TWO NOT COMPATIBLE TOKENS '" << b->getOPERATOR() << "' AND '" << (tokens[i + 2]->getType() == 5 ? ',' : ')') << "' MEET!";
										error(tmp.str(), i + 1);
										return false;
									}
									else {
										error("NOT REACHABLE!");
										return false;
									}
								}
								b.reset();
								tokens.erase(tokens.begin() + i + 1);
							}
							else {
								std::stringstream tmp;
								tmp << "EXPRESSION ENDS WITH A SIGN OPERATOR: '" << b->getOPERATOR() << "'!";
								error(tmp.str(), 0);
								return false;
							}
						}
						else {
							std::stringstream tmp;
							tmp << "OPERATOR '" << b->getOPERATOR() << "' IS NOT A SIGN OPERATOR / TWO NOT COMPATIBLE TOKENS '" << std::static_pointer_cast<FUNCTION, TOKEN>(tokens[i])->getFUNCTION() << "' AND '" << b->getOPERATOR() << "' MEET!";
							error(tmp.str(), i);
							return false;
						}
					}
				}
			}

		}
		//ADD BRACKETS AT THE END:
		int bracketsToAdd = 0;
		for (unsigned int i = 0; i < tokens.size(); i++) {
			if (tokens[i]->getType() == 10) {
				bracketsToAdd++;
			}
			if (tokens[i]->getType() == -10) {
				bracketsToAdd--;
			}
		}
		for (; bracketsToAdd > 0; bracketsToAdd--) {
			tokens.push_back(std::shared_ptr<CLOSING_BRACKET>(new CLOSING_BRACKET));
		}

		return true;
	}

	bool isLeftAssociative(std::shared_ptr<OPERATOR> op) {
		if (op->getOPERATOR() == '+' || op->getOPERATOR() == '-' || op->getOPERATOR() == '*' || op->getOPERATOR() == '/') {
			return true;
		}
		return false;
	}

	bool isPrecedenceLessOrEqual(std::shared_ptr<OPERATOR> op1, std::shared_ptr<OPERATOR> op2) {
		int precedence1;
		int precedence2;
		if (op1->getOPERATOR() == '+' || op1->getOPERATOR() == '-') {
			precedence1 = 2;
		}
		else if (op1->getOPERATOR() == '*' || op1->getOPERATOR() == '/') {
			precedence1 = 3;
		}
		else if (op1->getOPERATOR() == '^') {
			precedence1 = 4;
		}
		if (op2->getOPERATOR() == '+' || op2->getOPERATOR() == '-') {
			precedence2 = 2;
		}
		else if (op2->getOPERATOR() == '*' || op2->getOPERATOR() == '/') {
			precedence2 = 3;
		}
		else if (op2->getOPERATOR() == '^') {
			precedence2 = 4;
		}
		if (precedence1 <= precedence2) {
			return true;
		}
		return false;
	}

	bool shuntingYard() {
		std::stack<std::shared_ptr<TOKEN>> operatorStack;
		for (unsigned int i = 0; i < tokens.size(); i++) {
			if (tokens[i]->getType() == 1) {
				rpn.push_back(std::static_pointer_cast<SCALAR, TOKEN>(tokens[i]));
				continue;
			}

			if (tokens[i]->getType() == 2) {
				rpn.push_back(std::static_pointer_cast<VARIABLE, TOKEN>(tokens[i]));
				continue;
			}

			if (tokens[i]->getType() == 3) {
				operatorStack.push(std::static_pointer_cast<FUNCTION, TOKEN>(tokens[i]));
				continue;
			}

			if (tokens[i]->getType() == 5) {
				while (true) {
					if (operatorStack.empty()) {
						error("THERE IS AN INCORRECTLY PLACED ARGUMENT SEPARATOR OR THE CLOSING BRACKET ISN'T PRECEDED BY AN OPENING BRACKET!", i);
						return false;
					}
					if (operatorStack.top()->getType() == 10) {
						break;
					}
					else {
						if (operatorStack.top()->getType() == 0) {
							rpn.push_back(std::static_pointer_cast<OPERATOR, TOKEN>(operatorStack.top()));
						}
						else if (operatorStack.top()->getType() == 1) {
							rpn.push_back(std::static_pointer_cast<SCALAR, TOKEN>(operatorStack.top()));
						}
						else if (operatorStack.top()->getType() == 2) {
							rpn.push_back(std::static_pointer_cast<VARIABLE, TOKEN>(operatorStack.top()));
						}
						else if (operatorStack.top()->getType() == 3) {
							rpn.push_back(std::static_pointer_cast<FUNCTION, TOKEN>(operatorStack.top()));
						}
						else {
							error("NOT REACHABLE!");
							return false;
						}
						operatorStack.pop();
					}
				}
				continue;
			}

			if (tokens[i]->getType() == 0) {
				while (true) {
					if (operatorStack.empty()) {
						break;
					}
					if (operatorStack.top()->getType() == 0) {
						if (isLeftAssociative(std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i])) && isPrecedenceLessOrEqual(std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i]), std::static_pointer_cast<OPERATOR, TOKEN>(operatorStack.top()))) {
							rpn.push_back(std::static_pointer_cast<OPERATOR, TOKEN>(operatorStack.top()));
							operatorStack.pop();
						}
						else {
							break;
						}						
					}
					else {
						break;
					}
				}
				operatorStack.push(std::static_pointer_cast<OPERATOR, TOKEN>(tokens[i]));
				continue;
			}

			if (tokens[i]->getType() == 10) {
				operatorStack.push(std::static_pointer_cast<OPENING_BRACKET, TOKEN>(tokens[i]));
				continue;
			}

			if (tokens[i]->getType() == -10) {
				while (true) {
					if (operatorStack.empty()) {
						error("THE CLOSING BRACKET ISN'T PRECEDED BY AN OPENING BRACKET!", i);
						return false;
					}
					if (operatorStack.top()->getType() == 10) {
						operatorStack.pop();
						break;
					}
					else {
						rpn.push_back(operatorStack.top());
						operatorStack.pop();
					}
				}
				if (!operatorStack.empty()) {
					if (operatorStack.top()->getType() == 3) {
						rpn.push_back(std::static_pointer_cast<FUNCTION, TOKEN>(operatorStack.top()));
						operatorStack.pop();
					}
				}
				continue;
			}
			error("NOT REACHABLE!");
			return false;
		}
		while (true) {
			if (operatorStack.empty()) {
				break;
			}
			if (operatorStack.top()->getType() == 10) {
				error("THERE ARE MORE OPENING THAN CLOSING BRACKETS!");
				return false;
			}
			else {
				rpn.push_back(operatorStack.top());
				operatorStack.pop();
			}
		}
		return true;
	}

public:
	double calculate(std::vector<double> x) {
		double totalresult = NAN;
		std::stack<std::shared_ptr<TOKEN>> evaluationStack;
		std::vector<std::string> varnames;
		for (unsigned int i = 0; i < rpn.size(); i++) {
			if (rpn[i]->getType() == 1) {
				evaluationStack.push(std::static_pointer_cast<SCALAR, TOKEN>(rpn[i]));
				continue;
			}

			if (rpn[i]->getType() == 2) {
				unsigned int pos;
				std::string varname = std::static_pointer_cast<VARIABLE, TOKEN>(rpn[i])->getVARIABLE();
				bool registered = false;
				for (unsigned int i = 0; i < varnames.size(); i++) {
					if (varname == varnames[i]) {
						registered = true;
						pos = i;
						break;
					}
				}
				if (!registered) {
					varnames.push_back(varname);
					pos = varnames.size() - 1;
				}

				if (varnames.size() > x.size()) {
					std::stringstream tmp;
					tmp << "THERE ARE MORE VARIABLES THAN INPUT VALUES FOR THIS FUNCTION, VARAIBLES ARE: ";
					for (unsigned int i = 0; i < varnames.size(); i++) {
						tmp << "'" << varnames[i] << "', ";
					}
					tmp << "INPUT VALUES ARE: ";
					for (unsigned int i = 0; i < x.size(); i++) {
						tmp << "'" << x[i] << "', ";
					}
					tmp << "!";
					error(tmp.str());
					//return false;
				}
				evaluationStack.push(std::shared_ptr<SCALAR>(new SCALAR(x[pos])));
				continue;
			}

			if (rpn[i]->getType() == 0) {
				double result = 0;

				char op = std::static_pointer_cast<OPERATOR, TOKEN>(rpn[i])->getOPERATOR();
				
				//TEST ON TOKEN TO OPERATE:
				if (evaluationStack.empty()) {
					std::stringstream tmp;
					tmp << "OPERATOR '" << op << "' HAS NO TOKEN TO OPERATE!";
					error(tmp.str(), i);
					//return false;
				}
				//TEST ON SCALAR TO OPERATE:
				else if (evaluationStack.top()->getType() != 1) {
					std::stringstream tmp;
					tmp << "OPERATOR '" << op << "' HAS NO SCALAR TO OPERATE!";
					error(tmp.str(), i);
					//return false;
				}

				std::shared_ptr<SCALAR> right = std::static_pointer_cast<SCALAR, TOKEN>(evaluationStack.top());
				evaluationStack.pop();

				//TEST ON TOKENS TO OPERATE:
				if (evaluationStack.empty()) {
					std::stringstream tmp;
					tmp << "OPERATOR '" << op << "' NEEDS 2 TOKENS TO OPERATE!";
					error(tmp.str(), i);
					//return false;
				}
				//TEST ON SCALARS TO OPERATE:
				else if (evaluationStack.top()->getType() != 1) {
					std::stringstream tmp;
					tmp << "OPERATOR '" << op << "' NEEDS 2 SCALARS TO OPERATE!";
					error(tmp.str(), i);
					//return false;
				}

				std::shared_ptr<SCALAR> left = std::static_pointer_cast<SCALAR, TOKEN>(evaluationStack.top());
				evaluationStack.pop();

				if (op == '+') {
					result = left->getSCALAR() + right->getSCALAR();
				}
				else if (op == '-') {
					result = left->getSCALAR() - right->getSCALAR();
				}
				else if (op == '*') {
					result = left->getSCALAR() * right->getSCALAR();
				}
				else if (op == '/') {
					result = left->getSCALAR() / right->getSCALAR();
				}
				else if (op == '^') {
					result = pow(left->getSCALAR(), right->getSCALAR());
				}
				else {
					error("NOT REACHABLE!");
					//return false;
				}
				evaluationStack.push(std::shared_ptr<SCALAR>(new SCALAR(result)));
				continue;
			}

			if (rpn[i]->getType() == 3) {
				double result = 0;

				std::string f = std::static_pointer_cast<FUNCTION, TOKEN>(rpn[i])->getFUNCTION();

				//TEST ON INPUT VALUE:
				if (evaluationStack.empty()) {
					std::stringstream tmp;
					tmp << "FUNCTION '" << f << "' HAS NO INPUT VALUE!";
					error(tmp.str(), i);
					//return false;
				}
				//TEST ON SCALAR INPUT VALUE:
				else if (evaluationStack.top()->getType() != 1) {
					std::stringstream tmp;
					tmp << "FUNCTION '" << f << "' HAS NO SCALAR INPUT VALUE!";
					error(tmp.str(), i);
					//return false;
				}

				std::shared_ptr<SCALAR> right = std::static_pointer_cast<SCALAR, TOKEN>(evaluationStack.top());
				evaluationStack.pop();

				if (f == "exp") {
					result = std::exp(right->getSCALAR());
				}
				else if (f == "ln") {
					result = std::log(right->getSCALAR());
				}
				else if (f == "log10") {
					result = std::log10(right->getSCALAR());
				}
				else if (f == "log2") {
					result = std::log2(right->getSCALAR());
				}
				else if (f == "sin") {
					result = std::sin(right->getSCALAR());
				}
				else if (f == "cos") {
					result = std::cos(right->getSCALAR());
				}
				else if (f == "tan") {
					result = std::tan(right->getSCALAR());
				}
				else if (f == "arcsin") {
					result = std::asin(right->getSCALAR());
				}
				else if (f == "arccos") {
					result = std::acos(right->getSCALAR());
				}
				else if (f == "arctan") {
					result = std::atan(right->getSCALAR());
				}
				else if (f == "sinh") {
					result = std::sinh(right->getSCALAR());
				}
				else if (f == "cosh") {
					result = std::cosh(right->getSCALAR());
				}
				else if (f == "tanh") {
					result = std::tanh(right->getSCALAR());
				}
				else if (f == "arsinh") {
					result = std::asinh(right->getSCALAR());
				}
				else if (f == "arcosh") {
					result = std::acosh(right->getSCALAR());
				}
				else if (f == "artanh") {
					result = std::atanh(right->getSCALAR());
				}
				else if (f == "abs") {
					result = std::abs(right->getSCALAR());
				}
				else if (f == "min" || f == "max") {

					//TEST ON ENOUGH INPUT VALUES:
					if (evaluationStack.empty()) {
						std::stringstream tmp;
						tmp << "FUNCTION '" << f << "' HAS NOT ENOUGH INPUT VALUES (NEEDS 2)!";
						error(tmp.str(), i);
						//return false;
					}
					//TEST ON ENOUGH SCALAR INPUT VALUES:
					else if (evaluationStack.top()->getType() != 1) {
						std::stringstream tmp;
						tmp << "FUNCTION '" << f << "' HAS NOT ENOUGH SCALAR INPUT VALUES (NEEDS 2)!";
						error(tmp.str(), i);
						//return false;
					}

					std::shared_ptr<SCALAR> left = std::static_pointer_cast<SCALAR, TOKEN>(evaluationStack.top());
					evaluationStack.pop();

					if (f == "min") {
						result = std::min<double>(left->getSCALAR(), right->getSCALAR());
					}
					else if (f == "max") {
						result = std::max<double>(left->getSCALAR(), right->getSCALAR());
					}
					else {
						error("NOT REACHABLE!");
						//return false;
					}
				}
				else {
					std::stringstream tmp;
					tmp << "FUNCTION '" << f << "' IS NOT SUPPORTED!";
					error(tmp.str(), i);
					//return false;
				}

				evaluationStack.push(std::shared_ptr<SCALAR>(new SCALAR(result)));
				continue;
			}

			error("NOT REACHABLE!");
			//return false;
		}

		//TEST TOTAL RESULT ON SCALAR:
		if (evaluationStack.empty()) {
			error("NO RESULT!");
			//return false;
		}
		else if (evaluationStack.top()->getType() != 1) {
			error("RESULT IS NOT SCALAR!");
			//return false;
		}
		
		totalresult = std::static_pointer_cast<SCALAR, TOKEN>(evaluationStack.top())->getSCALAR();
		if (totalresult != totalresult) {
			std::stringstream tmp;
			tmp << "RESULT IS: '" << totalresult << "', FOR FUNCTION: '";
			for (unsigned int i = 0; i < tokens.size(); i++) {
				tmp << *(tokens[i]) << " ";
			}
			tmp << "'!";
			error(tmp.str());
		}
		return totalresult;
	}

	~INTERPRETER() {
		for (unsigned int i = 0; i < tokens.size(); i++) {
			tokens[i].reset();
		}
		for (unsigned int i = 0; i < rpn.size(); i++) {
			//delete rpn[i]; da in rpn nur zeiger von tokens sind;
		}
	}
};