#pragma once

#include <string>
#include <iostream>

class TOKEN {
	int type;
public:
	TOKEN(int type) : type(type) {}
	int getType() { return type; }

	virtual std::ostream & output(std::ostream & os) = 0;
	friend std::ostream & operator<<(std::ostream&, const TOKEN&);

	virtual ~TOKEN() {}
};

std::ostream & operator<<(std::ostream & os, TOKEN & obj) {
	return obj.output(os);
}

/*
TOKEN:			  |	TYPE:
------------------|-------
OPERATOR		  |	0
SCALAR			  |	1
VARIABLE		  |	2
FUNCTION		  |	3
ARGUMENT_SEPARATOR|	5
OPENING_BRACKET	  |	10
CLOSING_BRACKET	  |	-10
*/

class OPERATOR : public TOKEN {
	char o;
public:
	OPERATOR(char o) : o(o), TOKEN(0) {}
	char getOPERATOR() { return o; }

	std::ostream & output(std::ostream & os) { return os << o; }

	virtual ~OPERATOR() {}
};

class SCALAR : public TOKEN {
	double s;
public:
	SCALAR(double s) : s(s), TOKEN(1) {}
	double getSCALAR() { return s; }

	std::ostream & output(std::ostream & os) { return os << s; }

	virtual ~SCALAR() {}
};

class VARIABLE : public TOKEN {
	std::string v;
public:
	VARIABLE(std::string v) : v(v), TOKEN(2) {}
	std::string getVARIABLE() { return v; }

	std::ostream & output(std::ostream & os) { return os << v; }

	virtual ~VARIABLE() {}
};

class FUNCTION : public TOKEN {
	std::string f;
public:
	FUNCTION(std::string f) : f(f), TOKEN(3) {}
	std::string getFUNCTION() { return f; }

	std::ostream & output(std::ostream & os) { return os << f; }

	virtual ~FUNCTION() {}
};

class ARGUMENT_SEPARATOR : public TOKEN {
public:
	ARGUMENT_SEPARATOR() : TOKEN(5) {}

	std::ostream & output(std::ostream & os) { return os << ","; }

	virtual ~ARGUMENT_SEPARATOR() {}
};

class OPENING_BRACKET : public TOKEN {
public:
	OPENING_BRACKET() : TOKEN(10) {}

	std::ostream & output(std::ostream & os) { return os << "("; }

	virtual ~OPENING_BRACKET() {}
};

class CLOSING_BRACKET : public TOKEN {
public:
	CLOSING_BRACKET() : TOKEN(-10) {}

	std::ostream & output(std::ostream & os) { return os << ")"; }

	virtual ~CLOSING_BRACKET() {}
};
