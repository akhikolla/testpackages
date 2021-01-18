#ifndef POSITION_H
#define POSITION_H
class Position {
	public:
	double x;
	double y;

	Position(double x, double y);
	Position();
	Position(const Position& rhs);


	Position& operator+=(const Position& rhs);
	Position& operator-=(const Position& rhs);
	Position& operator*=(const double rhs);
};

Position operator+(Position lhs, const Position& rhs);
Position operator-(Position lhs, const Position& rhs);
Position operator*(const double lhs, Position& rhs);
#endif
