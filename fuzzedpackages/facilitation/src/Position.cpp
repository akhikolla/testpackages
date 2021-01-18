#include"Position.h"

Position operator+(Position lhs, const Position& rhs)
{
	lhs += rhs;
	return lhs;
}

Position operator-(Position lhs, const Position& rhs)
{
	lhs -= rhs;
	return lhs;
}

Position operator*(const double lhs, Position& rhs)
{
	rhs *= lhs;
	return rhs;
}

Position& Position::operator+=(const Position& rhs){
	x += rhs.x;
	y += rhs.y;
	return *this;
}
Position& Position::operator-=(const Position& rhs){
	x -= rhs.x;
	y -= rhs.y;
	return *this;
}
Position& Position::operator*=(const double rhs){
	x *= rhs;
	y *= rhs;
	return *this;
}

Position::Position(double x, double y): x(x),y(y){}
Position::Position(){}
Position::Position(const Position& rhs): x(rhs.x),y(rhs.y){}
