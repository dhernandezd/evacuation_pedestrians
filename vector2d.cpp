#include "vector2d.h"

/*!
 * \details Sets each element to zero.
 */
void Vec2D::setZero()
{
    X = 0.0;
    Y = 0.0;
}

/*!
 * \details Calculates the dot product of two vectors.
 * NB: this is a STATIC function!
 * \param[in] a     the first vector
 * \param[in] b     the second vector 
 * \return          the resulting scalar
 */
double Vec2D::dot(const Vec2D& a, const Vec2D& b)
{
    return a.X * b.X + a.Y * b.Y;
}

void Vec2D::normalise()
{
    double length2 = this->getLengthSquared();
    if (length2 == 0)
    {
        std::cerr<<"Normalizing a vector of length 0\n";
    }
    *this /= std::sqrt(length2);
}

void Vec2D::setLength(double length)
{
    this->normalise();
    *this *= length;
}

double Vec2D::cross(const Vec2D& a, const Vec2D& b)
{
    return  a.X * b.Y - a.Y * b.X;
}

double Vec2D::getDistance(const Vec2D& a, const Vec2D& b)
{
    return std::sqrt(getDistanceSquared(a, b));
}

double Vec2D::getLengthSquared() const
{
    return (X * X + Y * Y);
}

double Vec2D::getLength() const
{
    return std::sqrt(getLengthSquared());
}

double Vec2D::getComponent(const int index) const
{
    switch (index)
    {
        case 0:
            return X;
        case 1:
            return Y;
        default:
            std::cerr<<"[Vector::getComponent] Index is too high for a 2D vector (should be 0-1).\n";
            return 0;
    }
}

void Vec2D::setComponent(const int index, const double val)
{
    switch (index)
    {
        case 0:
            X = val;
            break;
        case 1:
            Y = val;
            break;
        default:
            std::cerr<<"[Vector::setComponent] Index is too high for a 2D vector (should be 0-1).\n";
    }
}

bool Vec2D::isEqualTo(const Vec2D& other, const double tol) const
{
    if ((Vec2D::getLengthSquared(*this - other)) <= tol * tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Vec2D Vec2D::getUnitVector(const Vec2D& a)
{
    double Length2 = a.getLengthSquared();
    if (Length2 != 0.0)
        return a / std::sqrt(Length2);
    else
        return Vec2D(0, 0);
}

double Vec2D::getLength(const Vec2D& a)
{
    return a.getLength();
}

/*!
 * \details Adds all elements of the vector to an output stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in] os    the output stream, 
 * \param[in] a     The vector of interest
 * \return          the output stream with vector elements added
 */
std::ostream& operator<<(std::ostream& os, const Vec2D& a)
{
    os << a.X << ' ' << a.Y;
    return os;
}

/*!
 * \details Reads all elements of a given vector from an input stream.
 * NB: this is a global function and a friend of the Vec3D class!
 * \param[in,out] is    the input stream
 * \param[in,out] a     the vector to be read in
 * \return              the input stream from which the vector elements were read
 */
std::istream& operator>>(std::istream& is, Vec2D& a)
{
    //TW: clearing the stream avoids the nasty problem that the failbit is set to true if numbers below DBL_MIN=1e-308 are read.
    is >> a.X; is.clear();
    is >> a.Y;
    return is;
}
