#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <cmath>
#include <sstream>
#include <iostream>
#include <cstdlib>

/*!
 * \class Vec2D
 * \brief here, it is included the necessary functions to operate with vectors
 */
class Vec2D
{
public:
    
    /*!
     * \brief the vector components
     */
    /*
     * \todo: Make these private.
     * \todo what is the idea of this constructor?
     * These should be private so we can implement things like a cvec etc.
     * Use getters / setters.
     */
//    Vec2D(int i);
    
    // private:
    double X, Y;
    
    /*!
     * \brief constructor
     */
    Vec2D()
    { setZero(); }
    
    /*!
     * \brief Alternative constructor, taking the three elements as arguments
     * \details Alternative constructor, lets you define all three elements.
     * \param[in] x     the x-component
     * \param[in] y     the y-component
     */
    Vec2D(const double x, const double y)
    {
        X = x;
        Y = y;
    }
    
    /*!
     * \brief Sets all elements to zero
     */
    void setZero();
    
    
    /*!
     * \brief Checks if ALL elements are zero
     */
    bool isZero() const
    { return X == 0.0 && Y == 0.0; }
    
    
    /*!
     * \brief Adds another vector
     * \details Adds vector to itself
     * \param[in] a     vector to be added
     * \return          resulting 3D vector
     */
    Vec2D operator+(const Vec2D& a) const
    {
        return Vec2D(X + a.X, Y + a.Y);
    }
    
    /*!
     * \brief Binary vector subtraction
     * \details Subtracts a vector from another vector
     * \param[in] a     vector to be subtracted
     * \return          resulting vector
     */
    inline Vec2D operator-(const Vec2D a) const
    {
        return Vec2D(X - a.X, Y - a.Y);
    };


    /*!
     * \brief Multiplies by a scalar
     * \details Multiplies each element with a scalar
     * \param[in] a     the scalar to be multiplied with
     * \return          the resulting vector
     */
    Vec2D operator*(const double a) const
    {
        return Vec2D(X * a, Y * a);
    }
    
    /*!
     * \brief Divides by a scalar
     * \details Divides each element by a scalar
     * \param[in] a     the scalar to be divided by
     * \return          resulting vector
     */
    Vec2D operator/(double a) const {
        return Vec2D(X / a, Y / a);
    }
    
    /*!
     * \brief Adds another vector
     * \details Adds a vector to itself
     * \param[in] a     vector to be added
     * \return          (reference to) itself, i.e. resulting vector
     */
    Vec2D& operator+=(const Vec2D& a)
    {
        X += a.X;
        Y += a.Y;
        return *this;
    }

    /*!
     * \brief Checks if all coordinates satisfy this>=a
     */
    bool operator>=(const Vec2D& a) const {
        return X>=a.X && Y>=a.Y;
    }

    bool operator<(const Vec2D& a) const {
        return X<a.X && Y<a.Y;
    }

    /*!
     * \brief Subtracts another vector
     * \details Subtracts a vector from itself
     * \param[in] a     vector to be subtracted
     * \return          (reference to) itself, i.e. resulting vector
     */
    Vec2D& operator-=(const Vec2D& a)
    {
        X -= a.X;
        Y -= a.Y;
        return *this;
    }

    /*!
     * \brief Multiplies by a scalar
    * \details Multiplies each element by a scalar
    * \param[in] a     scalar to be multiplied by
    * \return          (reference to) itself, i.e. resulting vector
    */
    Vec2D& operator*=(double a) {
        X *= a;
        Y *= a;
        return *this;
    }

    /*!
     * \brief Divides by a scalar
     * \details Divides each element by a scalar
     * \param[in] a     scalar to be divided by
     * \return          (reference to) itself, i.e. resulting vector
     */
    Vec2D& operator/=(const double a)
    {
        X /= a;
        Y /= a;
        return *this;
    }

    /*!
     * \brief Calculates the dot product of two Vec2D: \f$ a \cdot b\f$
     */
    static double dot(const Vec2D& a, const Vec2D& b);
    
    
    /*!
     * \brief Makes this Vec2D unit length
     */
    void normalise();
    
    /*!
     * \brief Make this Vec2D a certain length 
     */
    void setLength(double length);
    
    /*!
     * \brief Calculates the cross product of two Vec2D:  \f$ a \times b\f$. 
     *  Note: the result if a scalar denoting the z-component only since we are dealing only with 2D 
     *  systems
     */
    static double cross(const Vec2D& a, const Vec2D& b);
    
    /*!
     * \brief Calculates the distance between two Vec2D: \f$ \sqrt{\left(a-b\right) \cdot \left(a-b\right)} \f$
     * \details Calculates the square of the distance (i.e. the length of the difference)
     * between two vectors.
     * NB: this is a STATIC function!
     * \param[in] a     the first vector
     * \param[in] b     the second vector
     * \return          the square of the distance between the two arguments.
     */
    static double getDistance(const Vec2D& a, const Vec2D& b);
    
    /*!
     * \brief Calculates the squared distance between two Vec2D: \f$ \left(a-b\right) \cdot \left(a-b\right) \f$
     */
    static double getDistanceSquared(const Vec2D& a, const Vec2D& b) {
        const double X = a.X-b.X;
        const double Y = a.Y-b.Y;
        return (X * X + Y * Y);
        //return getLengthSquared(a - b);
    }

    
    /*!
     * \brief Calculates the length of a Vec2D: \f$ \sqrt{a\cdot a} \f$
     */
    static double getLength(const Vec2D& a);
    
    /*!
     * \brief Calculates the squared length of a Vec2D: \f$ a\cdot a \f$
     * \details Calculates the square of the length of a given vector.
     * NB: this is a STATIC function!
     * \param[in] a     the vector.
     * \return          the square of the length of the argument.
     */
    static double getLengthSquared(const Vec2D& a)
    {
        return (a.X * a.X + a.Y * a.Y);
    }
    
    /*!
     * \brief Calculates the length of this Vec2D: \f$ \sqrt{a\cdot a} \f$
     */
    double getLength() const;
    
    /*!
     * \brief Calculates the squared length of this Vec2D: \f$ a\cdot a \f$
     */
    double getLengthSquared() const;
    
    /*!
     * \brief Returns the requested component of this Vec2D
     */
    double getComponent(int index) const;
    
    /*!
     * \brief Sets the requested component of this Vec2D to the requested value
     */
    void setComponent(int index, double val);
    
    /*!
     * \brief RW reference to X
     */
    inline double& x()
    { return X; }
    
    /*!
     * \brief RO reference to X
     */
    inline double x() const
    { return X; }
    
    /*!
     * \brief RW reference to Y
     */
    inline double& y()
    { return Y; }
    
    /*!
     * \brief RO reference to Y
     */
    inline double y() const
    { return Y; }
    
    
    inline void setX(double x)
    { X = x; }
    
    inline void setY(double y)
    { Y = y; }
    
    inline double getX()
    { return X; }
    
    inline double getY()
    { return Y; }
    
    
    inline void set(double x, double y)
    {
        X = x;
        Y = y;
    }
    
    /*!
     * \brief Checks if the length this Vec2D is equal the length of other with a certain tolerance
     */
    bool isEqualTo(const Vec2D& other, double tol) const;
    
    /*!
     * \brief Returns a unit Vec2D based on a.
     */
    static Vec2D getUnitVector(const Vec2D& a);
    
    /*!
     * \brief Adds elements to an output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Vec2D& a);
    
    /*!
     * \brief Adds elements to an input stream
     */
    friend std::istream& operator>>(std::istream& is, Vec2D& a);
    
    /*!
     * \brief Reverts the direction of a vector
     */
    inline friend Vec2D operator-(const Vec2D& a) {
        return Vec2D(-a.X, -a.Y);
    }
    
    /*!
     * \brief Multiplies all elements by a scalar
    * \details Multiplies each element of a given vector (b) by a given scalar (a).
    * NB: this is a global function and a friend of the Vec2D class. Gets called when
    * a scalar multiplication of the form (double) * (Vec2D) is performed.
    * \param[in] a     the scalar
    * \param[in] b     the vector
    * \return          the resulting vector
    */
    friend Vec2D operator*(double a, const Vec2D& b) {
        return Vec2D(b.X * a, b.Y * a);
    }



};

#endif