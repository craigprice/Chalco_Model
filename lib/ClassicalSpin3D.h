
#ifndef CLASSICALSPIN3D_H
#define CLASSICALSPIN3D_H
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

class ClassicalSpin3D{
public:
    ClassicalSpin3D();

    
    /**
     *This creates a spin object with it's direction in the x_, y_, and z_
     *direction.
     */
    ClassicalSpin3D(double x_, double y_, double z_);
    
    /**
     *This prints out to cout the value of x and y and z.
     */
    void print() const;
    
    /**
     *This makes sure that the value stored in the x and y and z values are 
     *appropriate for a 3D classical spin with magnitude of 1.
     */
    void checkSpin();
    
    /**
     *Returns the value of the x component of the spin.
     */
    double getX() const;
    
    /**
     *Returns the value of the y component of the spin.
     */
    double getY() const;
    
    /**
     *Returns the value of the z component of the spin.
     */
    double getZ() const;
    
    void setPosX(unsigned int x_);
    void setPosY(unsigned int y_);
    void setPosS(unsigned int s_);
    void setPosZ(unsigned int z_);
    unsigned int getPosX() const;
    unsigned int getPosY() const;
    unsigned int getPosS() const;
    unsigned int getPosZ() const;
    void setX(double x_);
    void setY(double y_);
    void setZ(double z_);
    
    void doC2Rotation001();
    
    /**
     *Takes the dot product between this spin and the one passed to it.
     */
    double dotProd(const ClassicalSpin3D& spin) const;
    
    /**
     *Sets this spin to be oriented randomly in the unit spherical surface.
     */
    void setRandomOrientation();
    
    /**
     *flips this spin from it's current orientation to a random new orientation
     *Where the value of the angular rotation is at most range.
     */
    void flip(double range);
    void specifyRotation(double phi, double theta);
    
    /**
     *Returns this spin to the value that it had previously.
     */
    void reset();
    
    /**
     *Sets this spin's information to 0;
     */
    void clear();
private:
    //ClassicalSpin3D(const ClassicalSpin3D& s);
    //ClassicalSpin3D & operator = (const ClassicalSpin3D& c);
    
    /**
     *This function rotates the spin about an axis (specified by 'theta' and
     *'phi') by an angle of 'angle'.
     */
    void rotate(double theta, double phi, double angle);
    unsigned int xPos;
    unsigned int yPos;
    unsigned int sPos;
    unsigned int zPos;
    double x;
    double y;
    double z;
    double xOld;
    double yOld;
    double zOld;
};
inline double ClassicalSpin3D::getX() const{
    return x;
}
inline double ClassicalSpin3D::getY() const{
    return y;
}
inline double ClassicalSpin3D::getZ() const{
    return z;
}
inline void ClassicalSpin3D::setX(double x_) {
    x = x_;
}
inline void ClassicalSpin3D::setY(double y_) {
    y = y_;
}
inline void ClassicalSpin3D::setZ(double z_) {
    z = z_;
}
inline unsigned int ClassicalSpin3D::getPosX() const{
    return xPos;
}
inline unsigned int ClassicalSpin3D::getPosY() const{
    return yPos;
}
inline unsigned int ClassicalSpin3D::getPosS() const{
    return sPos;
}
inline unsigned int ClassicalSpin3D::getPosZ() const{
    return zPos;
}
inline void ClassicalSpin3D::setPosX(unsigned int x_) {
    xPos = x_;
}
inline void ClassicalSpin3D::setPosY(unsigned int y_) {
    yPos = y_;
}
inline void ClassicalSpin3D::setPosS(unsigned int s_) {
    sPos = s_;
}
inline void ClassicalSpin3D::setPosZ(unsigned int z_) {
    zPos = z_;
}
inline double ClassicalSpin3D::dotProd(const ClassicalSpin3D& spin) const{
    return (x * spin.getX() + 
            y * spin.getY() + 
            z * spin.getZ());
}


#endif
