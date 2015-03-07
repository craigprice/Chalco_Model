

/**
 ClassicalSpin3D.h
 Purpose: Spin object for the simulation of the lattice.
 
 @author Craig Price
 @version 1.0 2015/03/04
 */

#ifndef CLASSICALSPIN3D_H
#define CLASSICALSPIN3D_H
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

class ClassicalSpin3D{
public:
    
    /**
     This creates a spin object with its spin components = 0.
     */
    ClassicalSpin3D();

    /**
     This creates a spin object with its direction in the x_, y_, and z_
     direction.
     
     @param x_ the x component of the spin.
     @param y_ the y component of the spin.
     @param z_ the z component of the spin.
     */
    ClassicalSpin3D(double x_, double y_, double z_);
    
    /**
     This prints out to cout the value of x and y and z.
     */
    void print() const;
    
    /**
     This makes sure that the value stored in the x and y and z values are
     appropriate for a 3D classical spin with magnitude of 1.
     */
    void checkSpin();
    
    
    /**
     Takes the dot product between this spin and the one passed to it.
     
     @param spin the other vector to take the dot product with.
     @return the value of the dot product.
     */
    double dotProd(const ClassicalSpin3D& spin) const;
    
    /**
     Sets this spin to be oriented randomly in the unit spherical surface.
     */
    void setRandomOrientation();
    
    /**
     flips this spin from it's current orientation to a random new orientation
     Where the value of the angular rotation is at most range.
     
     @param range the polar angle theta that defines the polar spherical cap
     within which lies the possible orientations for the spin to point in.
     */
    void flip(double range);
    
    /**
     Rotate the spin.
     
     @param phi
     @param theta
     */
    void specifyRotation(double phi, double theta);
    
    /**
     Returns this spin to the value that it had previously.
     */
    void reset();
    
    /**
     Sets this spin's information to 0;
     */
    void clear();
    
    /**
     This function rotates the spin about an axis (specified by 'theta' and
     'phi') by an angle of 'angle'.
     
     @param theta
     @param phi
     @param angle
     */
    void rotate(double theta, double phi, double angle);
    
    unsigned int xPos; //position of the spin along the "a" axis
    unsigned int yPos; //position of the spin along the "b" axis
    unsigned int sPos; //sublatice that the spin belongs within
    unsigned int zPos; //position of the spin along the "c" axis
    double x; //x component of the 3D unit vector that is the spin
    double y; //y component of the 3D unit vector that is the spin
    double z; //z component of the 3D unit vector that is the spin
    double xOld; // the previous x component of the spin. (position before flip)
    double yOld; // the previous y component of the spin. (position before flip)
    double zOld; // the previous z component of the spin. (position before flip)
    
};

inline double ClassicalSpin3D::dotProd(const ClassicalSpin3D& spin) const{
    return (x * spin.x +
            y * spin.y +
            z * spin.z);
}


#endif
