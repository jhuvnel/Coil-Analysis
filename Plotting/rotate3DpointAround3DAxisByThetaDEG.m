
function qq=rotate3DpointAround3DAxisByThetaDEG(p,rr,thetadegrees)
%rotate a point p (or tip of a vector whose tail is at origin) in 3D 
%around a rotation axis rr by angle thetadegrees (in DEGREES) 
%NOTE: To use this routine to rotate a vector that is in LRZ coordinates
% around the pitch axis (the Y axis in XYZ coordinates), you should run the
% routine using a rotation axis of [-1 1 0] (which is the +Y axis as viewed
% by the LRZ coordinate frame, which "thinks" it is really xyz.
   qq = [0;0;0];
   r=rr/norm(rr); %ensure r is normalized
   costheta = cosd(thetadegrees);
   sintheta = sind(thetadegrees);

   qq(1) = qq(1) + (costheta + (1 - costheta) * r(1) * r(1)) * p(1);
   qq(1) = qq(1) + ((1 - costheta) * r(1) * r(2) - r(3) * sintheta) * p(2);
   qq(1) = qq(1) + ((1 - costheta) * r(1) * r(3) + r(2) * sintheta) * p(3);

   qq(2) = qq(2) + ((1 - costheta) * r(1) * r(2) + r(3) * sintheta) * p(1);
   qq(2) = qq(2) + (costheta + (1 - costheta) * r(2) * r(2)) * p(2);
   qq(2) = qq(2) + ((1 - costheta) * r(2) * r(3) - r(1) * sintheta) * p(3);

   qq(3) = qq(3) + ((1 - costheta) * r(1) * r(3) - r(2) * sintheta) * p(1);
   qq(3) = qq(3) + ((1 - costheta) * r(2) * r(3) + r(1) * sintheta) * p(2);
   qq(3) = qq(3) + (costheta + (1 - costheta) * r(3) * r(3)) * p(3);
end

%{

Details of this approach are detailed here:
 http://paulbourke.net/geometry/rotate/
and in the floowing 3 realizations in C and Python. 
You can also read about them in Haslwanter, Thomas Vision Research 3D
kinematics paper ~1997 and in Americo Migliaccio's PhD thesis.
http://www.engr.uvic.ca/~mech410/lectures/4_2_RotateArbi.pdf

%}

%{
===================
==========C code: http://paulbourke.net/geometry/rotate/example.c
====================
XYZ RotatePointAboutLine(XYZ p,double theta,XYZ p1,XYZ p2)
{
   XYZ u,q1,q2;
   double d;

   /* Step 1 */
   q1.x = p.x - p1.x;
   q1.y = p.y - p1.y;
   q1.z = p.z - p1.z;

   u.x = p2.x - p1.x;
   u.y = p2.y - p1.y;
   u.z = p2.z - p1.z;
   Normalise(&u);
   d = sqrt(u.y*u.y + u.z*u.z);

   /* Step 2 */
   if (d != 0) {
      q2.x = q1.x;
      q2.y = q1.y * u.z / d - q1.z * u.y / d;
      q2.z = q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   /* Step 3 */
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;

   /* Step 4 */
   q2.x = q1.x * cos(theta) - q1.y * sin(theta);
   q2.y = q1.x * sin(theta) + q1.y * cos(theta);
   q2.z = q1.z;

   /* Inverse of step 3 */
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   /* Inverse of step 2 */
   if (d != 0) {
      q2.x =   q1.x;
      q2.y =   q1.y * u.z / d + q1.z * u.y / d;
      q2.z = - q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   /* Inverse of step 1 */
   q1.x = q2.x + p1.x;
   q1.y = q2.y + p1.y;
   q1.z = q2.z + p1.z;
   return(q1);
}
%}
%{
=====================================
======A more text-efficient realization, using quaterions:
=====================================

typedef struct {
   double x,y,z;
} XYZ;

/*
   Rotate a point p by angle theta around an arbitrary axis r
   Return the rotated point.
   Positive angles are anticlockwise looking down the axis
   towards the origin.
   Assume right hand coordinate system.
*/
XYZ ArbitraryRotate(XYZ p,double theta,XYZ r)
{
   XYZ q = {0.0,0.0,0.0};
   double costheta,sintheta;

   Normalise(&r);
   costheta = cos(theta);
   sintheta = sin(theta);

   q.x += (costheta + (1 - costheta) * r.x * r.x) * p.x;
   q.x += ((1 - costheta) * r.x * r.y - r.z * sintheta) * p.y;
   q.x += ((1 - costheta) * r.x * r.z + r.y * sintheta) * p.z;

   q.y += ((1 - costheta) * r.x * r.y + r.z * sintheta) * p.x;
   q.y += (costheta + (1 - costheta) * r.y * r.y) * p.y;
   q.y += ((1 - costheta) * r.y * r.z - r.x * sintheta) * p.z;

   q.z += ((1 - costheta) * r.x * r.z - r.y * sintheta) * p.x;
   q.z += ((1 - costheta) * r.y * r.z + r.x * sintheta) * p.y;
   q.z += (costheta + (1 - costheta) * r.z * r.z) * p.z;

   return(q);
}
%}

%{
===========================
==== another version using a line segment between 2 points as input
===========================
/*
   Rotate a point p by angle theta around an arbitrary line segment p1-p2
   Return the rotated point.
   Positive angles are anticlockwise looking down the axis
   towards the origin.
   Assume right hand coordinate system.  
*/
XYZ ArbitraryRotate2(XYZ p,double theta,XYZ p1,XYZ p2)
{
   XYZ q = {0.0,0.0,0.0};
   double costheta,sintheta;
   XYZ r;

   r.x = p2.x - p1.x;
   r.y = p2.y - p1.y;
   r.z = p2.z - p1.z;
   p.x -= p1.x;
   p.y -= p1.y;
   p.z -= p1.z;
   Normalise(&r);

   costheta = cos(theta);
   sintheta = sin(theta);

   q.x += (costheta + (1 - costheta) * r.x * r.x) * p.x;
   q.x += ((1 - costheta) * r.x * r.y - r.z * sintheta) * p.y;
   q.x += ((1 - costheta) * r.x * r.z + r.y * sintheta) * p.z;

   q.y += ((1 - costheta) * r.x * r.y + r.z * sintheta) * p.x;
   q.y += (costheta + (1 - costheta) * r.y * r.y) * p.y;
   q.y += ((1 - costheta) * r.y * r.z - r.x * sintheta) * p.z;

   q.z += ((1 - costheta) * r.x * r.z - r.y * sintheta) * p.x;
   q.z += ((1 - costheta) * r.y * r.z + r.x * sintheta) * p.y;
   q.z += (costheta + (1 - costheta) * r.z * r.z) * p.z;

   q.x += p1.x;
   q.y += p1.y;
   q.z += p1.z;
   return(q);
}
%}

%{
=================================
==========Python:
==============================

## PointRotate.py Version 1.02
## Copyright (c) 2006 Bruce Vaughan, BV Detailing & Design, Inc.
## All rights reserved.
## NOT FOR SALE. The software is provided "as is" without any warranty.
#############################################################################
"""
    Return a point rotated about an arbitrary axis in 3D.
    Positive angles are counter-clockwise looking down the axis toward the origin.
    The coordinate system is assumed to be right-hand.
    Arguments: 'axis point 1', 'axis point 2', 'point to be rotated', 'angle of rotation (in radians)' >> 'new point'
    Revision History:
        Version 1.01 (11/11/06) - Revised function code
        Version 1.02 (11/16/06) - Rewrote PointRotate3D function

    Reference 'Rotate A Point About An Arbitrary Axis (3D)' - Paul Bourke        
"""
def PointRotate3D(p1, p2, p0, theta):
    from point import Point
    from math import cos, sin, sqrt

    # Translate so axis is at origin    
    p = p0 - p1
    # Initialize point q
    q = Point(0.0,0.0,0.0)
    N = (p2-p1)
    Nm = sqrt(N.x**2 + N.y**2 + N.z**2)
    
    # Rotation axis unit vector
    n = Point(N.x/Nm, N.y/Nm, N.z/Nm)

    # Matrix common factors     
    c = cos(theta)
    t = (1 - cos(theta))
    s = sin(theta)
    X = n.x
    Y = n.y
    Z = n.z

    # Matrix 'M'
    d11 = t*X**2 + c
    d12 = t*X*Y - s*Z
    d13 = t*X*Z + s*Y
    d21 = t*X*Y + s*Z
    d22 = t*Y**2 + c
    d23 = t*Y*Z - s*X
    d31 = t*X*Z - s*Y
    d32 = t*Y*Z + s*X
    d33 = t*Z**2 + c

    #            |p.x|
    # Matrix 'M'*|p.y|
    #            |p.z|
    q.x = d11*p.x + d12*p.y + d13*p.z
    q.y = d21*p.x + d22*p.y + d23*p.z
    q.z = d31*p.x + d32*p.y + d33*p.z

    # Translate axis and rotated point back to original location    
    return q + p1
    
## END PointRotate3D() ##########################
def test_PointRotate3D():
    from point import Point, PointLocate
    from member import Member, MemberLocate
    from param import dim_print, Prompt, yes_or_no, Warning, ClearSelection
    from macrolib.angle import dtor
    from macrolib.ExceptWarn import formatExceptionInfo
    from macrolib.PrintPtList import formatPtList
    a1 = 45.0
    while 1:
        ClearSelection()
        try:
            Warning("Rotate a Point About a Member Line - The rotation is REVERSED if the view " +\
                    "along the rotation axis is from the member LEFT end toward the member RIGHT end. " + \
                    "A right hand coordinate system is assumed. Positive angles are counter-clockwise " +\
                    "when looking down the axis toward the origin. The member LEFT end is the origin.")
            mem1 = MemberLocate("Select member to rotate about")
            pt1 = PointLocate("Pick point to rotate")
            a1 = Prompt(a1, "Rotation angle")
            pt2 = PointRotate3D(mem1.left.location, mem1.right.location, pt1, dtor(a1))
            Warning("The selected point is %s, %s, %s \nThe new point is %s, %s, %s" % \
                    (dim_print(pt1.x), dim_print(pt1.y), dim_print(pt1.z), dim_print(pt2.x), dim_print(pt2.y), dim_print(pt2.z)))
            print formatPtList("Original Point and Rotated Point Coordinates:", [pt1, pt2])
        except:
            Warning(formatExceptionInfo())
        if not yes_or_no("Continue?"):
                break
## END test_PointRotate3D() #####################
if __name__ == '__main__':
    try:
        test_PointRotate3D()
    finally:
        del test_PointRotate3D

%}
