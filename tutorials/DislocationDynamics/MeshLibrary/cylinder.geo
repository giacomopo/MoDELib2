lnum= 10;
lcar1 = 1/lnum;

R=0.5;  /* radius of cylinder (units of Burgers vector) */ 
H=4*R;   /* height of cylinder (units of Burgers vector) */
x0=0;  y0=0;  z0=0;    /* offset of cylinder axis */

np=40;
For k In {1:np}
    theta=(k-1)/np*2.0*Pi;
    x=R*Cos(theta)+x0;
    y=R*Sin(theta)+y0;
    Point(newp) = {x,y,z0,lcar1}; 
    Point(newp) = {x,y,H+z0,lcar1}; 
EndFor

For k In {1:np}
    nID = k+1;
    If (nID>np)
         nID = 1;
    EndIf
    Line(k) = {2*k-1, 2*nID-1};
    Line(k+np) = {2*k, 2*nID};
    Line(k+2*np) = {2*k-1, 2*k};
EndFor


For k In {1:np}
    nID = k+2*np+1;
    If ( nID>(3*np) )
        nID = 1+2*np;
    EndIf
    Line Loop(k) = {k+2*np, k+np, -nID, -k};
    Plane Surface(k) = {k};
EndFor

Line Loop(np+1) = {1:np};
Plane Surface(np+1) = {np+1};
Line Loop(np+2) = {np+1:2*np};
Plane Surface(np+2) = {np+2};

Surface Loop(1) = {1:np+2};
Volume(1) = {1};