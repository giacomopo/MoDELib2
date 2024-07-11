lnum=1200;

xlo=-0.5; xhi=0.5;
ylo=-0.5; yhi= 0.5;
zlo=-0.5; zhi= 0.5;

lcar1 = 100/lnum;


Point(newp) = {xlo, ylo, zlo,lcar1};    /* Point      1 */
Point(newp) = {xhi, ylo, zlo,lcar1};        /* Point      2 */
Point(newp) = {xhi, yhi, zlo,lcar1};   /* Point      3 */
Point(newp) = {xlo, yhi, zlo,lcar1};  /* Point      4 */
Point(newp) = {xlo, ylo, zhi,lcar1};    /* Point      1 */
Point(newp) = {xhi, ylo, zhi,lcar1};        /* Point      2 */
Point(newp) = {xhi, yhi, zhi,lcar1};   /* Point      3 */
Point(newp) = {xlo, yhi, zhi,lcar1};  /* Point      4 */


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(13) = {-1,-4,-3,-2};
Line Loop(14) = {5,6,7,8};
Line Loop(15) = {2,11,-6,-10};

Line Loop(16) = {3,12,-7,-11};
Line Loop(17) = {4,9,-8,-12};
Line Loop(18) = {1,10,-5,-9};

Plane Surface(19) = {13};
Plane Surface(20) = {14};
Plane Surface(21) = {15};
Plane Surface(22) = {16};
Plane Surface(23) = {17};
Plane Surface(24) = {18};

Surface Loop(25) = {19,20,21,22,23,24};

Volume(1) = {25};


