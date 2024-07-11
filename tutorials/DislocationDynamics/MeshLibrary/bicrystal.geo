lnum=100;

xlo=-1;  xmid= 0; xhi=1;
ylo=-0.5; yhi= 0.5;
zlo=-0.5; zhi= 0.5;

lcar1 = 100/lnum;


Point(newp) = {xmid,yhi,zhi,lcar1};    /* Point      1 */
Point(newp) = {xmid,yhi,zlo,lcar1};        /* Point      2 */
Point(newp) = {xlo,yhi,zhi,lcar1};   /* Point      3 */
Point(newp) = {xlo,ylo,zhi,lcar1};  /* Point      4 */
Point(newp) = {xmid,ylo,zhi,lcar1};   /* Point      5 */
Point(newp) = {xmid,ylo,zlo,lcar1};       /* Point      6 */
Point(newp) = {xlo,yhi,zlo,lcar1};       /* Point      7 */
Point(newp) = {xlo,ylo,zlo,lcar1};      /* Point      8 */

Point(newp) = {xhi,yhi,zhi,lcar1};    /* Point      9 */
Point(newp) = {xhi,yhi,zlo,lcar1};        /* Point      10 */
Point(newp) = {xhi,ylo,zhi,lcar1};   /* Point      11 */
Point(newp) = {xhi,ylo,zlo,lcar1};       /* Point      12 */


Line(1) = {3,1};
Line(2) = {3,7};
Line(3) = {7,2};
Line(4) = {2,1};
Line(5) = {1,5};
Line(6) = {5,4};
Line(7) = {4,8};
Line(8) = {8,6};
Line(9) = {6,5};
Line(10) = {6,2};
Line(11) = {3,4};
Line(12) = {8,7};

Line(13) = {1,9};
Line(14) = {2,10};
Line(15) = {5,11};
Line(16) = {6,12};
Line(17) = {9,11};
Line(18) = {10,12};
Line(19) = {9,10};
Line(20) = {11,12};


Line Loop(21) = {-6,-5,-1,11};
Plane Surface(22) = {21};
Line Loop(23) = {4,5,-9,10};
Plane Surface(24) = {23};
Line Loop(25) = {-3,-12,8,10};
Plane Surface(26) = {25};
Line Loop(27) = {7,12,-2,11};
Plane Surface(28) = {27};
Line Loop(29) = {-4,-3,-2,1};
Plane Surface(30) = {29};
Line Loop(31) = {8,9,6,7};
Plane Surface(32) = {31};

Line Loop(33) = {13,19,-14,4};
Plane Surface(34) = {33};
Line Loop(35) = {5,15,-17,-13};
Plane Surface(36) = {35};
Line Loop(37) = {-9,16,-20,-15};
Plane Surface(38) = {37};
Line Loop(39) = {14,18,-16,10};
Plane Surface(40) = {39};
Line Loop(41) = {17,20,-18,-19};
Plane Surface(42) = {41};



Surface Loop(43) = {22,32,-26,30,24,-28};
Surface Loop(44) = {34,36,38,40,-24,42};

Volume(1) = {43};
Volume(2) = {44};


