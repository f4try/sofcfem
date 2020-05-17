// Gmsh project created on Sun Apr 12 02:46:47 2020
SetFactory("OpenCASCADE");

// Rectangle(1) = {0, 0, 0, 2.5e-4, 2e-3, 0};
// Rectangle(2) = {2.5e-4, 0, 0, 1e-4, 2e-3, 0};
// Rectangle(3) = {3.5e-4, 0, 0, 2.5e-4, 2e-3, 0};
// Rectangle(4) = {-1e-4, 0, 0, 1e-4, 5e-4, 0};
// Rectangle(5) = {-1e-4, 1.5e-3, 0, 1e-4 , 5e-4, 0};
// Rectangle(6) = {6e-4, 0, 0, 1e-4, 5e-4, 0};
// Rectangle(7) = {6e-4, 1.5e-3, 0, 1e-4, 5e-4, 0};

ls = 0.5e-4;
Point(1) = {-1e-4, 0, 0, ls};
Point(2) = {2.5e-4, 0, 0, ls};
Point(3) = {3.5e-4, 0, 0, ls};
Point(4) = {7e-4, 0, 0, ls};

Point(5) = {-1e-4, 5e-4, 0, ls};
Point(6) = {0, 5e-4, 0, ls};
Point(7) = {6e-4, 5e-4, 0, ls};
Point(8) = {7e-4, 5e-4, 0, ls};

Point(9) = {-1e-4, 1.5e-3, 0, ls};
Point(10) = {0, 1.5e-3, 0, ls};
Point(11) = {6e-4, 1.5e-3, 0, ls};
Point(12) = {7e-4, 1.5e-3, 0, ls};

Point(13) = {-1e-4, 2e-3, 0, ls};
Point(14) = {2.5e-4, 2e-3, 0, ls};
Point(15) = {3.5e-4, 2e-3, 0, ls};
Point(16) = {7e-4, 2e-3, 0, ls};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Line(4) = {5, 6};
Line(5) = {7, 8};

Line(6) = {9, 10};
Line(7) = {11, 12};

Line(8) = {13, 14};
Line(9) = {14, 15};
Line(10) = {15, 16};

Line(11) = {1, 5};
Line(12) = {4, 8};

Line(13) = {6, 10};
Line(14) = {7, 11};

Line(15) = {9, 13};
Line(16) = {12, 16};

Line(17) = {2, 14};
Line(18) = {3, 15};

Curve Loop(1) = {1, 17, 8, 15,6,13,4,11};
Curve Loop(2) = {2, 17, 9,18};
Curve Loop(3) = {3, 18, 10, 16,7,14,5,12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Physical Curve(5) = {1, 2, 4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};