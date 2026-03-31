//+
Lx = DefineNumber[ 2, Name "Parameters/Lx" ];
//+
Ly = DefineNumber[ 5, Name "Parameters/Ly" ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Lx, 0, 0, 1.0};
//+
Point(3) = {Lx, Ly, 0, 1.0};
//+
Point(4) = {0, Ly, 0, 1.0};
//+
Physical Point("corners", 1) = {4, 3, 2, 1};
//+
Line(1) = {4, 1};
//+
Line(2) = {2, 3};
//+
Physical Curve("mirrors", 3) = {1, 2};
//+
Physical Curve(" mirrors", 3) -= {1, 2};
//+
Line(3) = {4, 3};
//+
Line(4) = {2, 1};
//+
Physical Curve("mirrors", 5) = {3};
//+
Physical Curve("mirrors", 5) += {4};
//+
Physical Curve("inlet", 6) = {1};
//+
Physical Curve("outlet", 7) = {2};
//+
Curve Loop(1) = {3, -2, 4, -1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("surface", 8) = {1};
