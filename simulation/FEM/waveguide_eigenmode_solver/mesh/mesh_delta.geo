// Gmsh project created on Tue Mar 31 09:21:51 2026
SetFactory("OpenCASCADE");
//+
Lx = DefineNumber[ 5, Name "Parameters/Lx" ];
//+
Ly = DefineNumber[ 1, Name "Parameters/Ly" ];
//+
dl = DefineNumber[ 0.1, Name "Parameters/dl" ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Lx/2, 0, 0, 1.0};
//+
Point(3) = {Lx/2, -dl, 0, 1.0};
//+
Point(4) = {Lx, -dl, 0, 1.0};
//+
Point(5) = {Lx, Ly+dl, 0, 1.0};
//+
Point(6) = {Lx/2, Ly+dl, 0, 1.0};
//+
Point(7) = {0, Ly, 0, 1.0};
//+
Point(8) = {Lx/2, Ly, 0, 1.0};
//+
Physical Point("corners", 1) = {7, 1, 2, 3, 4, 5, 8, 6};
//+
Line(1) = {7, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 8};
//+
Line(8) = {8, 7};
//+
Physical Curve("mirrors", 9) = {8, 7, 6, 4, 3, 2};
//+
Physical Curve("inlet", 10) = {1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Physical Surface("field", 11) = {1};
//+
Physical Curve("outlet", 12) = {5};
