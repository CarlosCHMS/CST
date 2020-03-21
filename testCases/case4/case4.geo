// Gmsh project created on Sat Mar 21 12:34:16 2020
SetFactory("OpenCASCADE");
//+
Circle(1) = {0.0, 0.0, 0.0, 1, 0, 2*Pi};
//+
Circle(2) = {0.0, 0.0, 0.0, 0.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("wall1") = {1};
//+
Physical Curve("wall2") = {2};
//+
Physical Surface("domain") = {1};
//+
Field[1] = Ball;
//+
Field[1].Radius = 1;
//+
Field[1].VIn = 0.05;
//+
Field[1].VOut = 0.05;
//+
Background Field = 1;
//+
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

