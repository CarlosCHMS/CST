//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Circle(5) = {0.5, 0.5, 0, 0.25, 0, 2*Pi};
//+
Curve Loop(3) = {4, 1, 2, 3};
//+
Curve Loop(4) = {5};
//+
Plane Surface(2) = {3, 4};
//+
Curve Loop(5) = {5};
//+
Plane Surface(3) = {5};
//+
Recursive Delete {
  Surface{1}; 
}
//+
Physical Curve("wall1") = {3, 2, 1, 4};
//+
Physical Surface("domain1") = {2};
//+
Physical Surface("domain2") = {3};
//+
Physical Curve("wall2") = {5};
//+
Characteristic Length {1, 2, 3, 4} = 0.03;
//+
Field[1] = Box;
//+
Field[1].VIn = 0.03;
//+
Field[1].VOut = 0.03;
//+
Field[1].XMax = 1;
//+
Field[1].YMax = 1;
//+
Background Field = 1;
//+
Curve{5} In Surface{2};
