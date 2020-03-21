//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 0.1, 0};
//+
Physical Curve("wall1") = {4};
//+
Physical Curve("wall2") = {2};
//+
Physical Curve("wall3") = {1};
//+
Physical Curve("wall4") = {3};
//+
Physical Surface("domain") = {1};
//+
Characteristic Length {1, 2, 3, 4} = 0.02;
