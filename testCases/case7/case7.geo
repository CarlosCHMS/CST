//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 0.4, 0.1, 0};
//+
Rectangle(2) = {0.4, 0, 0, 0.6, 0.1, 0};
//+
Characteristic Length {1, 2, 3, 4, 6, 7} = 0.02;
//+
Physical Curve("T1") = {4};
//+
Physical Curve("T2") = {6};
//+
Physical Curve("wall") = {3, 7, 5, 1};
//+
Physical Curve("intra") = {2};
//+
Physical Surface("domain1") = {1};
//+
Physical Surface("domain2") = {2};
//+
Coherence;
