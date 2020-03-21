# CSU
Repository of CSU (Code of Simulation for Unstructured meshes) project

The objective of the project this project is to create a CFD code using unstructure meshs that be simple enough to be used as laboratory for experiments of new algorithms for CFD. The cases are 2d and the meshes are obtained from converting .msh format to the .mshsim format that is read by CSU.

The .msh format is generated using the gmsh mesh generator. To know more about this powerfull mesher access:
https://gitlab.onelab.info/gmsh/gmsh/

The mesh conversion is made by the command: $ python meshConversion.py PathToYour.msh

The meshConversion.py is based on the gmsh API that can be downloaded in https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/api and must be included in the PATH.

The CSU code has a core that is write in c using the codeblock IDE. This choice keeps that code fast enough to the desired applications. The code also have the python code CSU.py that manage the data and core execution. In the CSU/testCases is possible find several test cases to use. The CSU input file is the .csi. You can run a case using the command:

python3 CSU.py PathToYour.csi

The CSU.py is developed in python 3 using the spyder IDE. It also runs in python 2.
