# 3D-MGTM-HS
## 3D Magnetic Gradient Tensor Modeling at High Susceptibility

by Ouyang Fang, and Longwei Chen

This repository contains the manuscripts and supplementary code for a paper about the open-source software package 3D-MGTM-HS. It also describes an optimal method for 3D magnetic gradient tensor modeling due to its applicability for high magentic susceptibility .

The 3D-MGTM-HS package is written in the FORTRAN programming language and has no external dependencies. 

The manuscript has been published in the journal Minerals on October 2021, and a relevant manuscript can also be found in the journal Geophysics on the January-February 2020 issue.

doi:10.3390/min11101129  (Minerals)

doi:10.1190/GEO2018-0851.1  (Geophysics)



Please cite it as:

Ouyang F., and L. Chen (2021), A FORTRAN program to model magnetic gradient tensor at high susceptibility using contraction integral equation method, MINERALS, 11(10), 1129, doi:10.3390/min11101129

Ouyang F., and L. Chen (2020), Iterative magnetic forward modeling for high susceptibility based on integral equation and Gauss-fast Fourier transform, GEOPHYSICS, 85(1): J1–J13, doi:10.1190/GEO2018-0851.1

You will find A "live" version of the FORTRAN source-code for 3D-MGTM-HS at https://github.com/Yonfou/3D-MGTM-HS. E-mail: FangOuyang92@163.com.

## Running the codes

We provide a FORTRAN program to compute the three magnetic field components and six magnetic gradient tensor components at high susceptibility. The codes consist of 14 subroutines (.f90) and 2 files (.txt): ModelFile.txt and Para.txt. 

In the algorithm, the model is defined as NX×NY×NZ cells. The grid associated with the model is defined as the center of these cells.

Modelfile.txt gives the the magnetic susceptibility at each grid point. In the Modelfile.txt, each line has 1 values (i.e. the magnetic susceptibility) and there are NX×NY×NZ lines in the file. X-coordinates are the inner loop, then the Y-coordinates, and the Z-coordinates are the outer-most loop. Note that the starting points of the X- and Y-coordinates are set to be ZEROS. UNIFORMLY-SPACED SAMPLING is used in ALL THREE directions.

Resultfile.txt saves the calculated magnetic anomalies at Z=Z0 (Z0 is the first point in the z direction). In the Resultfile.txt, each line has 12 values,  including the x-, y-, z-coordinates,  the three magnetic field components, and six magnetic gradient tensor components. There are NX×NY lines in this file.

To run the codes successfully, the Intel Math Kernel Library should be used. 

Here a Sphere model is provided to test the correctness of the codes. 

## License

All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See LICENSE.md for the full license text.

The manuscript text is not open source. The rights to the article content are reserved to the journal Minerals and Geophysics , where it has been accepted for publication.
