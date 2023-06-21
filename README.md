# Infinite-delay-2023
 
This repository contains codes to reproduce the tests in the manuscript:
Scarabel and Vermiglio (2023), "Equations with infinite delay: pseudospectral discretization for numerical stability and bifurcation in an abstract framework", submitted. 
Licensing information is available in the file LICENSE. 

The codes require the MATLAB packages from the references (publicly available online):
* Gautschi, "Gauss–Radau formulae for Jacobi and Laguerre weight functions", Math. Com- put. Simulation, 54 (2000), pp. 403–412, 1999 International Symposium on Computational Sciences, to honor John R. Rice (West Lafayette, IN).
* Weideman, Reddy, "A MATLAB differentiation matrix suite", ACM T.Math.Software, 26 (2000), 465–519.

The bifurcation analysis of nonlinear equations is performed with the package MatCont (v.7.1), publicly available online:
* Dhooge, Govaerts, Kuznetsov, Meijer, and Sautois, "New features of the software MatCont for bifurcation analysis of dynamical systems", Math. Comput. Model. Dyn. Syst., 14 (2008), 147–175. Available at: https://sourceforge.net/projects/matcont/

The codes for the linear tests are obtained in the code "linear_tests.m", which calls the functions:
* "construct_AN_DDE.m": code to construct the matrix AN for DDEs using quadrature
* "construct_AN_DDE_integral.m": code to construct the matrix AN for DDEs using the MATLAB built-in function integral
* "construct_AN_RE.m": code to construct the matrix AN for REs using quadrature
* "PSD_laguerre_standard_nodes.m": code to construct the quadrature nodes and weights and the differentiation matrix for given set of Laguerre nodes.

The codes for the nonlinear examples are collected in:
* "MC_beretta_breda_ext.m", which contains the MatCont continuation instructions to perform the bifurcation analysis of the system defined in "PS_beretta_breda.m"
* "MC_blowflies_inf_ext.m", which contains the MatCont continuation instructions to perform the bifurcation analysis of the system defined in "PS_blowflies_inf_ext.m"

The continuation requires the use of the package MatCont. 
