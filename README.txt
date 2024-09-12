The content of the folder
--------------------------
The folder contains the MATLAB source codes of the paper:
"Optimization of high-order harmonic generation by optimal control theory: Ascending a functional landscape in extreme conditions".
Detailed comments have been inserted into the codes.



The content of the subfolder OCT_HHG
-------------------------------------
The subfolder contains the main procedures which solve the OCT problems.

The following script files generate the results presented in the paper (the "main" programs in other computer languages):
harmonic13.m, harmonic14.m, harmonic15.m, harmonic17.m.

The computation was performed by the program OCfx_qn.m.

The folder contains also a slightly improved version, OCfx_qn1.m. This program unnecessarily converges to the same results.
The reason is that the optimization path in the current type of problems is extremely sensitive to roundoff errors,
even as tiny as the order of the machine precision. As a result, the two codes may converge to a different local maximum,
even though the computational process is nearly the same.

The script file example_OCfx_qn1.m demonstrates the application of OCfx_qn1.m for one of the problems.

The data file fieldw_results.mat contains the optimized fields in the frequency domain. The other output variables
of OCfx_qn1.m can be generated from the optimal fields by the script file get_output_params.m, without performing
the full optimization process.


The content of the subfolder quasiNewton
-----------------------------------------
The subfolder contains the optimization procedures of the BFGS method.

The quasi-Newton process is performed by the program quasiNewton.m. One of the input variables is a structure "options",
which contains the details of the optimization process, determined by the user. A detailed explanation of the different fields
in this structure is available in the procedure default_op_qn.m, which generates a default structure for general optimization
problems. This procedure is not used in the present case, but rather another procedure which generates the defaults for OCT
problems, optionsOCqn.m.


The content of the subfolder propagator
----------------------------------------
The folder contains the procedures which implement the propagation by the semi-global propagator (Refs. [55], [56] in the current
paper) in the context of OCT problems. The Arnoldi approach is employed for the computation of the function of matrix.

The program solveOCkr.m contains the skeleton of the propagation procedure.


The content of the subfolder absorbing_potential
-------------------------------------------------
The subfolder contains data files with the absorbing potential and other data required for the propagation.

coulomb_optV240.mat is a data file which contains the required unperturbed potential (Vabs240), the grid (x240),
kinetic energy vector in the p domain (K240), the ground state (fi0240), and the dipole moment which decays to zero slope at
the absorbing boundaries (xabs240). All these are required for the definition of the propagation and optimization problems. It
contains also some additional data.

Vabs.mat is a data file which contains the optimized complex absorbing potential. This potential was added at the boundaries of
the variable V0240 from coulomb_optV240.mat for the generation of the variable Vabs240 from coulomb_optV240.mat.


The content of the subfolder spectral_tools
--------------------------------------------
This subfolder contains the required spectral tools which participate in the computational process, and related programs.
Most of these tools are required only in the context of the present propagation method for ideal efficiency and precision,
and not for other propagation methods.


Don't forget to add the folder and the subfolders to the MATLAB path before running the programs.

Do not hesitate to contact me if required.

Ido Schaefer
