The source code for this project are in c_src, matlab_src and utilities.

./bin contains .mexa64 binaries, compiled from the code in c_src. These binaries are used in the matlab code.

./utilities contains some utility functions I found from the web. 
arrow.m draws arrows in a plot, which is useful for visualizing camera orientations
rodrigues.m converts a Rodrigues rotation vector into a rotation  matrix

./test contains some functions to test matlab and c code. 
checkIsForest.m checks whether a graph is a forest
checkIsLegitimateObsMask.m is currently unused, but may be used later.
debug.m is a general script that are useful for executing short debugging functions. The contents of debug.m changes constantly.

./c_src contains the source code for the binaries in ./bin and some other C/C++ soure code
calculate_residual.c and calculate_one_residual.c calculates the error residual f of the objective function. f'f is the objective function value.
calculate_jacobian_croute.c and calculate_one_jacobian.c calculates the jacobian matrix of the system at each LM iteration.
C_diagonal_inverse.c computes the block inverse of the C matrix in the Jacobian.
computeSkeletal.cpp computes a skeletal graph for an input bundle data. This is used to preprocess bundle data and drop observations.
build_b_large_scale.c, calculate_jacobian_crout_large_scale.c, read_bundle_output.c, pcg_croute, eval_Axb_large_scale.c, build_b_large_scale.c and schur_eval.c are currently unused, but may be used in the future.

./matlab_src contains the Matlab source code for this project.
build_camera_matrix.m builds either a block-jacobi, MLST or skeletal graph on
the connectivity of the camera graph
calculate_jacobian.m calculates the jacobian at each LM iteration
LM.m is the main function. It performs LM iterations.
my_pcg.m is a slight modification of Matlab's standard PCG route so that some
variables can be output.
plot_camera_evolution.m plots the camera positions and orientations at each
LM iteration.
read_bundle_data.m reads in bundle data from ../data/
read_bundle_output.m reads the output file of LM iterations for analysis
purpose
BundleAdjustment.m, build_skeletal_on_obs_graph.m, bundlerBatchjob.m, calculate_jacobian_backup.m, calculate_residual_backup.m, plot_camera_graph.m, read_normrs.m, rebuild_camera_matrix.m are currently
unused, but may be used in the future.
