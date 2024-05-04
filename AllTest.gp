/*
A note about running GP modules:
All paths are provided relative to the root directory unit-certification,
including dependencies within submodules.
If you wish to run anything, it is easiestto do so from this directory,
AllTest.py illustrates how to run test files from this directory.

If you only want to test a subset, either make a new file similar to this one,
or comment out the unwanted tests from the list below.

See https://linux.die.net/man/1/gp for available flags
If you have problems with memory, use the -s flag to specify stack size,
validated for 8GB stack size
*/


{
read("test/TestingUtility.gp");

read("test/test_vector_methods.py");
read("test/test_log_lattice_functions.py");

read("test/test_lower_bound.py");


read("test/test_neighbours.py");



read("test/test_pmax_normal.py");

read("test/test_reduction.py");
read("test/test_pmax_log.py");



read("test/test_compact_representation.py");

read("test/test_bsgs.py");
read("test/test_minima.py");

read("test/test_hybrid.py");

\\read("test/test_mlll.py");
\\read("test/test_pmax_log_pohst_example.py");

}
