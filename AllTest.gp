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

Excerpt from Pari/GP Manual
By default, arguments are passed by value, not as variables:
modifying a function’s argument in the function body is allowed,
but does not modify its value in the calling scope.
In fact, a copy of the actual parameter is assigned to the formal parameter
when the function is called. (This is not literally true:
a form of copy-on-write is implemented so an object is not duplicated unless
modified in the function.) If an argument is prefixed by a tilde ~ in the
function declaration and the call, it is passed by reference.
(If either the declaration or the call is missing a tilde, we revert to a call by value.)

What does this mean?
- Copy-on-write implies that unless a variable is modified in the function, it will not be copied
- You must prefix the variable with a ~ in both the function declaration and the call to get
  pass by reference behaviour
*/


{
/*
Debugging notes:
currently the original implementation is on msx that uses the non-normalized log vector


*/
read("test/TestingUtility.gp");


\\read("test/test_examples.py");

read("test/test_compact_representation.py");

read("test/test_mlll.py");
read("test/test_vector_methods.py");
read("test/test_log_lattice_functions.py");

read("test/test_lower_bound.py");
read("test/test_minima.py");
read("test/test_neighbours.py");

read("test/test_reduction.py");



read("test/test_bsgs.py");
read("test/test_pmax_log.py");

read("test/test_pmax_normal.py");
read("test/test_hybrid.py");








\\read("test/test_pmax_log_pohst_example.py");

}
