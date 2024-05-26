read("src/PmaxLog.py")
read("src/PmaxNormal.py")
read("src/bounds.gp")
read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 15G);
\\ Global variables
eps = 10^(-80);      \\ error tolerance
sqrt2 = sqrt(2);
setrand(121213);

DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;





OUTFILE1 = "data/lpohst-3-0.txt";



\\read("input/test-poly-1-1.gp");  ;
read("input/test-poly-3-0.gp");  signature_string = "3-0"; 
\\read("input/test-poly-4-0.gp");  ;
\\read("input/test-poly-2-1.gp");  ;
\\read("input/test-poly-0-2.gp");  ;

\\read("input/test-poly-1-2.gp");  ;
\\read("input/test-poly-3-1.gp");  ;
\\read("input/test-poly-5-0.gp");  ;
\\read("input/test-poly-0-3.gp");  ;
\\read("input/test-poly-2-2.gp");  ;

\\read("input/test-poly-4-1.gp");  ;
\\read("input/test-poly-0-4.gp");  ;
\\read("input/test-poly-1-3.gp");  ;

SMALLEXAMPLE = 0;
{
    start = 1;
    end   = 10;
    step  = 1;
    loop_ranges = [start, end, step];
    pmax_log_experiment(signature_string, loop_ranges, []);

}
