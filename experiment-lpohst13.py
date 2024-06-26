read("src/PmaxLog.py")
read("src/PmaxNormal.py")
read("src/bounds.gp")
read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 15G);
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);
setrand(121213);

DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;


\\INPUT_FILE = "input/extra3-polynomials-1-3";
\\OUTPUT_FILE = "data/pmax-extra-";
INPUT_FILE = "input/test-poly-1-3.gp";
OUTPUT_FILE = "data/pmax-large-";
\\ if the input file and output file strings are removed, then default files
\\ will be used
{
    sigstring = "1-3";
    OUTPUT_FILE = concat(OUTPUT_FILE, sigstring);
    start = 23;
    end   = 23;
    step  = 1;
    loop_ranges = [start, end, step];
    pmax_log_experiment(sigstring, loop_ranges, [INPUT_FILE, OUTPUT_FILE]);
    \\pmax_log_experiment(sigstring, loop_ranges, []);
}
