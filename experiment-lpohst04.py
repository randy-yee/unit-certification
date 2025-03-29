read("src/PmaxLog.py")
read("src/PmaxNormal.py")
read("src/bounds.gp")
read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 25G);
\\ Global variables
eps = 10^(-80);      \\ error tolerance
sqrt2 = sqrt(2);
setrand(121213);

DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;


INPUT_FILE = "input/experiment-polynomials-0-4";
OUTPUT_FILE = "data/pmax-Bsize-";
\\ if the input file and output file strings are removed, then default files
\\ will be used
{
    sigstring = "0-4";
    OUTPUT_FILE = concat(OUTPUT_FILE, sigstring);
    start = 20;
    end   = 25;
    step  = 1;
    loop_ranges = [start, end, step];

    \\pmax_bound(sigstring, loop_ranges, [INPUT_FILE, OUTPUT_FILE]);
    \\\ # generate data for specified input and output files
    pmax_log_experiment(sigstring, loop_ranges, [INPUT_FILE, OUTPUT_FILE]);

    \\\# default version
    \\pmax_log_experiment(sigstring, loop_ranges, []);
}
