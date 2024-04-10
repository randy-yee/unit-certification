
read("experiment-bsgs31.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 15G);
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);


DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;


{
start = 1;
final = 9;
steps = 1;
aux = [0];
run_bsgs_experiment(sigstring, [start,final,steps], b_ranges, aux);

}
