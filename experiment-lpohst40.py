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





OUTFILE1 = "data/lpohst-4-0.txt";

read("input/test-poly-4-0.gp");  ;


\\INPUT_FILE = "input/extra-polynomials-4-0";
\\OUTPUT_FILE = "data/pmax-extra-";
INPUT_FILE = "input/test-poly-4-0.gp";
OUTPUT_FILE = "data/pmax-large-";
\\ if the input file and output file strings are removed, then default files
\\ will be used

SMALLEXAMPLE = 0;
{
    sigstring = "4-0";
    OUTPUT_FILE = concat(OUTPUT_FILE, sigstring);
    start = 19;
    end   = 25;
    step  = 1;
    loop_ranges = [start, end, step];
    pmax_log_experiment(sigstring, loop_ranges, [INPUT_FILE, OUTPUT_FILE]);
    \\pmax_log_experiment(sigstring, loop_ranges,[]);

/*
for(i=1,1,
    \\
    \\ INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    data[i] = specialField40;
    K = nfinit(data[i][1]);
    \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;
    lglat = process_complex_loglattice(K ,data[i][3]);
    cpct_units = cpct_from_loglattice(K, lglat, eps);
    reg1 = unscaled_determinant(K, lglat, []);


    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X = prec_rigorous(poldegree(K.pol), log(abs(K.disc)), log(infinity_norm(sumv)) ,abs(reg1));
    default(realprecision, ceil(X));

    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value", ceil(X));


    lglat_new = lglat; \\modpair1[2] =1;
    unitvector_cpct = cpct_from_loglattice(K, lglat_new, eps);                  \\ computation of compact reps
    tbefore = getabstime();

    \\ -1 indicates to use the usual j-value, the last argument says to limit the size of the lowerbound unit search area
    t_x = getabstime(); indexbound = get_index_bound2(K, lglat_new, eps,-1, 1000000); t_y = getabstime(); boundtime = (t_y-t_x)/60000.0;
    write(OUTFILE1, "Index bound: ", indexbound, ".   bound calc time: ", precision(boundtime,15)  );

    logout = log_pohst_pari(K,lglat_new,unitvector_cpct, indexbound, eps);
    tafter = getabstime();
    outreg = unscaled_determinant(K,logout);
    write(OUTFILE1, "Output Regulator: ", precision(outreg,10 ), "\n  lpohst time ",precision((tafter-tbefore),10), "ms. In mins: ", precision((tafter-tbefore)/60000.0 ,15));


);
*/
}
