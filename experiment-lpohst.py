read("src/PmaxLog.py")
read("src/PmaxNormal.py")
read("src/bounds.gp")

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 8000000000);
\\ Global variables
eps = 10^(-80);      \\ error tolerance
sqrt2 = sqrt(2);
setrand(121213);

DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;





OUTFILE1 = "data/tmp-experiment-data-lpohst-3-0.txt";



\\read("input/test-poly-1-1.gp");  ;
read("input/test-poly-3-0.gp");  ;
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
for(i=11,12,
    \\
    \\ INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    K = nfinit(data[i][1]);
    \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;
    lglat = process_complex_loglattice(K ,data[i][3]);
    cpct_units = cpct_from_loglattice(K, lglat, eps);
    reg1 = unscaled_determinant(K, lglat);


    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X = prec_rigorous(poldegree(K.pol), log(abs(K.disc)), log(infinity_norm(sumv)) ,abs(reg1));
    default(realprecision, ceil(X));

    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value", ceil(X));
    \\
    \\  This is a good spot to modify the log lattice to test sublattice performance.
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    latticetype = 0;
    \\[lglat_new,modpair1]=compute_sublattice(lglat, OUTFILE1, latticetype);
    \\print(precision(lglat_new,10), "\nMODPAIR", modpair1);

    \\lglat_new = lglat; lglat_new[,1] = 3*lglat_new[,1]; print(precision(lglat_new,10), "\nMODPAIR", modpair1);


    \\write(OUTFILE1,"\n--- Field number ", i, " " , K.pol, "\nModified lglat ", precision(lglat_new,10));
    \\inputreg = unscaled_determinant(K,lglat_new);
    \\write(OUTFILE1," Input Regulator: ", precision(inputreg,10));

    lglat_new = lglat; \\modpair1[2] =1;
    unitvector_cpct = cpct_from_loglattice(K, lglat_new, eps);                  \\ computation of compact reps
    tbefore = getabstime();

    \\ -1 indicates to use the usual j-value, the last argument says to limit the size of the lowerbound unit search area
    t_x = getabstime(); indexbound = get_index_bound2(K, lglat_new, eps,-1, 1000000); t_y = getabstime(); boundtime = (t_y-t_x)/60000.0;
    write(OUTFILE1, "Index bound: ", indexbound, ".   bound calc time: ", precision(boundtime,15)  );

    logout = log_pohst_pari(K,lglat_new,unitvector_cpct, indexbound, eps);
    tafter = getabstime();
    outreg = unscaled_determinant(K,logout);
    \\write(OUTFILE1,"Output Regulator: ", precision(outreg,10 ), "  quot: ", precision(inputreg/outreg,10), "YN? ",norml2(outreg*quot - inputreg) < eps, ". Ratios: ", (modpair1[2]-inputreg/outreg)< eps);
    write(OUTFILE1, "Output Regulator: ", precision(outreg,10 ), "\n  lpohst time ",precision((tafter-tbefore),10), "ms. In mins: ", precision((tafter-tbefore)/60000.0 ,15));


);

}
