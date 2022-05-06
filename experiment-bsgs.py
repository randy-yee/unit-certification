
read("src/BabyStepGiantStep.py")
read("src/CompactRepresentation.py");
read("src/bounds.gp")
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 8000000000);
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);


DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;


OUTFILE1 = "data/tmp-experiment-data-bsgs-3-0.txt";



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

{
for(i=9, 10,

    \\
    \\ INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    K = nfinit(data[i][1]);
    \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;
    lglat = process_complex_loglattice(K ,data[i][3]);

    reg1 = unscaled_determinant(K, lglat);


    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K.pol), log(abs(K.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K.pol), log(abs(K.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));

    default(realprecision, ceil(REQ_BSGS));
    REQ_COMPARE = ceil((poldegree(K.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);

    \\
    \\  This takes the log lattice and modifies it so that we get a sublattice.
    \\  The modification depends on the 'latticetype', see function compute_sublattice
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\latticetype = 0;
    \\[lglat_new,modpair1]=compute_sublattice(lglat, OUTFILE1, latticetype);
    \\print(precision(lglat_new,10), "\nMODPAIR", modpair1);


    \\ This line provides an example of modifying the log lattice to have a
    \\ particular index divisor
    \\lglat_new = lglat; lglat_new[,1] = 15*lglat_new[,1]; print(precision(lglat_new,10), "\nMODPAIR");
    lglat_new = lglat;

    \\inputreg = unscaled_determinant(K,lglat_new);
    print("input determinant ", precision(unscaled_determinant(K,lglat_new),10));
    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value", ceil(REQ_BSGS));

    cpct_units = cpct_from_loglattice(K, lglat_new, eps);

    B = 1;          \\ 1 means you scan the whole region

    t9 = getabstime();
    bsgs_output= bsgs(K,cpct_units, B, 1/2, eps,REQ_BSGS,OUTFILE1);
    t10 = getabstime();
    bsgs_out_lattice = log_lattice_from_compact_set(K,bsgs_output);
    print(precision(unscaled_determinant(K, bsgs_out_lattice),10));
    print(precision(reg1,10));
    write(OUTFILE1, "\n  Total BSGS time ",precision(t10-t9,10), "  In mins: " ,precision((t10-t9)/60000.0,10) );



);

}
