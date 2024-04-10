read("src/PmaxNormal.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Set PRECISION:
\p 500;
default(parisizemax, 15G);
setrand(121213);

\\ Global variables
eps = 10^(-80);      \\ error tolerance, used in the pth root search
sqrt2 = sqrt(2);
DEBUG_MAIN = 1;
DEBUG_IDEAL_SEARCH = 0;

OUTFILE1 = "data/tmp-experiment-data-normal-3-0.txt";



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
for(i=1,10,
    K = nfinit(data[i][1]); \\K1 = bnfinit(data[i][1],1);
    \\unit_matrix1 = units_to_matrix(K, K1.fu);
    \\power_units_LLL = unit_matrix;


    lglat = process_complex_loglattice(K ,data[i][3]);
    \\lglat = compute_sublattice(lglat, OUTFILE1, 2)[1];
    cpct_units = cpct_from_loglattice(K, lglat, eps);
    unitmatrix2 = compact_reconstruct(K, cpct_units[1][1],cpct_units[1][2]);
    reg1 = unscaled_determinant(K, lglat);

    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X = prec_rigorous(poldegree(K.pol), log(abs(K.disc)), log(infinity_norm(sumv)) ,abs(reg1));
    default(realprecision, ceil(X));

    for(i=2, length(cpct_units),
        matconcat([unitmatrix2, compact_reconstruct(K, cpct_units[i][1],cpct_units[i][2])] );
    );
    power_units_LLL = Mat(unitmatrix2);
    \\ power_units_LLL =example_gen(K, K1, 31, 1);
    \\power_units_LLL = LLL_reduce_units(K,power_units_LLL);

    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value", ceil(X));


    inreg = unscaled_determinant(K, lglat);
    \\ CONSTRUCT K TO BE USED TO COMPUTE LOWER BOUND AND INDEX BOUND
    Gv = K[5][2]*nfalgtobasis(K, K.zk[length(K.zk)]); Gv = Gv~*Gv;
    \\print("Gv = ", precision(Gv,10));
    pohst_k = 2*Gv;
    \\lbound = lower_regbound(K,  pohst_k,eps);
    \\index_bound = ceil(inreg/(lbound));

    t_x = getabstime(); index_bound = get_index_bound2(K, lglat, eps,-1, 1000000); t_y = getabstime(); boundtime = (t_y-t_x)/60000.0;
    write(OUTFILE1, "Index bound: ", index_bound, ".   bound calc time: ", precision(boundtime,15)  );


    \\
    pohst_output = pohst_check_normal(K,power_units_LLL, index_bound,eps);
    t2 = getabstime();
    write(OUTFILE1, "\n  normal pmax time ",precision(t2-t_x ,10), ". In mins:  ", precision((t2-t_x)/60000.0, 15)  );

    \\print("Pohst_out regulator ", precision(log_determinant(K, pohst_output),10));
);

}
