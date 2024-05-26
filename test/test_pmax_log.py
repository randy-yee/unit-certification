{
read("src/PmaxLog.py");

default(realprecision, 1000);
}

\\ provide polynomial and optional input lattice
\\ return as list [Number field, input_regulator, bnfregulator ]
\\ if no input lattice provided, input_regulator = bnfregulator
testing_setup_pmax_log(pol, input_lattice = [], eps)={
    my(K1,K2, unit_rank,input_regulator, noinput =false);

    \\ initialize number field and basic invariants
    K1 = nfinit(pol);
    K2 = bnfinit(pol);
    n = poldegree(K1.pol);
    unit_rank = K1.r1+K1.r2 -1;
    \\ if no input lattice is provided, use the bnfinit one
    if(input_lattice == [],
        input_lattice = get_log_lattice_bnf(K2);
        noinput = true;
    );
    GP_ASSERT_TRUE(unit_rank > 0);
    input_regulator = unscaled_determinant(K1, input_lattice);
    ratio = input_regulator/K2.reg;
    GP_ASSERT_NEAR( abs(ratio -round(ratio)) ,0,eps);

    if (noinput,
        return([K1, input_regulator, K2.reg ,input_lattice]);
    ,
        return([K1, input_regulator, K2.reg]);
    );
}


{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,
        output_regulator,pmax_log_output,
        lglat, eps = 10^(-20)
    );


    K1= nfinit(x^3 - x^2 - 3872*x - 91215);
    K2= bnfinit(x^3 - x^2 - 3872*x - 91215);
    n = poldegree(K1.pol);
    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);
    GP_ASSERT_NEAR(reg1, K2.reg, eps);

    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K1.pol), log(abs(K1.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K1.pol), log(abs(K1.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    default(realprecision, ceil(REQ_BSGS));

    \\(n^2+2)log(inf(sumv)) + 2n^2+5
    REQ_COMPARE = ceil((n^2 +2)*log(infinity_norm(sumv))+2*n^2 +5);
    eps = 2^(-REQ_COMPARE);
    print("eps: ", eps);
    indexbound = get_index_bound2(K1, lglat, eps,-1, 1000000);
    cpct_units = cpct_from_loglattice(K1, lglat, eps);

    pmax_log_output = log_pohst_pari(K1, lglat,cpct_units, indexbound, eps);
    output_regulator = unscaled_determinant(K1, pmax_log_output);
    GP_ASSERT_NEAR(reg1,output_regulator, eps );
}

{
    my(K1, K2, O_K, n, r, cpct_units, delta_K, start_regulator, bnf_regulator,
        output_regulator,pmax_log_output, sumv,
        lglat, eps = 10^(-20)
    );
    [K1, start_regulator, bnf_regulator, lglat] =
        testing_setup_pmax_log(x^3 - x^2 - 3872*x - 91215, [], eps);

    lglat[,1] *= 2*3*5*7*11*13*17;
    start_regulator *= 2*3*5*7*11*13*17;

    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K1.pol), log(abs(K1.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K1.pol), log(abs(K1.disc)),abs(start_regulator),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    default(realprecision, ceil(REQ_BSGS));

    \\(n^2+2)log(inf(sumv)) + 2n^2+5
    REQ_COMPARE = ceil((n^2 +2)*log(infinity_norm(sumv))+2*n^2 +5);
    eps = 2^(-REQ_COMPARE);
    print("eps: ", eps);
    indexbound = get_index_bound2(K1, lglat, eps,-1, 1000000);
    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    pmax_log_output = log_pohst_pari(K1, lglat,cpct_units, indexbound, eps);
    output_regulator = unscaled_determinant(K1, pmax_log_output);
    GP_ASSERT_NEAR(bnf_regulator,output_regulator, eps );
}

print("Pmax-log tests finished");
