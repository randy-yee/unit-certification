{
read("src/PmaxLog.py");

default(realprecision, 1000);
}

{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        output_regulator,pmax_log_output,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - x^2 - 3872*x - 91215);
    K2= bnfinit(x^3 - x^2 - 3872*x - 91215);

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
    REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);

    indexbound = get_index_bound2(K1, lglat, eps,-1, 1000000);
    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    pmax_log_output = log_pohst_pari(K1, lglat,cpct_units, indexbound, eps);
    output_regulator = unscaled_determinant(K1, pmax_log_output);
    GP_ASSERT_NEAR(reg1,output_regulator, eps );
}

print("pmax log tests finished");
