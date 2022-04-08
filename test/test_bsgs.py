{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
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

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    bsgs_out= bsgs(K1,cpct_units, B, 1/2, eps,REQ_BSGS);
    bsgs_out_lattice = log_lattice_from_compact_set(K1,bsgs_out);
    GP_ASSERT_MAT_NEAR(lglat,bsgs_out_lattice, eps );
}
print("Testing baby step giant step functions complete");
