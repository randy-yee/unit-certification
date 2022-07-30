{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

{   \\\ TEST step incrementation
    my(avec, v, directions, counter = 1);
    avec = [10,10,10];
    directions = [1,1,1];
    v = [0,0,1];
    place_marker = 3;
    while (v != [0,0,0],
        place_marker = increment_with_place_marker(avec, v);
        updateDirections(directions, place_marker);
        GP_ASSERT_EQ(v[1]*avec[2]*avec[3]+v[2]*avec[3]+v[3] == counter);
        counter+=1;
    );
    GP_ASSERT_EQ(counter, 1000);
}

SCREEN("commented test");
/*
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
*/
{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - 67*x^2 + 2032*x - 2053);
    K2= bnfinit(x^3 - 67*x^2 + 2032*x - 2053);
    n = poldegree(K1.pol);
    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);

    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    L = Mat(-6970.84528270648362176158817493947745897457384631054091381075996);
    glegs = Mat(34.8542264135324181088079408746973872948728692315527045690537);
    [L, hashboys] = babystock_scan_jump(y, L, glegs, K1, eps);
    \\GP_ASSERT_EQ(length(Mat(hashboys)~), 209);
    \\babystock_scan_jump

    inc = get_giant_step_increment_vectors_compact(K1, glegs, n, eps);


    incremental_giant_steps(K1, L, glegs, ~hashboys, [10], eps);

}
/*
{

    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    \\ D = 3638703101
    K1= nfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    K2= bnfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    bsgs_output= bsgs(K1,cpct_units, B, 25, eps,20,"alltest.txt");

    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    glegs = [4.619061865536216760, -16.675036477848972758, 14.169130776523246111;
    3.705556112720428534, 14.697056881773186001, -26.315297544331964588;
    -5.731205932306944782, -9.809986206492638008, -9.633237149198505804];

}
*/
print("Testing baby step giant step functions complete");
