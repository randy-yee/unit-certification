{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

/*
func(~map)={
    mapput(map, 1,20);
    print(mapget(map, 1), "  ", Mat(map));
}
{
    testmap = Map();
    func(~testmap);
    print(Mat(testmap)); breakpoint();
}
*/

{   \\\ TEST step incrementation
    print("Test 1 - incrementer ");
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

{
    print("Test 2 - BSGS on complex cubic ");

    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - x^2 - 3872*x - 91215);
    K2= bnfinit(x^3 - x^2 - 3872*x - 91215);

    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);
    GP_ASSERT_NEAR(reg1, K2.reg, eps);
    scaled_lglat =lglat*2;
    sumv = scaled_lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K1.pol), log(abs(K1.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K1.pol), log(abs(K1.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    default(realprecision, ceil(REQ_BSGS));
    REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);
    scanRadius = 1;

    cpct_units = cpct_from_loglattice(K1, scaled_lglat, eps);
    bsgs_out= bsgs(K1,cpct_units, B, sqrt(abs(matdet(scaled_lglat))), scanRadius, eps,REQ_BSGS);
    bsgs_out_lattice = log_lattice_from_compact_set(K1,bsgs_out);
    GP_ASSERT_NEAR(reg1, unscaled_determinant(K1, bsgs_out_lattice), eps );
    print(precision(lglat,10));
}

{
    print("\nTest 3 - BSGS on complex cubic ");
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    scan_ball_radius = 0.75;
    K1= nfinit(x^3 - 67*x^2 + 2032*x - 2053);
    K2= bnfinit(x^3 - 67*x^2 + 2032*x - 2053);
    n = poldegree(K1.pol);
    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);

    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    L = Mat(-6970.84528270648362176158817493947745897457384631054091381075996);
    glegs = Mat(34.8542264135324181088079408746973872948728692315527045690537);
    hashmap1 = Map();
    hashmap2 = Map();
    my(temp1, temp);
    \\[L, temp] = babystock_scan_jump(y, L, glegs, ~hashmap1, K1, scan_ball_radius, eps);
    [L1, temp1] = incremental_baby_steps(y, L, glegs, ~hashmap2, K1, scan_ball_radius, eps);
    print("Compare babystock set sizes: ", matsize(Mat(hashmap1)), "  ",matsize(Mat(hashmap2))  );

    \\babystock_scan_jump

    inc = get_giant_step_increment_vectors_compact(K1, glegs, n, eps);


    incremental_giant_steps(K1, L1, glegs, hashmap2, [10], eps);
    print("TEST 3 Succeeds")
}


{

    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    scanRadius =1;
    \\ D = 3638703101
    K1= nfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    K2= bnfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    lglat = get_log_lattice_bnf(K2);
    reg1 = unscaled_determinant(K1, lglat);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    bsgs_output= bsgs(K1,cpct_units, B, 25, scanRadius, eps,20,"alltest.txt");

    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    glegs = [4.619061865536216760, -16.675036477848972758, 14.169130776523246111;
    3.705556112720428534, 14.697056881773186001, -26.315297544331964588;
    -5.731205932306944782, -9.809986206492638008, -9.633237149198505804];

}

print("Testing baby step giant step functions complete");
