{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

{
    print("--Testing minima functions:");
    print("--Test 1: Gram-schmidt and reduced lattice test");
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    scanRadius =1;
    \\ D = 3638703101
    K1= nfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    K2= bnfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);

    lglat = get_log_lattice_bnf(K2);
    reg1 = get_abs_determinant(lglat);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);

    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    glegs = [4.619061865536216760, -16.675036477848972758, 14.169130776523246111;
    3.705556112720428534, 14.697056881773186001, -26.315297544331964588;
    -5.731205932306944782, -9.809986206492638008, -9.633237149198505804];
    extra_coords = vector(length(glegs), i , extra_log_coordinate(K1.r1, K2.r2, glegs[,i]));
    glegs = matconcat([glegs ; extra_coords]);

    giant_position = 5*glegs[,1] +2*glegs[,2] + 3*glegs[,3];
    divisor = jump_compact(matid(length(K1.zk)), giant_position, K1, length(K1.zk), eps);

    idealLattice = embed_real(K1, K1[5][1]*divisor[1]);
    lll_lat = idealLattice*qflll(idealLattice);
    ortho_basis = gram_schmidt(lll_lat);
    for(i=1, length(ortho_basis),
        for (j =i+1, length(ortho_basis),
            GP_ASSERT_NEAR((ortho_basis[,i]~ * ortho_basis[,j]), 0, eps);
        );
    );
    get_enumeration_bounds(poldegree(K1.pol), lll_lat);

    check_ideal_reduced(K1, divisor[1]);
}

{
    eps = eps = 10^(-20);
    K1 = nfinit(x^3 - x^2 - 10*x + 13);
    ideal = [1, 0,  694/1483; 0, 1, 1201/1483;0, 0,1/1483];

    print("ideal test check on ideal with norm ", idealnorm(K1, ideal));
    GP_ASSERT_EQ(check_ideal_reduced(K1, ideal), checkred_old(ideal, K1, eps));


    K1 = nfinit(x^4 - 43*x^3 + 158*x^2 - 1484*x + 17534);
    ideal = [1, 0,   0,  2166/146303; 0, 1,   0, 77122/146303; 0, 0, 1/2, 72941/292606;0, 0,   0,1/292606];
    print(check_ideal_reduced(K1, ideal));
    print(checkred_old(ideal, K1, eps));

    K2 = nfinit(x^4 - 43*x^3 + 158*x^2 - 1484*x + 17534);
    ideal2 = [1, 0, 1/2, 3874251/7640905; 0, 1, 1/2, 7099593/7640905; 0, 0, 1/2, 7192711/15281810; 0, 0, 0, 1/15281810];
    print("check reduced ideal: ", check_ideal_reduced(K2, ideal2));
    print("check red old:       ", checkred_old(ideal2, K2, eps));
    /*
    ideal check mismatch:
    x^4 - 43*x^3 + 158*x^2 - 1484*x + 17534  [1, 0, 0, 2015465/2521144; 0, 1, 0, 3107/1260572; 0, 0, 1, 177261/315143; 0, 0, 0, 1/5042288]  1/3138550867693340381917894711603833208051177722232017256448
    ideal check mismatch:
    x^4 - 43*x^3 + 158*x^2 - 1484*x + 17534  [1, 0, 1/2, 3874251/7640905; 0, 1, 1/2, 7099593/7640905; 0, 0, 1/2, 7192711/15281810; 0, 0, 0, 1/15281810]  1/3138550867693340381917894711603833208051177722232017256448

    ideal check mismatch:
    x^4 - 43*x^3 + 158*x^2 - 1484*x + 17534  [1, 2/5, 2/5, 143573/157975; 0, 1/5, 0, 21174/157975; 0, 0, 1/10, 3277/157975; 0, 0, 0, 1/315950]  1/3138550867693340381917894711603833208051177722232017256448

    x^5 - 2*x^4 + 4*x^3 - 19*x^2 + 34*x - 17  [1, 0, 0, 1/7, 54/77; 0, 1, 0, 0, 23/77; 0, 0, 1, 4/7, 36/77; 0, 0, 0, 1/7, 3/77; 0, 0, 0, 0, 1/77]  1/43556142965880123323311949751266331066368
    */
}
