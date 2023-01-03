{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

{
    print("Test 1 - Gram-schmidt and reduced lattice test");
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
