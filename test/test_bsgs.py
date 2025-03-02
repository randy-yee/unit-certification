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
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    scanRadius =0.5;
    \\ D = 3638703101
    K1= nfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);
    K2= bnfinit(x^4 - 41*x^3 + 587*x^2 - 3427*x + 6773);

    lglat = get_log_lattice_bnf(K2);
    reg1 = get_abs_determinant(lglat);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    totaltime = 0;
    start_time = getabstime();

    bsgs_output= bsgs(K1,cpct_units, B, 18, scanRadius, eps,20,"alltest.txt");
    end_time = getabstime();
    totaltime +=(end_time- start_time);
    print("BSGS time: " ,totaltime, " Expected ", 2600);
    GP_ASSERT_WITHIN_RATIO(totaltime, 2600, 0.15);
    \\\# babystock 18 with scan radius 0.5 runs in about 2600 seconds.
    y = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    glegs = [4.619061865536216760, -16.675036477848972758, 14.169130776523246111;
    3.705556112720428534, 14.697056881773186001, -26.315297544331964588;
    -5.731205932306944782, -9.809986206492638008, -9.633237149198505804];
    print("Testing baby step giant step functions complete");
}

{
    read("test/babystock_data.txt");
    my(K1, K2, O_K, n, r, cpct_units, delta_K, B, disc,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - 85*x^2 + 2750*x - 21391);
    K2= bnfinit(x^3 - 85*x^2 + 2750*x - 21391);
    logdisc = log(abs(K1.disc));
    print("Test 2 - BSGS: Sublattice on complex cubic ", precision(K2.reg,10));
    lglat = get_log_lattice_bnf(K2);
    reg1 = get_abs_determinant(lglat);
    GP_ASSERT_NEAR(reg1, K2.reg, eps);

    prime_scale = 2;
    while(prime_scale < 5,

        scaled_lglat =lglat*prime_scale;
        sumv = scaled_lglat[,1];

        for(j=2, length(lglat), sumv+=lglat[,j]);
        X1 = prec_baby(poldegree(K1.pol), logdisc, infinity_norm(sumv));
        X2 = prec_giant(poldegree(K1.pol), logdisc,abs(reg1),infinity_norm(sumv) );
        req_baby1 = REQ_BABY(K1, get_abs_determinant(scaled_lglat), column_sum(scaled_lglat));
        req_giant1 = REQ_GIANT(K1, get_abs_determinant(scaled_lglat), column_sum(scaled_lglat));
        req_baby2 = REQ_BABY(K1, get_abs_determinant(lglat), column_sum(lglat));
        req_giant2 = REQ_GIANT(K1, get_abs_determinant(lglat), column_sum(lglat));

        print("OLD ", ceil(X1), "  ", ceil(X2));
        print("New scaled: ", req_baby1, "  ", req_giant1);
        print("New non-scaled: ", req_baby2, "  ", req_giant2);
        REQ_BSGS = ceil(max(ceil(req_baby1),ceil(req_giant1)));
        default(realprecision, ceil(REQ_BSGS));
        REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
        eps1 = 2^(-REQ_COMPARE);
        scanRadius = 1;

        cpct_units = cpct_from_loglattice(K1, scaled_lglat, eps);
        bsgs_out= bsgs(K1,cpct_units, B, sqrt(abs(matdet(scaled_lglat))), scanRadius, eps,REQ_BSGS);
        bsgs_out_lattice1 = log_lattice_from_compact_set(K1,bsgs_out);
        GP_ASSERT_NEAR(reg1, get_abs_determinant(bsgs_out_lattice1), eps );
        print("returned log lattice: ", precision(bsgs_out_lattice1,10));
        prime_scale = nextprime(prime_scale+1);

    );
    print("Test 2 Succeeds");
}

{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - x^2 - 3872*x - 91215);
    K2= bnfinit(x^3 - x^2 - 3872*x - 91215);
    print("Test 3 - BSGS: Sublattice on real cubic ", precision(K2.reg,10));
    lglat = get_log_lattice_bnf(K2);

    prime_scale = 2;
    \\# while loop shortened for testing. Change back to 14
    while(prime_scale < 6,
        reg1 = get_abs_determinant(lglat);
        GP_ASSERT_NEAR(reg1, K2.reg, eps);
        scaled_lglat =lglat;
        scaled_lglat[,1] *=prime_scale;
        sumv = scaled_lglat[,1];
        for(j=2, length(lglat), sumv+=lglat[,j]);
        X1 = prec_baby(poldegree(K1.pol), log(abs(K1.disc)), infinity_norm(sumv));
        X2 = prec_giant(poldegree(K1.pol), log(abs(K1.disc)),abs(reg1),infinity_norm(sumv) );
        req_baby1 = REQ_BABY(K1, get_abs_determinant(scaled_lglat), column_sum(scaled_lglat));
        req_giant1 = REQ_GIANT(K1, get_abs_determinant(scaled_lglat), column_sum(scaled_lglat));
        REQ_BSGS = ceil(max(ceil(req_baby1),ceil(req_giant1)));
        default(realprecision, ceil(REQ_BSGS));
        REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
        eps = 2^(-REQ_COMPARE);
        scanRadius = 1;

        cpct_units = cpct_from_loglattice(K1, scaled_lglat, eps);
        bsgs_out= bsgs(K1,cpct_units, B, sqrt(abs(matdet(scaled_lglat))), scanRadius, eps,REQ_BSGS);
        bsgs_out_lattice = log_lattice_from_compact_set(K1,bsgs_out);
        GP_ASSERT_NEAR(reg1, get_abs_determinant(bsgs_out_lattice), eps );

        prime_scale = nextprime(prime_scale+1);
    );
    print("end test 3");
}

{

    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,disc,
        lglat, eps = 10^(-20)
    );
    B = 1;
    K1= nfinit(x^3 - 85*x^2 + 2750*x - 21391);
    K2= bnfinit(x^3 - 85*x^2 + 2750*x - 21391);
    logdisc = log(abs(K1.disc));
    print("Test 3 - BSGS: Sublattice on complex cubic ", precision(K2.reg,10));
    lglat = get_log_lattice_bnf(K2);
    reg1 = get_abs_determinant(lglat);
    GP_ASSERT_NEAR(reg1, K2.reg, eps);

    prime_scale = 2;
    while(prime_scale < 20,

        scaled_lglat =lglat*prime_scale;
        sumv = scaled_lglat[,1];

        for(j=2, length(lglat), sumv+=lglat[,j]);
        X1 = prec_baby(poldegree(K1.pol), logdisc, infinity_norm(sumv));
        X2 = prec_giant(poldegree(K1.pol), logdisc,abs(reg1),infinity_norm(sumv) );
        REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
        default(realprecision, ceil(REQ_BSGS));
        REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
        eps = 2^(-REQ_COMPARE);
        scanRadius = 1;

        cpct_units = cpct_from_loglattice(K1, scaled_lglat, eps);
        bsgs_out= bsgs(K1,cpct_units, B, sqrt(abs(matdet(scaled_lglat))), scanRadius, eps,REQ_BSGS);
        bsgs_out_lattice = log_lattice_from_compact_set(K1,bsgs_out);
        GP_ASSERT_NEAR(reg1, get_abs_determinant(bsgs_out_lattice), eps );
        print("returned log lattice: ", precision(bsgs_out_lattice,10));
        prime_scale = nextprime(prime_scale+1);

    );
    print("Test 2 Succeeds");
}
