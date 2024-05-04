{
read("src/FPLowerBound.py");
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("src/VectorMethods.gp");
}

generate_nf_data(poly)=
{
    my(K1 = bnfinit(poly));
    \\print("Test case degree: ", length(K1.zk), "  Sig: ", K1.r1, ",",K1.r2 );
    return([K1.nf, K1.reg]);
}

\\ note that Pohst recommends this value for K, which is based on the length
\\ of the longest vector in the integral basis.
\\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
compute_pohst_K(K1)=
{
    my(varG);
    varG = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    varG = varG~*varG;
    return(2*varG);
}
{
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^3 - x^2 - 29*x + 4;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    GP_ASSERT_NEAR(lbound, 12.318863258, eps);
}

{
    my(K1, reg, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^3 - 44*x - 94;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000); \\autoselect j

    GP_ASSERT_LT(lbound, reg);
    GP_ASSERT_NEAR(lbound, 10.3007897632, eps);\\Magma's bound
}

{
    my(K1, reg, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^4 - x^3 - 6*x^2 -2*x +4;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);

    \\ Note for this example, the bound is not the same as magma, but at least
    \\ it is lower.

    for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));

    \\GP_ASSERT_NEAR(lbound, 3.2009782, eps);
}

{
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    GP_ASSERT_NEAR(lbound, 22.991060641, eps); \\ this is Magma's bound

}

{
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    GP_ASSERT_NEAR(lbound, 53.962667818, eps); \\ this is Magma's bound
}

{
    my(K1, reg, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^4 - 41*x^2 + 359;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, reg);
    GP_ASSERT_NEAR(lbound, 4.52820523359411, eps); \\ this is Magma's bound
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    poly2 = x^4 - x^3 - 11*x^2 - x + 17;
    [K1, reg] = generate_nf_data(poly2);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));

    GP_ASSERT_NEAR(lbound, 10.99192204, eps); \\ this is Magma's bound
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    poly3 = x^4 - x^3 - 8*x^2 + 5*x + 9;
    [K1, reg] = generate_nf_data(poly3);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    GP_ASSERT_NEAR(lbound, 15.7244996475, eps); \\ this is Magma's bound
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    poly4 = x^4 - 2*x^3 - 15*x^2 + x + 5;
    [K1, reg] = generate_nf_data(poly4);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    GP_ASSERT_NEAR(lbound, 10.6957730479, eps); \\ this is Magma's bound
}

{
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^4 - x^3 - 2*x^2 + 24*x - 25;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    GP_ASSERT_NEAR(lbound, 7.004833791, eps); \\ this is Magma's bound
}
{
    print("test lower bound complex cubics:");
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^3 - x^2 + 43*x - 162;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\GP_ASSERT_NEAR(lbound, 4.25012110660, eps); \\ this is Magma's bound

    /*
Field Polynomial: x^3 - x^2 + 43*x - 162. Discriminant: -100003
lower bound: 4.250121106604499314 REG  83.11685151358242906
MAGMA COMPARISIONS:	4.2501211
Field Polynomial: x^3 + 35*x - 92. Discriminant: -100007
lower bound: 4.421135429906412046 REG  177.58209602136350454
MAGMA COMPARISIONS:	4.42113542
Field Polynomial: x^3 - x^2 + 50*x - 356. Discriminant: -100011
lower bound: 4.030525204275756910 REG  77.36579348618080501
MAGMA COMPARISIONS:	4.030525204
Field Polynomial: x^3 - x^2 + 967*x + 6192. Discriminant: -13182675
lower bound: 5.704327046802538089 REG  492.7712274420387190
MAGMA COMPARISIONS:	5.7043270468
*/
}

{
    print("test lower bound complex quartics:");
    \\ in this signature, j = 0 seems to be correct

    my(K1, K2, O_K, n, eps = 10^(-8));
    K1 = nfinit(x^4 - x^3 + 49*x^2 - 6*x + 684);
    K2 = bnfinit(x^4 - x^3 + 49*x^2 - 6*x + 684);
    deg = length(K1.zk);
    print("Test case degree: ", deg, "  Sig: ", K1.r1, ",",K1.r2 );
    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    print(precision(lbound,10), "  ", precision(K2.reg,10));
    print("Warning, this bound does not match Magma");
    \\GP_ASSERT_NEAR(lbound, 4.66728844, eps); \\ this is Magma's bound
}

{
    print("test complex sextic 1:"); \\ in this signature, j = 0 seems to be correct

    my(K1, K2, O_K, n, eps = 10^(-7));
    K1 = nfinit(x^6 - x^5 + 2*x^4 - x^3 - x^2 + 2*x + 1);
    K2 = bnfinit(x^6 - x^5 + 2*x^4 - x^3 - x^2 + 2*x + 1);
    deg = length(K1.zk);
    print("Test case degree: ", deg, "  Sig: ", K1.r1, ",",K1.r2 );
    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(K2.reg,10));

    GP_ASSERT_NEAR(lbound, 1.9003695232, eps); \\ this is Magma's bound

    print("test complex sextic 2:"); \\ in this signature, j = 0 seems to be correct

    my(K1, K2, O_K, n, eps = 10^(-8));
    K1 = nfinit(x^6 - x^5 + 2*x^4 - x^3 + x^2 - 2*x + 1);
    K2 = bnfinit(x^6 - x^5 + 2*x^4 - x^3 + x^2 - 2*x + 1);
    deg = length(K1.zk);
    print("Test case degree: ", deg, "  Sig: ", K1.r1, ",",K1.r2 );
    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(K2.reg,10));

    GP_ASSERT_NEAR(lbound, 0.669980365772, eps); \\ this is Magma's bound
}
{
    print("test complex octics:"); \\ in this signature, j = 0 seems to be correct

    my(K1, kconstant, lbound, deg, eps = 10^(-9));

    poly = x^8 - 4*x^7 + 7*x^6 - 6*x^5 + 3*x^4 - 3*x^3 + 3*x^2 + 1;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));
    GP_ASSERT_NEAR(lbound, 4.593866855758570, eps); \\ this is Magma's bound, match at j=6
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^8 - 5*x^7 + 6*x^6 + 5*x^5 - 7*x^4 - 8*x^3 + 3*x^2 + 5*x + 2;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));
    GP_ASSERT_NEAR(lbound, 5.706608556182, eps); \\ this is Magma's bound, match at j=6
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    poly = x^8 - 6*x^7 + 22*x^6 - 55*x^5 + 103*x^4 - 145*x^3 + 149*x^2 - 97*x + 32;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);

    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 45.565399085719, eps); \\ this is Magma's bound, match at j=6


}

{
    \\ Test against Magma's bound for (5,0), match at j=3
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^5 - 22*x^4 + 175*x^3 - 621*x^2 + 963*x - 534;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 46.1844491760820, eps);

    poly = x^5 - 30*x^4 + 325*x^3 - 1574*x^2 + 3370*x - 2605;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 150.0778326985420, eps);
}
{
    \\ test (3,1) -- Match magma at j = 3
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^5 - 15*x^4 + 80*x^3 - 186*x^2 + 173*x - 31;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 17.569196455888, eps);
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    pol2 = x^5 - 27*x^4 + 246*x^3 - 924*x^2 + 1540*x - 939;
    [K1, reg] = generate_nf_data(pol2);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 113.5553670748, eps);
}
{
    \\ test (1,2) -- Match magma at j = 3
    my(K1, kconstant, lbound, deg, eps = 10^(-9));
    poly = x^5 - 2*x^4 + 4*x^3 - 19*x^2 + 34*x - 17;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    \\for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 7.6667538040973, eps);
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    poly = x^5 - 5*x^4 + 36*x^3 - 109*x^2 + 271*x - 429;
    [K1, reg] = generate_nf_data(poly);
    kconstant = compute_pohst_K(K1);
    deg = length(K1.zk);
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);

    GP_ASSERT_LT(lbound, reg);
    
    for(i = 0, deg-2,print("Set j = ",i, ":  ", precision(lower_regbound(K1, kconstant, eps, i, 2500000000), 10)));
    \\print(precision(lbound,10), "  ", precision(reg,10));

    GP_ASSERT_NEAR(lbound, 38.43652225500929, eps);

}

/*
Code snippet pastable in Magma's calculator:
Use Verbose options to get info on K/number of enumerated elements

R<x> := PolynomialRing(Rationals());
f2 :=x^5 - 5*x^4 + 36*x^3 - 109*x^2 + 271*x - 429;
K := NumberField(f2); Discriminant(K)/6724;
Signature(K);
Lreg := RegulatorLowerBound( K); Lreg;
Regulator(K);
*/

print("testing lower bound complete");
