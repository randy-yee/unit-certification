{
read("src/FPLowerBound.py");
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("src/VectorMethods.gp");
}

{
    my(K1, K2, O_K, n, eps = 10^(-9));
    K1 = nfinit(x^3 - x^2 - 29*x + 4);
    K2 = bnfinit(x^3 - x^2 - 29*x + 4);

    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\GP_ASSERT_NEAR(lbound, 12.318863258, eps);

}

{
    my(K1, K2, O_K, n, eps = 10^(-9));
    K1 = nfinit(x^3 - 44*x - 94);
    K2 = bnfinit(x^3 - 44*x - 94);

    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\GP_ASSERT_NEAR(lbound, 10.3007897632, eps);

}

{
    my(K1, K2, O_K, n, eps = 10^(-9));
    K1 = nfinit(x^4 - x^3 - 6*x^2 -2*x +4);
    K2 = bnfinit(x^4 - x^3 - 6*x^2 -2*x +4);

    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\ Note for this example, the bound is not the same as magma, but at least
    \\ it is lower.
    \\GP_ASSERT_NEAR(lbound, 3.2009782, eps);

}

{
    my(K1, K2, O_K, n, eps = 10^(-9));
    K1 = nfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);
    K2 = bnfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);

    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\ Note for this example, the bound is not the same as magma, but at least
    \\ it is lower.
    \\GP_ASSERT_NEAR(lbound, 22.99106064, eps);

}

{
    my(K1, K2, O_K, n, eps = 10^(-9));
    K1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    K2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    my(Gv, kconstant, lbound);
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]);
    Gv = Gv~*Gv;
    kconstant = 2*Gv;
    \\ note that Pohst recommends this value for K, which is based on the length
    \\ of the longest vector in the integral basis.
    \\ see Pohst 1994, page 104, and also Pohst-Zassenhaus chapter
    lbound = lower_regbound(K1, kconstant, eps, -1, 2500000000);
    GP_ASSERT_LT(lbound, K2.reg);
    \\ Note for this example, the bound is not the same as magma, but at least
    \\ it is lower.
    \\GP_ASSERT_NEAR(lbound, 53.962667818, eps);

}

print("testing lower bound complete");
