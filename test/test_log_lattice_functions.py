{
 read("src/LogarithmLattice.py");
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("test/TestingUtility.gp");
\\ read("src/VectorMethods.gp");

}



{ \\ test cases for testequalOK
    my(G1, G2, ideal1, ideal2, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    n = poldegree(G1.pol);
    ideal1 = matid(n);
    GP_ASSERT_TRUE(testequalOK(ideal1, G1));
    ideal2 = ideal1;
    ideal2[,2] += 2*nfalgtobasis(G1, G1.zk[3]);
    ideal2[,3] = ideal2[,3] + 3*nfalgtobasis(G1, G1.zk[1]) + 5*nfalgtobasis(G1, G1.zk[4]);
    GP_ASSERT_TRUE(testequalOK(ideal2, G1));

    \\G2 = nfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);

}
{ \\ test cases for minkowski_absval
    my(v1, v2, v3,G1, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    v1 = nfeltembed(G1, G2.fu[1]);
    v2 = nfeltembed(G1, G2.fu[2]);
    v3 = nfeltembed(G1, G2.fu[3]);

    GP_ASSERT_VEC_NEAR(log(minkowski_absval(v1,G1.r1)), [13.10366530105, 2.932931306,-14.936957398, -1.099639208 ] ,eps);
    GP_ASSERT_VEC_NEAR(log(minkowski_absval(v2,G1.r1)), [3.83927141135,-20.380049263,  1.184114143, 15.356663708] ,eps);
    GP_ASSERT_VEC_NEAR(log(minkowski_absval(v3,G1.r1)), [20.68882227729, 0.206698451, 20.858979979,-41.754500708] ,eps);

}

{ \\ test cases for minkowski_absval
    my(v1, v2, v3,G1, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    v1 = nfalgtobasis(G1, G2.fu[1]);
    v2 = nfalgtobasis(G1, G2.fu[2]);
    v3 = nfalgtobasis(G1, G2.fu[3]);

    GP_ASSERT_VEC_NEAR(valuationvec(G1,v1), [490737.82411335640698, 18.78260745615617,3.258080610544176 E-7, 0.33299120237759 ] ,eps);
    GP_ASSERT_VEC_NEAR(valuationvec(G1,v2), [46.49158885340465,1.409473983671455 E-9,  3.267790743845232, 4669972.3504438980915] ,eps);
    GP_ASSERT_VEC_NEAR(valuationvec(G1,v3), [966142867.58891563714, 1.229611728039522, 1145354438.81852736414298,7.3493807553061848 E-19] ,eps);

    GP_ASSERT_VEC_NEAR(valuationvec(G1,G2.fu[1]), [490737.82411335640698, 18.78260745615617,3.258080610544176 E-7, 0.33299120237759 ] ,eps);
    GP_ASSERT_VEC_NEAR(valuationvec(G1,G2.fu[2]), [46.49158885340465,1.409473983671455 E-9,  3.267790743845232, 4669972.3504438980915] ,eps);
    GP_ASSERT_VEC_NEAR(valuationvec(G1,G2.fu[3]), [966142867.58891563714, 1.229611728039522, 1145354438.81852736414298,7.3493807553061848 E-19] ,eps);

}

{ \\ test cases for absoluteval_nvec
    my(v1, v2, v3,G1, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    v1 = nfalgtobasis(G1, G2.fu[1]);
    v2 = nfalgtobasis(G1, G2.fu[2]);
    v3 = nfalgtobasis(G1, G2.fu[3]);

    GP_ASSERT_VEC_NEAR(absoluteval_nvec(G1,v1), [490737.82411335640698, 18.78260745615617,3.258080610544176 E-7, sqrt(0.33299120237759),sqrt(0.33299120237759)  ] ,eps);
    GP_ASSERT_VEC_NEAR(absoluteval_nvec(G1,v2), [46.49158885340465,1.409473983671455 E-9,  3.267790743845232, sqrt(4669972.3504438980915),sqrt(4669972.3504438980915)] ,eps);
    GP_ASSERT_VEC_NEAR(absoluteval_nvec(G1,v3), [966142867.58891563714, 1.229611728039522, 1145354438.81852736414298,sqrt(7.3493807553061848 E-19),sqrt(7.3493807553061848 E-19)] ,eps);
}

{ \\ test cases for unsquare_log_embeddings
    my(v1, v2, v3, G1, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    v1 = nfalgtobasis(G1, G2.fu[1]);
    v2 = nfalgtobasis(G1, G2.fu[2]);
    v3 = nfalgtobasis(G1, G2.fu[3]);

    GP_ASSERT_VEC_NEAR(log(abs(G1[5][1]*v1)), real(unsquare_log_embeddings(G1, G2[3][,1])) ,eps);
    GP_ASSERT_VEC_NEAR(log(abs(G1[5][1]*v2)), real(unsquare_log_embeddings(G1, G2[3][,2])) ,eps);
    GP_ASSERT_VEC_NEAR(log(abs(G1[5][1]*v3)), real(unsquare_log_embeddings(G1, G2[3][,3])) ,eps);

}

{ \\ test cases for units_to_matrix and unsquare_log_embeddings
    my(G1, G2, Umat, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    \\ take the units and convert to coefficient matrix. Then use log_determinant
    \\ to see if we can obtain the regulator
    Umat = units_to_matrix(G1, G2.fu);
    GP_ASSERT_NEAR(log_determinant(G1, Umat), G2.reg,eps);

    G1 = nfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);
    G2 = bnfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);
    Umat = units_to_matrix(G1, G2.fu);
    GP_ASSERT_NEAR(log_determinant(G1, Umat), G2.reg,eps);

}

{ \\ test cases for get_log_lattice_bnf, process_complex_loglattice and unscaled determinant
    my(G1, G2, Umat, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    logarithm_lattice = get_log_lattice_bnf(G2);
    GP_ASSERT_NEAR(unscaled_determinant(G1, logarithm_lattice), G2.reg,eps);

    G1 = nfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);
    G2 = bnfinit(x^6 - 9*x^5 + 40*x^4 - 95*x^3 + 132*x^2 - 101*x + 31);
    logarithm_lattice = get_log_lattice_bnf(G2);
    GP_ASSERT_NEAR(unscaled_determinant(G1, logarithm_lattice), G2.reg,eps);

    GP_ASSERT_NEAR(unscaled_determinant(G1,process_complex_loglattice(G1, G2[3])), G2.reg, eps);
}


{ \\ test cases for unsquare_log_embeddings
    my(v1, v2, v3, G1, eps = 10^(-9));
    G1 = nfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);
    G2 = bnfinit(x^5 - 15*x^4 + 56*x^3 - 65*x^2 + 48*x - 15);

    v1 = nfeltembed(G1, G2.fu[1]);
    v2 = nfeltembed(G1, G2.fu[2]);
    v3 = nfeltembed(G1, G2.fu[3]);

    GP_ASSERT_VEC_NEAR(get_real_vec(G1, v1), [490737.82411335640698, 18.7826074561561, 3.2580806105441 E-7, -0.429016631231, 0.6942097196683] ,eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, v2), [-46.4915888534046532, 1.409473983671 E-9, -3.26779074384, -1696.519312386000, 2542.0005750568],eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, v3), [-966142867.588915637142527, -1.22961172803952, 1145354438.818527364142, 7.34729004979520 E-10, 9.64390709361 E-10],eps);

    M = G1[5][1];
    Mreal = embed_real(G1,M);

    GP_ASSERT_VEC_NEAR(get_real_vec(G1, nfeltembed(G1, G1.zk[1])), Mreal[,1] ,eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, nfeltembed(G1, G1.zk[2])), Mreal[,2] ,eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, nfeltembed(G1, G1.zk[3])), Mreal[,3] ,eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, nfeltembed(G1, G1.zk[4])), Mreal[,4] ,eps);
    GP_ASSERT_VEC_NEAR(get_real_vec(G1, nfeltembed(G1, G1.zk[5])), Mreal[,5] ,eps);

}

{ \\ test cases for truncatereal
    my(eps1, realnum);
    realnum = 1/3.0;
    eps = 10^(-20);
    GP_ASSERT_NEAR(truncatereal(realnum, 20), 0.333333333333333333330,eps);

}


{ \\ test cases for getidealdenom
    my(mat);
    mat = matid(8);

    GP_ASSERT_EQ(get_ideal_denom(1/21*mat), 21);
    mat = 2*mat;
    GP_ASSERT_EQ(get_ideal_denom(1/20*mat), 10);

}

/*
Untested functions:

checkred
checkred_old
ideal_contains1
is_vec_in_lattice
*/

print("Testing logarithm lattice functions complete")
