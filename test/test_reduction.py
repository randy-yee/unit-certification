{
\\ Test cases for reduction algorithm
read("src/CompactRepresentation.py");
read("src/BSGSHelper.py")
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("src/VectorMethods.py");
\\ read("src/Neighbours.py");

}

\\ test cases for compact_rep_buchmann and compact_reconstruct
{
    print("TEST: Reduction 1: ");
    my(G1, G2, O_K, n, eps = 10^(-9));
    G1 = nfinit(x^3 - x^2 - 3872*x - 91215);
    G2 = bnfinit(G1.pol);
    n = poldegree(G1.pol);
    O_K = matid(n);

    
    idealI = [1, 19/27, 25/81; 0, 1/27, 1/81; 0, 0, 1/27];
    vec_u = [131.21434213576893100, 0.17518922387204805097, 31.713109882894650631];
    u_inv = [0.0076211181165340061708, 5.7081136493324741635, 0.031532700630516777569];
    reddiv_output1 = reddiv_compact(idealI, vec_u, G1, G1[5][1]);
    beta1 = reddiv_output1[4];
    logbeta1 = log(abs(nfeltembed(G1, beta1)));

    reddiv_output2 = reddiv_compact_invert_u(idealI, vec_u, G1, G1[5][1]);
    beta2 = reddiv_output1[4];
    logbeta2 = log(abs(nfeltembed(G1, beta2)));
    GP_ASSERT_TRUE(beta1 == beta2);

    print("TEST2: Reductions in compact rep:" );
    unit1 = [35.6255880235477673352175797786846437308, 31.9587525114019920708788256282890935825, -67.584340534949759406096405406973737313]~;
    cpct_rep = compact_rep_full_input(G1, unit1, O_K, eps);
    GP_ASSERT_TRUE(normlp(log_from_cpct(G1, cpct_rep)~ - unit1) < 0.00001);
    print("TEST: Reduction 1: Complete  " , precision(log_from_cpct(G1, cpct_rep),10));
}
/*
\\ test cases for compact_rep_buchmann and compact_reconstruct
{
    print("Test Reduction 1: ");
    my(G1, G2, O_K, n, eps = 10^(-9));
    G1 = nfinit(x^3 - x^2 - 3872*x - 91215);
    G2 = bnfinit(G1.pol);
    n = poldegree(G1.pol);
    O_K = matid(n);

    idealI = matid(n);
    vec_u = [0.3284733319156815252927014, 0.3683539376353983827877492, 8.264842226327170307157500];
    u_inv = [3.044387178002925617315149, 2.714780263839104862703106, 0.1209944452193604657678592];
    reddiv_output1 = reddiv_compact(idealI, vec_u, G1, G1[5][1]);
    beta1 = reddiv_output1[4];
    logbeta1 = log(abs(nfeltembed(G1, beta1)));

    reddiv_output2 = reddiv_compact_invert_u(idealI, vec_u, G1, G1[5][1]);
    beta2 = reddiv_output1[4];
    logbeta2 = log(abs(nfeltembed(G1, beta2)));
    GP_ASSERT_TRUE(beta1 == beta2);
}
*/
