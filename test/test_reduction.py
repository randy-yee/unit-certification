{
\\ Test cases for reduction algorithm
read("src/CompactRepresentation.py");
read("src/BSGSHelper.py")
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("src/VectorMethods.py");
\\ read("src/Neighbours.py");

}

{
    print("TEST--Reduction 1: ");
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
    print("TEST--Reduction 1: Complete  " , precision(log_from_cpct(G1, cpct_rep),10));
}

{
    print("TEST--SCAN 1: ");
    EXTRA_INFO = 0;
    # setup
    eps = 10^(-9);
    my(G1, G2, O_K, n, eps = 10^(-9));
    G1 = nfinit(x^3 - x^2 - 3872*x - 91215);
    G2 = bnfinit(G1.pol);
    n = poldegree(G1.pol);
    O_K = matid(n);
    reduced_ideal = O_K;

    \\# enumeration on O_K
    ideal_lattice = G1[5][1]*reduced_ideal;
    real_ideal_lattice = embed_real(G1, ideal_lattice);
    LLL_reduced = real_ideal_lattice*qflll(real_ideal_lattice);
    gram_mat = LLL_reduced~ * LLL_reduced;
    scan_elements = qfminim(gram_mat,(100*1.5)^2,,2)[3];

    \\# Using the method below, this element was confirmed a minimum of O_K
    \\# it's l2 norm in Minkowski space is <1149, infinity norm is < 108
    test_minimum = [-35,0,1]~;
    test_min_minkowski = nfeltembed(G1, test_minimum);
    if(EXTRA_INFO,
        print("Norms: l2 : ", ceil(norml2(test_min_minkowski)), "inf : ", ceil(normlp(test_min_minkowski)));
    );

    \\# we establish the pretend target to be kind of close to the log vector of
    \\# our test element
    \\# target ~= [-0.28603985039, 4.674390364, 0.73465946472]
    target = get_normalized_log_vector(G1, test_minimum)+[-0.2, 0, 0];

    minimum_counter = 0;
    norm_deltaK = ceil(((2/Pi)^(G1.r2))*abs(G1.disc)^(1/2)*idealnorm(G1,reduced_ideal));
    for(i=1, length(scan_elements),
        if (is_minimum_alternate(reduced_ideal, scan_elements[,i], G1, eps),
            minimum_counter+=1;
            \\#uncomment these if you wish to view the minima and their distance from target in log space
            \\print(i, "  ", scan_elements[,i], "  ",  ". Distance: ",
            \\    precision(norml2(get_normalized_log_vector(G1, scan_elements[,i])-target),10) ;);
        )
    );
    \\# Sept-28-2024 asserts
    GP_ASSERT_TRUE(length(scan_elements) == 506);
    GP_ASSERT_TRUE(minimum_counter == 10);

    \\# compute a nearby reduced ideal defined by test minimum
    test_ideal = idealdiv(G1, reduced_ideal, test_minimum);                     \\#I_x
    test_min_minkowski = nfeltembed(G1, test_minimum);
    log_vector_tm = get_normalized_log_vector(G1, test_minimum);                \\#m_x

    \\# compute the ideal scaling factor
    reduce_target = exponentiate_logvec(G1.r1+G1.r2-1, log_vector_tm);
    if(EXTRA_INFO,
        print("Test minimum and reduced ideal: ", test_minimum, "  ", test_ideal);
        print("   Log vector: ", precision(log_vector_tm,8));
        print("target vector: ", precision(target,8));
    );

    scan_target = log_vector_tm - target;                        \\#m_x - x
    ideal_scale_vector = exponentiate_logvec(G1.r1+G1.r2-1, unsquare_log_embeddings(G1, scan_target));
    scaled_ideal = mulvec(G1[5][1]*test_ideal, ideal_scale_vector);
    scaled_ideal = embed_real(G1, scaled_ideal);
    LLL_change_of_basis = qflll(scaled_ideal);
    scaled_ideal_LLL = scaled_ideal*LLL_change_of_basis;
    gram_matrix = scaled_ideal_LLL~*scaled_ideal_LLL;

    \\# how far we want to scan, plus compare old and new bounds
    rho = 0.3;
    new_bound = exp(rho)*sqrt(n); \\#(this is the square of the actual bound

    vecholder = scaled_ideal_LLL[,1];
    oldbound = sqrt(n)*exp(2*rho)*sqrt(norml2(vecholder)); \\# not squared

    print("Old bound = ", round(oldbound^2), "  new bound = ", round(new_bound^2));
    scan_elements = qfminim(gram_matrix,new_bound^2,,2)[3];
    scan_elements = test_ideal*LLL_change_of_basis*scan_elements;

    found_minima = List();
    CORRECTNESS_FLAG = 0;
    for(i=1, length(scan_elements),
        eltnorm = abs(nfeltnorm(G1,scan_elements[,i] ));
        if(eltnorm<=norm_deltaK,
            newIdeal = idealdiv(G1, test_ideal, scan_elements[,i]);
            if(checkred_old(newIdeal, G1, eps),
                if(EXTRA_INFO,
                    print(i, "  I_x minimum:", scan_elements[,i], "\nOK-minimum(after mult): ", nfeltmul(G1, scan_elements[,i], test_minimum) ,". Distance: ",
                        precision(norml2(get_normalized_log_vector(G1, scan_elements[,i])+log_vector_tm-target),10) ;);
                    print("Minkowski vector: ", precision(nfeltembed(G1, scan_elements[,i]),10));
                    \\print("Log vector: ", get_normalized_log_vector(G1, scan_elements[,i])+log_vector_tm-target);
                );
                ok_min = nfeltmul(G1, scan_elements[,i], test_minimum);
                if(ok_min == test_minimum, CORRECTNESS_FLAG = 1);
                listput(~found_minima, [scan_elements[,i]]);
            );
        );
    );
    GP_ASSERT_TRUE(CORRECTNESS_FLAG);
    print(found_minima);

    print("SCAN test end");
    /*
    \\# Comments:
    \\# This is a minimum. We check that it gives a non-trivial reduced ideal here
    \\# we then compute the log vector, and compute inputs to reddiv_compact
    \\# interestingly, even if the input vector is exactly the exponentiated
    \\# log vector of a known minimum, reddiv doesn't give pick the right ideal
    \\# we find that
    \\# 2-norm of the element corresponding to [1,0,0] is: 1.4179417013682331946
    \\# 2-norm of the element corresponding to [-35,0,1] is  3.000000000000000000 (this is expected)
    test_minimum = [-35,0,1]~;      \\#
    test_ideal = idealdiv(G1, reduced_ideal, test_minimum);         \\#I_x
    print(test_minimum, "  ", test_ideal);
    test_min_minkowski = nfeltembed(G1, test_minimum);
    log_vector_tm = get_normalized_log_vector(G1, test_minimum);    \\#m_x

    reduce_target = exponentiate_logvec(G1.r1+G1.r2-1, log_vector_tm);
    print(log_vector_tm);
    print("target vector: ", target);
    print("\nNEW RED IDEAL: ", reddiv_compact_invert_u(reduced_ideal,reduce_target,~G1,G1[5][1])[1], "  ", beta);
    */
}

/*
\\ test cases for compact_rep and compact_reconstruct
{
    print("Test--Reduction 2: ");
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
