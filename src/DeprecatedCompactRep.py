/*
-- find_coprime_divisor_qfminim
-- find_coprime_divisor_lazy
-- find_coprime_divisor_neighbors



*/
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ The section of functions below is related to UNUSED functionality
\\ for ensuring that the denominators in a compact representation
\\ are coprime toa particular integer value
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUTS:
\\ - G a number field
\\ - y an ideal and
\\ - u a real vector. the pair (y,u) is a divisor
\\ - sv is the numerical vector from an LLL basis
\\ - beta is the corresponding coefficient vector
\\ - eps is an error
\\ OUTPUTS:
\\ - a reduced ideal new_y
\\ - beta the element b satisfying (1/beta)y = new_y
\\ - lmu a log vector representing a distance from the original.
minima_from_LLL(G, y, u, sv, beta, eps)={
    my(new_y, nu, lmu, insidepoints,ip_ctr);

    \\UPDATE the ideal, and distance tracking elements nu, lmu,
    \\ SCAN NORMED BODY OF 1 FOR ELEMENTS
    new_y = idealdiv(G,y,beta);
    nu=abs(sv)~;                                                                \\ nu is a t_VEC
    lmu = log(nu)-log(u);                                                       \\ dimension r t_VEC
    \\insidepoints = Mat(nfalgtobasis(G,1));
    insidepoints = [nfalgtobasis(G,1)];
    while(length(insidepoints) != 0,                                            \\ loop until ideal is reduced (until N(1) is trivial)
        ip_ctr = 1;
        while(type(nfbasistoalg(G,insidepoints[ip_ctr]))==t_FRAC , ip_ctr+=1;);   \\ make sure not to choose a rational number.

        beta  = nfeltmul(G, beta, insidepoints[ip_ctr]);
        new_y = idealdiv(G, new_y, insidepoints[ip_ctr]);
        new_y = idealhnf(G,new_y);

        lmu += log(abs(nfeltembed(G, insidepoints[ip_ctr])));
        nu = pointwise_vector_mul(nu, abs(nfeltembed(G,insidepoints[ip_ctr])));
        insidepoints = cubescan_unequal(new_y, valuationvec(G,1), G, eps);

    );
    \\print("NU CHECK: ", precision(abs(pointwise_vector_mul(nfeltembed(G, beta),u) ),10) );
    \\print("NU CHECK: ", precision(nu,10));
    return([new_y, beta, nu, lmu]);

};


\\\ Old sub function used in finding compact representations thta avoid a certain prime p in the denominator
p_avoid_reddiv_compact(~y,~u,~G,~M1, p_avoid=1)={
    my(y1, ideal_uY, numerical_mat_Y, red1, shortest_vec, nu, lmu,
        ideal_denom,vec_ctr,beta_found =0,
        comp = 2^(-ceil((poldegree(G.pol)^2 +2)*log(infinity_norm(u))+2*poldegree(G.pol)^2 +5))
    );

    numerical_mat_Y = M1*y;                                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(~numerical_mat_Y, ~u);                                    \\ ideal u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = y*LLL_change_of_basis;                                        \\ LLL basis, coords wrt the integral basis

    beta= LLLcoeffmat[,1];                                                      \\ beta holds coordinates of mu wrt the integral basis
    \\ need to scan to make sure the first basis vector is a shortest one
    real_mat_uY = embed_real(~G, ~LLL_numerical_uY);
    enum_result = qfminim(real_mat_uY~*real_mat_uY, norml2(real_mat_uY[,1])-comp,,2 );

    true_shortest = qfminim(real_mat_uY~*real_mat_uY,,,2 );
    /* NOTE THIS CHECK IS SLOW IN FIELDS WITH LARGE DEGREE (see pohst example)
    \\
    */
    \\boolA = (enum_result[1]!=2 && !is_minimum(LLLcoeffmat,beta , G, comp));

    boolB = (length(enum_result[3])!=0 && !is_minimum(LLLcoeffmat,beta , G, comp));

    if(boolB,
        short_index =1;
        short_length = norml2(real_mat_uY*enum_result[3][,1]);
        for(j=1, length(enum_result[3]),
            iter_length = norml2(real_mat_uY*enum_result[3][,j]);
            if(iter_length < short_length,
                short_index = j; short_length = iter_length;
            );
        );
        beta = LLLcoeffmat*enum_result[3][,short_index];
        if(!is_minimum(LLLcoeffmat, beta, G, comp),
            print("elements found in normed body of supposed minimum!!!!");
            breakpoint()
        );
        shortest_vec = LLL_numerical_uY*enum_result[3][,short_index];
    , \\else
        shortest_vec = LLL_numerical_uY[,1];                                        \\ vector of complex embeddings for the shortest vector of u*y, denoted mu

    );

    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);                      \\ the reduced ideal y / mu, in hnf form

    \\\#Ran into trouble using abs(shortest_vec) with massive precision loss
    \\\# instead use alternate formula u*psi(beta)
    \\nu=abs(shortest_vec)~;                                                      \\ nu is a t_VEC of dimension r, (complex coordinates are not squared)
    nu = pointwise_vector_mul(abs(M1*beta),u)~;
    \\GP_ASSERT_VEC_NEAR(nu,abs(shortest_vec), comp  );
    lmu = log(nu)-log(u);                                                       \\ expect equal to log(nfeltembed(G, beta))
    \\GP_ASSERT_VEC_NEAR(lmu,log(nfeltembed(G, beta) ),10);

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ If p_avoid is not equal to 1, then we need to find an ideal with coprime denominator
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(p_avoid !=1,

        \\ USES LAZY METHOD, CHECK IF ANY OTHER LLL BASIS ELTS ARE OKAY TO USE
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        ideal_denom = get_ideal_denom(new_y);
        if( gcd(p_avoid, ideal_denom)!= 1,
            [new_y, beta, nu, lmu, ideal_denom]=find_coprime_divisor_lazy(G, y, u, new_y, p_avoid, LLLcoeffmat, LLL_numerical_uY , eps);
            breakpoint();
        );

        \\ USES QFMINIM TO TRY TO FIND AN IDEAL WITH COPRIME DENOMINATOR
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        if(gcd(p_avoid,ideal_denom)!=1,
            my(rmat = embed_real(G,LLL_numerical_uY));                          \\ nxn real matrix of the LLL reduced ideal u*y

            [beta_found, new_y, nu, beta] = find_coprime_divisor_qfminim(G, y, LLL_numerical_uY, LLLcoeffmat, rmat, p_avoid, EXPANSION_LIMIT, eps );
            if(beta_found ==1,
                lmu = log(nu)-log(u);
                ideal_denom = get_ideal_denom(new_y);
            , \\else
                ideal_denom = p_avoid;
            );
        );

        \\ FINALLY, USE NEIGHBOURS, WHICH IS EXHAUSTIVE IF PREVIOUS METHODS HAVE FAILED
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        if(gcd(p_avoid, ideal_denom)!= 1,
            print("Using neighbours to find coprime denominator");
            [new_y, nu, lmu, beta, ideal_denom] = find_coprime_divisor_neighbors(G,y, u, LLLcoeffmat, ideal_denom, p_avoid,eps);
        );
    );
    return([new_y, nu, lmu, beta]);                                                     \\
}


\\ subalgorithm for finding a coprime divisor by just checking the other elements of the LLL basis
find_coprime_divisor_lazy(G, y, u, new_y, p_avoid, LLLcoeffmat, LLL_numerical_uY ,eps)={
    my(vec_ctr =1, shortest_vec, beta, nu, lmu, ideal_denom);
    ideal_denom = get_ideal_denom(new_y);
    while(gcd(p_avoid, ideal_denom)!= 1 && vec_ctr < length(LLLcoeffmat),
        if(vec_ctr ==1, print("-- NON-COPRIME DENOMS, CHECK LLL BASIS"));
        vec_ctr +=1;

        shortest_vec = LLL_numerical_uY[,vec_ctr];                                      \\ change mu to reflect choice of beta
        beta = LLLcoeffmat[,vec_ctr];                                               \\ select beta as corresponding LLL vector

        [new_y, beta, nu, lmu]=minima_from_LLL(G, y, u, shortest_vec, beta, eps);   \\ update all necessary variables to track the ideal
        ideal_denom = get_ideal_denom(new_y);                                       \\ get new ideal denominator.
    );
    return([new_y, beta, nu, lmu, ideal_denom]);
}



find_coprime_divisor_neighbors(G, y, u, LLLcoeffmat, ideal_denom, p_avoid,eps)={
            my(
                beta =LLLcoeffmat[,1],
                nhood_ctr = 0,
                beta_neighbours = Set([beta]),
                nu, lmu, new_ideal,
                unsearched_neighbours=[],
                new_neighbours;
            );
            print("-- USING NEIGHBORS METHOD ");

            new_neighbours = NEIGHBORS(G, y, beta, eps, p_avoid);
            unsearched_neighbours = concat(unsearched_neighbours, new_neighbours);

            while(gcd(p_avoid, ideal_denom) != 1,
                nhood_ctr+=1;

                if(nhood_ctr > length(new_neighbours),
                    print(" - NEXT NEIGHBOUR LAYER NEEDED");

                    new_neighbours = Set();
                    while( length(unsearched_neighbours) != 0,
                        beta = unsearched_neighbours[1];
                        unsearched_neighbours = unsearched_neighbours[^1];

                        if(setsearch(beta_neighbours,beta) == 0,
                            new_neighbours = setunion(new_neighbours, NEIGHBORS(G, y, beta, eps, p_avoid));
                            beta_neighbours = setunion(beta_neighbours, [beta]);
                        );
                    );

                    unsearched_neighbours = concat(unsearched_neighbours,new_neighbours);
                    nhood_ctr=1;
                );

                if(nhood_ctr > length(new_neighbours),
                    ideal_denom = p_avoid;
                ,
                    new_ideal = idealdiv(G, y, new_neighbours[nhood_ctr]); new_ideal = idealhnf(G,new_ideal);
                    ideal_denom = get_ideal_denom(new_ideal);
                );
            );
            beta = new_neighbours[nhood_ctr];
            nu = abs(pointwise_vector_mul(nfeltembed(G,new_neighbours[nhood_ctr]), u));
            lmu = log(abs(nfeltembed(G, new_neighbours[nhood_ctr])));

            return([new_ideal, nu, lmu, beta, ideal_denom]);

}

find_coprime_divisor_qfminim(G, coeff_y, numerical_LLL_uy, coeff_LLL_uy, real_uy, p_avoid, expansion_limit, eps)={
    my(
        onevec=valuationvec(G,1),
        qfminim_radius = norml2(real_uy[,1]),
        checked_elements=Map(),
        beta_found = 0,
        expansion_ctr = 1,
        scan_elts, scan_ctr
    );

    while(beta_found != 1 && expansion_ctr < expansion_limit,
        scan_elts = qfminim(real_uy~*real_uy, qfminim_radius+eps,,2)[3];
        scan_ctr = 1;
        expansion_ctr += 1;

        print(" - Number of qfminim increases ", expansion_ctr);
        \\print("shortList: ", scan_elts, "\n", coeff_LLL_uy*scan_elts); \\ get elements in terms of the integral basis now

        while(scan_ctr <= length(scan_elts),
            \\print("current-element ",scan_elts[,scan_ctr]);
            if(mapisdefined(checked_elements, scan_elts[,scan_ctr]),
                \\ Do nothing
            ,\\else
                mapput(checked_elements, scan_elts[,scan_ctr], 1);
                beta = coeff_LLL_uy*scan_elts[,scan_ctr];               \\ in terms of the integral basis
                shortest_vec = numerical_LLL_uy*scan_elts[,scan_ctr];
                new_ideal = idealdiv(G, coeff_y, beta);
                ideal_denom = get_ideal_denom(new_ideal);
                \\tempvar = cubescan_unequal(new_ideal, onevec, G,eps);if(length(tempvar) < 1,print("CUBESCAN: ", tempvar););

                if(gcd(p_avoid,ideal_denom) ==1 && length(cubescan_unequal(new_ideal, onevec, G,eps)) == 0,
                    print("Counters: ", scan_ctr, " ", expansion_ctr);
                    print(" - acceptable beta found via qfminim");
                    nu=abs(shortest_vec)~;
                    return([1, new_ideal, nu, beta]);
                );
            );
            scan_ctr += 1;
        );
        qfminim_radius = qfminim_radius*2;
    );
    if(beta_found == 0, ideal_denom = p_avoid);                         \\ if qfminim failed, ensure we enter neighbours

    return([0,0,0,0]);
}
