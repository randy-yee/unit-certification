\\
/*

COMPACT REPRESENTATION FUNCTIONS:

- idealsquare
- minima_from_LLL
- create_target

- reddiv_compact
-- find_coprime_divisor_qfminim
-- find_coprime_divisor_lazy
-- find_coprime_divisor_neighbors
-- get_LLL_basis_change

- giantstep
-- get_doubling_number
-- double_and_reduce
-- double_and_reduce_divisor_loop

- compact_reconstruct
- log_from_cpct
- intermediate_check

- compact_rep_buchmann
-- compute_initial_s
-- get_alpha_and_d
-- update_tracking_values
-- cpct_rep_final_collect

- cpct_from_loglattice
- cpct_denom_switch
- update_expmat
- build_unit
- cpct_rep_modp
- cpct_modp_from_factors
*/


\\ IMPORTS
read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/Neighbours.py");
DEBUG_CPCT=0;

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

GIANT STEP functions

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
/******************************************************************************/
/*THE GIANT STEP--RENE'S REDUCTION JUMP ALGORITHM: find a shortest vector of the lattice L = v*I */
/*Compute the square of a frac. ideal I = y where y is matrix with fractional entries
  Input: y is a matrix representing a fractional ideal,
    G is the numberfield produced by nfinit */

\\  Squares the ideal y in the number field G
\\ INPUT:
\\ - y an ideal in terms of integral basis of G
\\ - G a number field
idealsquare(y,G)={
    my( y_hnf=idealhnf(G,y) );
    my(y_square=idealpow(G,y_hnf,2));
    return(y_square);
};

\\subroutine of the compact_rep_buchmann function
\\ s_term is only length G.r1 + G.r2-1, and we need to obtain an extra coordinate
\\ After that we exponentiate each coordinate
create_target(G, s_term, inverse= 0)={
    my(exponentiated_s, new_coord=0);
    for(i=1, length(s_term),
        if(i <=G.r1, new_coord-=s_term[i], new_coord-=2*s_term[i]);
    );
    if( (length(s_term)+1) > G.r1, new_coord = new_coord/2);
    exponentiated_s = concat(s_term, new_coord);
    if(inverse ==1, exponentiated_s = -exponentiated_s);
    exponentiated_s = vector(length(exponentiated_s), i, exp(exponentiated_s[i]));
    return(exponentiated_s);
}

\\ Input: number field G and an ideal represented by columns of embeddings
\\ of the basis elements
\\ Output: the change of basis matrix that transforms the representation
\\ to an LLL reduced representation
get_LLL_basis_change(G, ideal_numerical)={
    real_embedded_mat = embed_real(G, ideal_numerical);
    LLL_change_of_basis = qflll(real_embedded_mat);
    return(LLL_change_of_basis);
}

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


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


/******************************************************************************/
/*REDUCE THE DIVISOR (y, u) where y is a fractional ideal and u is a totally
\\ positive vector in R^n
\\ Given a divisor (y,u), finds a reduced divisor (new_y, nu) close to it.
\\ Necessarily, new_y is a reduced ideal, nu corresponds to a minimum of new_y */
\\ INPUT:
\\  - y is a coefficient matrix of a fractional ideal,
\\  - u is a totally positive vector in R^(r1+r2). Later converted to a vector in R^n.
\\  - G is the number field
\\  - M1 is the (r1+r2 x n) matrix of conjugates of elements of the integral basis
\\  - p_avoid is optional, will force a choice of beta so that the ideal new_y has denom coprime to p_avoid */
\\ OUTPUT:
\\ - new_y, the coefficient matrix of a reduced ideal,
\\ - nu, corresponds to u*beta, where beta is a minimum such that new_y = (1/beta)*y is reduced
\\ - lmu is the distance from the original divisor i.e. y/u as compared to new_y/nu. equal to log(nu)-log(u) = log(beta)

/******************************************************************************/

reddiv_compact(y,u,G,M1, p_avoid=1)={

    my(y1, ideal_uY, numerical_mat_Y, red1, shortest_vec, nu, lmu,
        ideal_denom,vec_ctr,beta_found =0,
        comp = 2^(-ceil((poldegree(G.pol)^2 +2)*log(infinity_norm(u))+2*poldegree(G.pol)^2 +5)),
        EXPANSION_LIMIT = 9                                                     \\ used to limit the number of qfminim radius doublings
    );

    numerical_mat_Y = M1*y;                                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(numerical_mat_Y, u);                                      \\ ideal u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = y*LLL_change_of_basis;                                        \\ LLL basis, coords wrt the integral basis

    beta= LLLcoeffmat[,1];                                                      \\ beta holds coordinates of mu wrt the integral basis

    \\ need to scan to make sure the first basis vector is a shortest one
    real_mat_uY = embed_real(G, LLL_numerical_uY);
    enum_result = qfminim(real_mat_uY~*real_mat_uY, norml2(real_mat_uY[,1])-comp,,2 );

    /* NOTE THIS CHECK IS SLOW IN FIELDS WITH LARGE DEGREE (see pohst example)
    \\
    */
    if(length(enum_result[3])!=0 && !is_minimum(LLLcoeffmat,beta , G, comp),
        short_index =1;
        short_length = norml2(real_mat_uY[,1]);

        for(j=1, length(enum_result[3]),
            iter_length = norml2(real_mat_uY*enum_result[3][,j]);
            if(iter_length < short_length,
                short_index = j; short_length = iter_length;
            );
        );
        beta = LLLcoeffmat*enum_result[3][,short_index];

    );


    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);                      \\ the reduced ideal y / mu, in hnf form

    shortest_vec = LLL_numerical_uY[,1];                                        \\ vector of complex embeddings for the shortest vector of u*y, denoted mu

    nu=abs(shortest_vec)~;                                                      \\ nu is a t_VEC of dimension r, (complex coordinates are not squared)

    lmu = log(nu)-log(u);                                                       \\ expect equal to log(nfeltembed(G, beta))

    if(DEBUG_REDDIV && 0,
        \\print("REDDIV: shortvec=", precision(abs(shortest_vec),10));
        \\print("equal sv? ", precision(abs(pointwise_vector_mul(nfeltembed(G,beta), u)),10));
    );

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


\\ subalgorithm for giantstep aka jump algorithm.
\\ obtain the number of doublings needed to approximate the divisor
\\ v is a log vector, length r1+r2-1, disc, and field_degree are the discriminant and number field degree
\\ eps is the allowed error
get_doubling_number(v, disc, field_degree, eps)={
    my( maxv = max(normlp(v), abs( sumvec(v) )) );
    if(maxv > eps^2,
        return( max(0,floor(log(field_degree*maxv/log(disc))/log(2))+1) )       \\ the condition in Schoof alg 10.8 for t
    ,\\else
        return(0));
}


\\ This is the inner loop operation for the giantstep function
double_and_reduce(G, div_ideal, div_u,log_distance, avp =1)={
    my(u_square, log_mu_k,beta);
    div_ideal = idealsquare(div_ideal, G);
    u_square = pointwise_vector_mul(div_u,div_u)~;
    [div_ideal,div_u,log_mu_k,beta] = reddiv_compact(div_ideal, u_square, G, G[5][1],avp);
    \\psmiu=psmiu+2^(double_num-k)*log_mu_k;
    log_distance = 2*log_distance+log_mu_k;

    if(DEBUG_CPCT, print("log_beta = ", precision(log_mu_k,10)));
    return([div_ideal, div_u, log_distance,beta]);
}

\\ INPUTS:
\\ - G a number field
\\ - double_num is the number of times needed to double to get an approximate divisor (see get_doubling_number)
\\ - div_ideal and div_u are components of a divisor (ideal, real vector)
\\ - log_distance is the log distance of the initial divisor from 0
\\ - outputs a divisor close to D^(2^t), where t is the doubling number
double_and_reduce_divisor_loop(G, double_num, div_ideal, div_u, log_distance)={
    my(beta);
    for(k=1,double_num,
        [div_ideal, div_u, log_distance,beta]=double_and_reduce(G, div_ideal, div_u,log_distance);

    );
    return([div_ideal, div_u, log_distance]);
}


\\ INPUTS:
\\ - y is an ideal
\\ - v is a totally positive real vector of dimenension r+s-1
\\ - G is the number field
\\ - n is the degree of the field
\\ - eps is the acceptable error
\\ OUTPUT:
\\ - idealB a reduced ideal along with u
/******************************************************************************/
giantstep(y,v,G,n,eps)={
    if(type(v) == "t_COL", v = v~);
    my(
        disc = abs(G.disc),
        t = get_doubling_number(v, disc, n, eps),
        idealB, u, log_distance, beta,
        shrunken_v, shrunken_target
    );



    shrunken_v=2^(-t)*v;                                                        \\ this is y_sigma
    shrunken_target = create_target(G, shrunken_v, 1);
    [idealB,u,log_distance,beta]=reddiv_compact(y,shrunken_target,G,G[5][1]);   \\ (Ideal, minimum, log distance)

    [idealB, u, log_distance] = double_and_reduce_divisor_loop(G, t, idealB, u, log_distance);

    return( [idealB,u,vector_approximate(log_distance,eps)] );
}
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

giantstep_high_precision(y,v,G,n,eps, complex = 0)={
if(type(v) == "t_COL", v = v~);
    my(
        disc = abs(G.disc),
        t = get_doubling_number(v, disc, n, eps),
        idealB, u, log_distance, beta,
        shrunken_v, shrunken_target,
        complex_log
    );
    shrunken_v=2^(-t)*v;                                                        \\ this is y_sigma
    shrunken_target = create_target(G, shrunken_v, 1);
    [idealB,u,log_distance,beta]=reddiv_compact(y,shrunken_target,G,G[5][1]);   \\ (Ideal, minimum, log distance)

    if(complex ==1,
        complex_log = log(nfeltembed(G,beta));
    );
    for(k=1,t,
        shrunken_v*=2;
        idealB = idealsquare(idealB, G);

        u_square = pointwise_vector_mul(u,u)~;

        \\debug_compare(log(u_square), shrunken_v- 2*log_distance[1..G.r1+G.r2-1]);
        \\debug_compare(u_square, abs(create_target(G, shrunken_v- 2*log_distance[1..G.r1+G.r2-1], 1)) );
        [idealB,u,log_mu_k,beta] = reddiv_compact(idealB, u_square, G, G[5][1],avp);

        log_distance = 2*log_distance+log_mu_k;
        if(complex == 1, complex_log = 2*complex_log +log(nfeltembed(G,beta)) ;);
    );
    if(complex ==1,
        \\print("giant step outputs: \n", precision(log_distance,10), "   ", precision(complex_log,10));
        return([idealB,u,complex_log]); );
    return( [idealB,u,vector_approximate(log_distance,eps)] );
}

\\END GIANT STEP functions



\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Converts the compact representation form of a number field element to standard form
\\ INPUT:
\\ - G a number field
\\ - vec_alpha, a vector of number field elements
\\ - vec_d a vector of integers
\\ OUTPUT:
\\ - A number field element corresponding to \prod_{j=1}^k (alpha_j/d_j)^{2^{k-j}}
\\\\\\\\\\\\\\\\\\\\\\\\\
compact_reconstruct(G, vec_alpha, vec_d)={
my(alphafinal = 1, intermediate);
for(j=1, length(vec_alpha),
    intermediate = nfeltdiv(G, vec_alpha[j], vec_d[j]);
    \\print("************** Loop ", j, " **************\n  alpha/d: ", intermediate);

    intermediate = nfeltpow(G, intermediate, 2^(length(vec_alpha)-j));
    alphafinal = nfeltmul(G,alphafinal, intermediate);
    \\print("(alpha/d)^",2^(length(vec_alpha)-j), " ", intermediate, "\nintermediate product: ", alphafinal);
  );
  return(alphafinal);
}

complex_log_from_cpct(G, cpct_rep)={
    my(intermediate, lvec_complex = vector(G.r1+G.r2,k,0) ) ;
    for(j=1, length(cpct_rep[1]),
        lvec_complex*=2;
        intermediate = log(nfeltembed(G, nfeltdiv(G, cpct_rep[1][j], cpct_rep[2][j])));
        lvec_complex+= intermediate;
      );
      return(lvec_complex);
}

log_from_cpct(G, cpct_rep)={
    my(intermediate, lvec_complex = vector(G.r1+G.r2,k,0) ) ;
    for(j=1, length(cpct_rep[1]),
        lvec_complex*=2;
        intermediate = log(abs(nfeltembed(G, nfeltdiv(G, cpct_rep[1][j], cpct_rep[2][j]))));
        lvec_complex+= intermediate;
        \\print("************** Loop ", j, " **************\n  alpha/d: ", intermediate);
      );
      return(lvec_complex);
}
log_lattice_from_compact_set(G, cpct_reps)={
    my(logarithm_lattice, urank = G.r1+G.r2-1);
    logarithm_lattice = Mat(log_from_cpct(G, cpct_reps[1])[1..urank]~);
    for(i=2, length(cpct_reps),
        logarithm_lattice = matconcat([logarithm_lattice,log_from_cpct(G, cpct_reps[i])[1..urank]~ ] );
    );
    return(logarithm_lattice);
}
intermediate_check(G, alphavec, dvec, ideal1)={
    my(reconstruction, intermediate_norm, intermediate, finalnorm =1);
    for(j=1, length(alphavec),

        intermediate = nfeltdiv(G, alphavec[j], dvec[j]);
        \\print("alpha/d: ", intermediate);
        intermediate = nfeltnorm(G, intermediate);
        \\intermediate = intermediate^(2^(length(alphavec)-j));
        \\print("(alpha/d)^",2^(length(alphavec)-j), " ", intermediate);

        finalnorm = finalnorm*finalnorm*intermediate;
    );
    \\reconstruction = compact_reconstruct(G, alphavec, dvec); print(nfeltnorm(G,reconstruction));
    intermediate_norm = idealnorm(G, ideal1);
    return([intermediate_norm, finalnorm]);
}


\\ subalgorithm of compact representation to get the initial value of s_term,
\\ which is the input log embedding divided by an appropriate power of 2
\\ Also returns the value kprime, which determines the number of doublings
compute_initial_s(alpha, kbound)={
    my(s_term,
        kprime = 1,
        alpha_inf_norm = infinity_norm(alpha)/2
    );
    while(alpha_inf_norm > kbound,
      alpha_inf_norm = alpha_inf_norm/2;
      kprime += 1;
    );
    s_term = alpha/(2^kprime);
    if(DEBUG_CPCT >1,
        print("CPCT REP:\n -- kprime = ", kprime);
        print("\nW's boundary: ", precision(kbound,10), ". |alpha|= ",precision(abs(alpha), 10));
    );
    return([s_term,kprime]);
}

\\ subalgorithm of compact representation function for computing alpha and d
get_alpha_and_d(G, idealB, beta)={
    my(ideal_denom, alpha_i);
    ideal_denom = get_ideal_denom(idealB);
    alpha_i = nfeltdiv(G, ideal_denom, nfbasistoalg(G, beta));
    return([alpha_i, ideal_denom]);
}

\\ subalgorithm of compact representation function for updating the values s_term and desired_beta_log
update_tracking_values(s_term, log_rho)={
    my(double_s = 2*s_term);
    return([double_s, double_s- 2*log_rho]);
}

cpct_rep_final_collect(G, idealB, log_beta, desired_betalog, eps)={
    my(neighbours_output, checkvec, boundary_vector, ctr=1, unit_rank = G.r1+G.r2-1);
    if(DEBUG_CPCT >0, print("USING COLLECT\nFinal Beta is not equal to the desired value. Using NEIGHBOURS Function"););

    boundary_vector = vector(unit_rank, j, 4* sqrt(abs(G.disc)));
    boundary_vector =concat(Bvec,(2*4*abs(G.disc))^(G.r1+G.r2/2));

    neighbours_output = COLLECT(G, idealB, boundary_vector, eps);               \\ is a list of distinct neighbors of 1, as column vectors
    if(DEBUG_CPCT >0, print("OUTPUT of COLLECT ", neighbours_output););

    while(!samevecs(desired_betalog, checkvec, eps),
        checkvec = vector(length(desired_betalog),j, log(abs( nfeltembed(G,neighbours_output[ctr]) ))[j]);
        if(DEBUG_CPCT >0, print(ctr ," minima ", precision(checkvec,10)););
        checkvec += log_beta[1..unit_rank];

        ctr += 1;
        if(ctr > length(neighbours_output), print("No neighbours of 1 satisfy the condition. Error."); return(-1));
    );

    idealB = idealdiv(G, idealB, neighbours_output[ctr]);
    beta = nfeltmul(G, beta, neighbours_output[ctr]);
    log_beta = checkvec;

    return([idealB, log_beta, beta]);
}
/******************************************************************************/
\\ INPUT:
\\ - alpha should be a row vector of the log embedding of a unit
\\ - alphaOK is the Z basis for principal ideal generated by alpha (O_K if alpha a unit)
\\ - G a number field
\\ - eps some error
\\ - avp a prime for which we want all denominators to be coprime to
\\ - if arith_prog is nonzero, then makes it so if a denom >=2, then it is not 1 mod p.
\\ - see Thiel's description of compact representation
/******************************************************************************/
compact_rep_buchmann(G, alpha, alphaOK , eps, avp=1)={

  my(
    unit_rank = G.r1 +G.r2 -1,
    kprime = 1,                                               \\ determines no. steps the algorithm runs for
    kbound,                                                   \\ holds the boundary of the area W
    idealB = alphaOK,                                         \\ variable for the ideal B
    beta = 1,                                                 \\ holds the beta
    d_vec = [1],                                              \\ holds the d_i values
    ideal_denom = 1,                                          \\ used to determine d_i
    alpha_vec= List([1]),                                     \\ holds the alpha_i
    alpha_i = 1,                                              \\ used to determine alpha_i
    log_rho = vector(length(alpha), j, 0),                    \\ holds log of rho_i
    s_term,
    target,
    desired_betalog,
    exp_s
  );
  if(type(alpha) == "t_COL", alpha = alpha~);
  kbound = log(sqrt( abs(G.disc) ))/log(2)/2;                 \\ defines the boundary of the area W

  \\ MAIN LOOP: following the algorithm of Thiel, under the assumption alpha is a unit.
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\ Note: In this case, gamma = 1, and beta1 = 1 also. This means we can skip iteration 1

  [s_term, kprime] = compute_initial_s(alpha, kbound);                                     \\ This is technically the value -s1, since Log(rho) = -alpha
  exp_s = create_target(G, 2*s_term);
  [idealB, target, log_beta, beta] = reddiv_compact(idealB, exp_s, G, G[5][1],avp);        \\ reddiv_compact computes A_2, target, Log(beta_2), beta2

  s_term = -2*s_term;                                                                      \\ this is s2
  log_rho = log_beta;                                                                      \\ rho2 = beta2

  [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta);
  d_vec = concat(d_vec,ideal_denom);                                                       \\ Get d_i= d(A) and append to tracking vector
  alpha_vec = concat(alpha_vec, alpha_i);                                                  \\ compute alpha_i = d/beta and append to tracking vector

  if(DEBUG_CPCT >0,
    print("log beta and s: ", precision(log_beta,10), "   ", precision(s_term,10));
    print(" - ROUND ", 2, ": beta close enough?: ", check_closeness(log_beta[1..unit_rank], s_term, kbound), ": NORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
    );

  \\ TRIVIAL CASE HANDLING
  if (kprime == 1, return([ alpha_vec,d_vec ]));

  for(i=3, kprime,
      [s_term, desired_betalog] = update_tracking_values(s_term, log_rho[1..unit_rank]);  \\s_i and the target beta
      [idealB, target, log_rho, beta] = double_and_reduce(G, idealB, target,log_rho,avp);

      log_beta= log(abs(G[5][1]*beta));

      [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta);
      d_vec = concat(d_vec,ideal_denom);
      alpha_vec = concat(alpha_vec, alpha_i);

      if(DEBUG_CPCT,
        print("logbeta and s: ", precision(log_beta,10), "   ", precision(s_term,10));
        print(" - ROUND ", i, ": beta close enough?: ", check_closeness(log_beta[1..unit_rank], desired_betalog, kbound), ".\nNORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
      );

  ); \\ end main for loop

    \\\  ENTER THE FINAL ITERATION OF THE ALGORITHM, WHICH IS SPECIAL
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    s_term = 2*s_term;                                                          \\ this is s_{k'+1} = log rho
    desired_betalog = (-alpha) - (2*log_rho[1..unit_rank]);                     \\ this equation is "rho - 2*rho_{k'}"
    [idealB, target, log_rho, beta] = double_and_reduce(G, idealB, target,log_rho);
    log_beta= log(abs(G[5][1]*beta))~;

    if(DEBUG_CPCT >0,
        print(" - PRIOR TO COLLECT:\n -  rho_i ", precision(log_rho,10), " -  alpha ", precision(-alpha, 10));
        debug_print(" -  log of needed minimum ",-alpha-log_rho[1..unit_rank]);
    );

    \\ LAGRANGE PART: idealB is reduced and we have a minimum beta which may or may not have the correct log vector
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(!samevecs(abs(desired_betalog), abs(log_beta[1..unit_rank]),eps),
        [idealB, log_beta, beta] = cpct_rep_final_collect(G, idealB, log_beta, desired_betalog, eps);
    );

    [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta);
    d_vec = concat(d_vec,ideal_denom);
    alpha_vec = concat(alpha_vec, alpha_i);

    if(DEBUG_CPCT >0, print(" - ROUND ", kprime+1, ": NORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB)););
    if(DEBUG_CPCT, print(" -- final rho_i ", precision(2*log_rho[1..length(log_beta)] + log_beta,10)););

    return( [alpha_vec, d_vec] );
} \\ end compact_rep_buchmann





\\ INPUT:
\\ - lglat a log lattice corresponding to an independent set of units
\\ - G a number field,
\\ - eps some precision value
\\ OUTPUT:
\\ - a vector of compact representations for each of the units in lglat
cpct_from_loglattice(G, lglat, eps, avp=1)={
    my(c_rep_list = [],
        aOK = matid(poldegree(G.pol)),
        cr_i = [];
    );
    for (i=1, length(lglat),
        cr_i = compact_rep_buchmann(G,lglat[,i]~,aOK, eps, avp);
        c_rep_list = concat(c_rep_list, [cr_i]);
    );
    return(c_rep_list);
}

\\ Given a list of cpct reps and a prime p, check each to make sure the denominators are coprime to p. If not, compute a new representation
\\ INPUT:
\\ - G a number field
\\ - lglat a log lattice corresponding to an independent set of units
\\ - cpct_list a corresponding list of compact reps of the units
\\ - eps some precision value
\\ - p a prime number
\\ OUTPUT:
\\ - a vector of compact representations for each of the units in lglat whose denominators are all coprime to p
cpct_denom_switch(G, lglat, cpct_list, eps, p)={
    my(dlist, dflag);
    for(i = 1, length(cpct_list),
        dlist = cpct_list[i][2];
        dflag = 1;
        for(j =1, length(dlist),
            if(gcd(dlist[j],p)!= 1, dflag = 0);
        );
        if(dflag ==0,
            print("Nontrivial denominator switch: ", p); breakpoint();
            cpct_list[i] = compact_rep_buchmann(G, lglat[,i]~, matid(poldegree(G.pol)), eps, p));
    );
    return(cpct_list);
};

\\ INPUT:
\\ - lglat a log lattice corresponding to an independent set of units
\\ - G a number field,
\\ - eps some precision value
\\ OUTPUT:
\\ - a vector of compact representations for each of the units in lglat
update_expmat(eta_exp_mat, updatevec,k,p )={
    for(j = k+1, length(eta_exp_mat),
        eta_exp_mat[,j] = eta_exp_mat[,j] + updatevec[j-k]*eta_exp_mat[,k];         \\ add the correct multiple of the kth column
    );

    eta_exp_mat = eta_exp_mat %p ;                                  \\ reduce coeffs mod p
    \\ ensure exponents lie in the range [-p/2, p/2]
    for(bb =1, matsize(eta_exp_mat)[1],
        for(aa =1, matsize(eta_exp_mat)[2],
          if( eta_exp_mat[aa,bb] > p/2, eta_exp_mat[aa,bb] = eta_exp_mat[aa,bb] -p );
    ));
    return(eta_exp_mat);
}

\\ INPUT:
\\ - units: a log matrix of units
\\ - exponent_mat: an matrix of exponents,
\\ - k the index of the element we're trying to compute
\\ OUTPUT:
\\ - a column of exp mat yields a linear combination
\\   of the columns of the log matrix
build_unit(units, exponent_mat,k)={
    my(test_eta_k = exponent_mat[1,k]*units[,1]);
    for(i=2, length(units),
        test_eta_k = test_eta_k + exponent_mat[i,k]*units[,i];
    );
    return(test_eta_k);
}


\\ INPUT:
\\ - G a number field
\\ - res_field a residue field generated from a prime ideal
\\ - cpct is a compact representation, a 2 elt array, where the first elt is the set of alpha, the second is the set of d
\\   such that the original element is given by the product over i of (alpha[i]/d[i])^{2^{k-i}}
\\ OUTPUT:
\\ - the embedding of the cpct element into the specified residue field.
cpct_rep_modp(G,cpct, res_field)={
    my(alphafinal = 1, intermediate, twopow,
        vec_d = cpct[2],
        vec_alpha = cpct[1]
    );

    for(i=1, length(vec_d), if(gcd(vec_d[i], res_field[3][1])!=1, print("ERROR: CPCT REP is not suitable to be reduced mod the ideal"); return(-1); ));

    twopow = 2^(length(vec_alpha)-1);

    for(j=1, length(vec_alpha),
        intermediate = nfeltdiv(G, vec_alpha[j], vec_d[j]);
        intermediate = nfmodpr(G, intermediate, res_field);
        intermediate = intermediate^twopow;
        alphafinal = alphafinal*intermediate;
        twopow = twopow/2;
    );
    return(alphafinal);
}

\\ Obtain the compact representation of an element which is given in terms of an exponent vector of length r, and r compact representations
cpct_modp_from_factors(G, cpct_list, expvec, res_field)={
    my(finalproduct, intermediate);
    finalproduct = nfmodpr(G, 1, res_field);
    for(i =1, length(cpct_list),
        intermediate = cpct_rep_modp(G,cpct_list[i], res_field);                 \\ this should embed the ith cpct rep from supplied list
        intermediate = intermediate^expvec[i];                                   \\ raise it to the ith exponent from the list
        finalproduct *= intermediate;
    );
    return(finalproduct);
}

compute_subgroup(G, unimat, modpair, extype=0)={
    my(coord = modpair[1], pow = modpair[2]);
    unimat[,coord] = nfeltpow(G, unimat[,coord], pow);
    if(extype ==0,
        return(unimat);
    );
    coord2 = modpair[3]; pow = modpair[4];
    if(extype ==1,
        unimat[,coord2] = nfeltpow(G, unimat[,coord2], pow);
        return(unimat);
    );
    if(extype ==2,
        unimat[,coord] = nfeltmul(G, unimat[,coord], nfeltpow(G, unimat[,coord2], pow));
        return(unimat);
    );
}
