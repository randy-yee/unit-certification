\\
/*

COMPACT REPRESENTATION FUNCTIONS:

- idealsquare:
- create_target
- exponentiate_logvec
- extra_log_coordinate
- get_LLL_basis_change

- reddiv_compact
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


\\ # Squares the ideal y in the number field G
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
    new_coord = extra_log_coordinate(G.r1, G.r2, s_term);

    exponentiated_s = concat(s_term, new_coord);
    if(inverse ==1, exponentiated_s = -exponentiated_s);
    exponentiated_s = unsquare_log_embeddings(G, exponentiated_s);
    exponentiated_s = vector(length(exponentiated_s), i, exp(exponentiated_s[i]));
    return(exponentiated_s);
}

\\ assuming we have r+1 terms already
\\ r - the unit rank
\\ s_term - a vector of length (r+1) corresponding to log vector of an element
\\ inverse - flag indicating whether to compute the inverse element, default 0
exponentiate_logvec(r, s_term, inverse = 0)={
    my(exponentiated_s, a = (-1)^inverse);
    GP_ASSERT_TRUE(length(s_term) == r+1);
    exponentiated_s = vector(r+1, i, exp(a*s_term[i]));
    return(exponentiated_s);
}

\\ function to compute the extra r+1 coordinate, under the assumption that
\\ s_term is the logarithm vector of a unit (with factors of 2 on the complex coords)
extra_log_coordinate(r1, r2, s_term)={
    my(new_coord = 0);
    GP_ASSERT_TRUE(length(s_term) == r1+r2-1);
    for(i=1, length(s_term),
        new_coord-=s_term[i];
    );
    return(new_coord);
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
/******************************************************************************/
/*THE GIANT STEP--RENE'S REDUCTION JUMP ALGORITHM: find a shortest vector of the lattice L = v*I */

/******************************************************************************/
/*REDUCE THE DIVISOR (y, u) where y is a fractional ideal and u is a totally
\\ positive vector in R^(r+1)
\\ Given a divisor (y,u), finds a reduced divisor (new_y, nu) close to it.
\\ new_y is a reduced ideal, nu corresponds to a minimum beta of y such that
\\ new_y = y/beta
*/
\\ INPUT:
\\  - y is a coefficient matrix of a fractional ideal,
\\  - u is a totally positive vector in R^(r1+r2).
\\  - G is the number field
\\  - M1 is the (r1+r2 x n) conjugate embedding matrix of the integral basis
\\  - p_avoid optional, DEPRECATED originally will forced a choice of beta so that the ideal new_y has denom coprime to p_avoid */
\\ OUTPUT:
\\ - new_y, the coefficient matrix of a reduced ideal,
\\ - nu, corresponds to u*beta, where beta is a minimum such that new_y = (1/beta)*y is reduced
\\ - lmu is equal to log(nu)-log(u) = log(beta). In some sense meaures distance btwn (y,u) as compared to (new_y,nu).

/******************************************************************************/
\\
reddiv_compact_invert_u(~y,~u,~G,~M1, p_avoid=1)={
    my(ideal_uY, numerical_mat_Y, red1, shortest_vec, nu, lmu,
        ideal_denom,vec_ctr,
        comp = 2^(-ceil((poldegree(G.pol)^2 +2)*log(infinity_norm(u))+2*poldegree(G.pol)^2 +5))
    );
    \\ inverts so that we compute the lattice u^(-1)I
    u = vector(length(u), i, 1/u[i]);
    numerical_mat_Y = M1*y;                                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(~numerical_mat_Y, ~u);                                    \\ ideal u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = y*LLL_change_of_basis;                                        \\ LLL basis, coords wrt the integral basis

    \\ need to scan to make sure the first basis vector is a shortest one
    real_mat_uY = embed_real(~G, ~LLL_numerical_uY);
    \\print("norm comparisons: ", precision(norml2(ideal_uY*[1,0,0]~),10), "  ", precision(norml2(ideal_uY*[-35,0,1]~),10));
    true_shortest = qfminim(real_mat_uY~*real_mat_uY,,,2);
    \\print("1 ", true_shortest[3][,1], "  ",precision(norml2(LLL_numerical_uY*true_shortest[3][,1]),10));
    beta = LLLcoeffmat*true_shortest[3][,1];
    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);                      \\ the reduced ideal y / mu, in hnf form

    nu = pointwise_vector_mul(abs(M1*beta),u)~;
    \\nu = pointwise_vector_mul(abs(get_normalized_embedding_vector(G,beta)),u)~;
    lmu = log(nu)-log(u);                           \\ this is log of beta
    lmu = double_complex_coordinates(G.r1, lmu);

    \\ # this comparison doesn't work anymore. need to double log coords
    \\lmu = -log(abs(LLL_numerical_uY*true_shortest[3][,1]))~;
    \\print("compare: ", precision(log(abs(nfeltembed(G,beta))),10), "  ", precision(lmu,10), "\n");

    \\nu*=idealnorm(G, new_y)^(-1/poldegree(G.pol));   \\ account for ideal norm (see schoof formula for v)
    nu = vector(length(nu), i , 1/nu[i]);           \\ embedding of u/beta

    return([new_y, nu, lmu, beta]);
}

reddiv_compact(~y,~u,~G,~M1, p_avoid=1)={
    my(ideal_uY, numerical_mat_Y, red1, shortest_vec, nu, lmu,
        ideal_denom,vec_ctr,
        comp = 2^(-ceil((poldegree(G.pol)^2 +2)*log(infinity_norm(u))+2*poldegree(G.pol)^2 +5))
    );

    numerical_mat_Y = M1*y;                                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(~numerical_mat_Y, ~u);                                    \\ lattice u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = y*LLL_change_of_basis;                                        \\ LLL basis, coords wrt the integral basis

    \\# enumerate to find the shortest vector
    real_mat_uY = embed_real(~G, ~LLL_numerical_uY);
    true_shortest = qfminim(real_mat_uY~*real_mat_uY,,,2);
    beta = LLLcoeffmat*true_shortest[3][,1];

    \\# compute the ideal (1/beta)y
    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);

    nu = pointwise_vector_mul(abs(M1*beta),u)~;
    lmu = log(nu)-log(u);
    lmu = double_complex_coordinates(G.r1, lmu);

    \\print("compare: ", precision(log(abs(nfeltembed(G,beta))),10), "  ", precision(unsquare_log_embeddings(G,lmu),10), "\n");

    return([new_y, nu, lmu, beta]);
}



\\\ # The reduction algorithm as described in Jacobson-Tran-Yee
\\ Input:
\\ - nf_K is the number field
\\ - phi_Omega is the precomputed Minkowski matrix for the integral basis
\\ - idealI is the input ideal, given as a coefficient matrix wrt integral basis
reduction(~nf_K, ~phi_Omega, ~idealI, ~u)=
{
    my(
        nf_K_degree = poldegree(nf_K.pol),
        minkowski_I,
        lattice_uI,
        LLL_transform_mat,
        minkowski_uI_LLL,
        transformed_I,
        real_lattice_uI,
        shortest_vec, beta, nu, lmu, idealJ,
        comp = 2^(-ceil((nf_K_degree^2 +2)*log(infinity_norm(u))+2*nf_K_degree^2 +5))
    );

    minkowski_I = phi_Omega*idealI;                                             \\ minkowski embedding matrix of idealI
    lattice_uI = mulvec(~minkowski_I, ~u);                                      \\ the lattice u*idealI in minkowski space
    LLL_transform_mat = get_LLL_basis_change(nf_K, lattice_uI);               \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    minkowski_uI_LLL = lattice_uI*LLL_transform_mat;                          \\ Obtain LLL basis of u*idealI in numerical form (possibly complex)
    transformed_I = idealI*LLL_transform_mat;                                   \\ LLL basis, coords wrt the integral basis

    \\ # Embed to real space (R^n) and use enumeration to ensure shortest vector
    real_lattice_uI = embed_real(~nf_K, ~minkowski_uI_LLL);
    shortest_vec = qfminim(real_lattice_uI~*real_lattice_uI,,,2);

    \\ # compute coordinates of minimum beta (c_beta) and the reduced ideal I/beta
    beta = transformed_I*shortest_vec[3][,1];
    idealJ = idealdiv(nf_K, idealI, beta);

    \\ # compute the distance from the original target vector
    nu = pointwise_vector_mul(abs(phi_Omega*beta),u)~;
    lmu = log(nu)-log(u);

    \\GP_ASSERT_VEC_NEAR(nu,abs(shortest_vec), comp  );
    \\print(precision(log(abs(nfeltembed(nf_K,beta))),10), "  ", precision(lmu,10), "\n");

    lmu = double_complex_coordinates(nf_K.r1, lmu);

    return([idealJ, nu, lmu, beta]);
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

    \\print("REMOVE!!!");
    [div_ideal,div_u,log_mu_k,beta] = reddiv_compact(div_ideal, u_square, G, G[5][1],avp);

    \\[div_ideal,div_u,log_mu_k,beta] = reddiv_compact_invert_u(div_ideal, u_square, G, G[5][1],avp);

   \\ print("INVERTED ", precision([div_ideal,div_u,log_mu_k,beta],10));
    \\print("REDDIV 2" , precision(tempvec,10));
    /*
    if (beta != tempvec[4],
        beta1 = beta;
        logbeta1 = log(abs(nfeltembed(G, beta1)));
        beta2 = tempvec[4]; \\ non inverted u
        logbeta2 = log(abs(nfeltembed(G, beta2)));
        print("    inverted distance: l_inf: ", precision(normlp(log(u_square)-logbeta1),10), " l_2: ", precision(norml2(log(u_square)-logbeta1),10));
        print("non-inverted distance: l_inf: ", precision(normlp(log(u_square)-logbeta2),10),  " l_2: ", precision(norml2(log(u_square)-logbeta2),10));
        breakpoint();
        \\ need to figure out why log nu and these outputs are different.
        \\ remember that you altered reddiv_compact and double reduce and compact_rep_full
    );
    */
    log_distance = 2*log_distance+log_mu_k;
    return([div_ideal, div_u, log_distance,beta]);
}

\\ This is the inner loop operation for the giantstep function
double_and_reduce_invert(G, div_ideal, div_u,log_distance, avp =1)={
    my(u_square, log_mu_k,beta);
    div_ideal = idealsquare(div_ideal, G);
    u_square = pointwise_vector_mul(div_u,div_u)~;

    [div_ideal,div_u,log_mu_k,beta] = reddiv_compact_invert_u(div_ideal, u_square, G, G[5][1],avp);
    log_distance = 2*log_distance+log_mu_k;

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
    shrunken_target = create_target(G, shrunken_v, 1);                          \\ adds the (r+s)th coord, so that the sum of the log embeddings is 0

    \\print(length(shrunken_v), "  ", G.r1+G.r2);
    print("WARNING: Giant step is being deprecated. Please ensure appropriate use");

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
    for(k=1, t,
        shrunken_v*=2;
        idealB = idealsquare(idealB, G);

        u_square = pointwise_vector_mul(u,u)~;

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

\\ INPUTS:
\\ - y is an ideal
\\ - v is a totally positive real vector of dimenension r+s-1
\\ - G is the number field
\\ - n is the degree of the field
\\ - eps is the acceptable error
\\ - complex is a flag for if you want complex logs instead.
\\ OUTPUT:
\\ - a triplet consisting of a reduced ideal J and a minimum beta satisfying
\\ - J = (1/beta)*y, and also u which indicates the distance from the vector v
\\ Here beta is a compact representation

/******************************************************************************/
jump_compact(y,v,G,n,eps, complex = 0)={
if(type(v) == "t_COL", v = v~);
    my(
        disc = abs(G.disc),
        new_t, new_k = 1,
        idealB, u, log_distance, beta,
        shrunken_v, shrunken_target,
        complex_log,
        alpha_vec= List([1]),                               \\ holds alpha_i
        d_vec = [1],                                        \\ holds d_i
        alpha_i, ideal_denom
    );
    shrunken_v = unsquare_log_embeddings(G, v);
    parF = ((2/Pi)^(G.r2))*sqrt(abs(G.disc));               \\\ used to define the shrinking bound
    new_k = floor(log( normlp(shrunken_v)*poldegree(G.pol)/ log(parF)) / log(2) )+1;
    new_t = max(0, new_k);
    shrunken_v=2^(-new_t)*shrunken_v;                                                        \\ this is y_sigma

    if(length(v) == G.r1+G.r2,
        shrunken_target = exponentiate_logvec(G.r1+G.r2-1, shrunken_v, 1);
        \\shrunken_target = exponentiate_logvec(G.r1+G.r2-1, shrunken_v);
    ,
        \\shrunken_target = create_target(G, shrunken_v, 1);
        shrunken_target = create_target(G, shrunken_v);
        print("Jump was only given an r-vec, make sure it's a unit!");
    );

    \\[idealB,u,log_distance,beta]=reddiv_compact_invert_u(y,shrunken_target,G,G[5][1]);
    [idealB,u,log_distance,beta]=reddiv_compact(y,shrunken_target,G,G[5][1]);

    \\# Get denominator of the new ideal and alpha_i = d/beta;
    [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 0);
    d_vec = concat(d_vec,ideal_denom);
    alpha_vec = concat(alpha_vec, alpha_i);

    if(complex ==1,
        breakpoint();
        complex_log = log(nfeltembed(G,beta));
    );
    for(k=1, new_t,
        shrunken_v*=2;
        idealB = idealsquare(idealB, G);

        u_square = pointwise_vector_mul(u,u)~;

        \\[idealB,u,log_mu_k,beta] = reddiv_compact_invert_u(idealB, u_square, G, G[5][1],avp);
        [idealB,u,log_mu_k,beta] = reddiv_compact(idealB, u_square, G, G[5][1],avp);

        GP_ASSERT_VEC_NEAR(log_mu_k,get_normalized_log_vector(G, beta),10);
        [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 0);
        d_vec = concat(d_vec,ideal_denom);
        alpha_vec = concat(alpha_vec, alpha_i);

        log_distance = 2*log_distance+log_mu_k;
        if(complex == 1, complex_log = 2*complex_log +log(nfeltembed(G,beta)) ;);
    );
    if(complex ==1,
        \\print("giant step outputs: \n", precision(log_distance,10), "   ", precision(complex_log,10));
        return([idealB,u,complex_log]);
    );
    myCompactRep = [alpha_vec, d_vec];
    \\debug_print("A ", log_from_cpct(G, myCompactRep)); debug_print("B ", log_distance);
    return( [idealB, u, myCompactRep] );
}

\\ Input: G - a number field, and a compact cpct_rep
\\ Return: The inverse element
invert_compact(~G, ~cpct_rep)=
{
    nf_argcheck(G);
    my(n_terms = length(cpct_rep[1]),
        inverse_elt,
        denom_lcm,
        inverse_rep = [vector(n_terms,t, 0),vector(n_terms,t, 0)]);
    for(k=1, length(cpct_rep[1]),
        inverse_elt = nfeltdiv(G, cpct_rep[2][k], cpct_rep[1][k]);
        inverse_rep[1][k] = numerator(inverse_elt);
        inverse_rep[2][k] = denominator(inverse_elt);
    );
    inverse_rep[1] = List(inverse_rep[1]);
    return(inverse_rep)
}

\\ Input: a number field, two compact representations
\\ return the product as a compact representation
\\ Function does not check for redundant terms (leading powers of 1)
mul_compact(~G, ~cpct1, ~cpct2)=
{
    nf_argcheck(G);
    cpct_rep_argcheck(cpct1);
    cpct_rep_argcheck(cpct2);

    my(num_terms, len1, len2, product_numerator=List(), product_denominator=[],
        numer, denom, common, numBound, denomBound);
    numBound = G.disc^(3/4*(G.r1+G.r2+2));
    denomBound = G.disc^(1/2);

    if( length(cpct1[1]) >= length(cpct2[1]),
        len1 = length(cpct1[1]);
        len2 = length(cpct2[1]);
        product_numerator=List(cpct1[1][1..(len1-len2)]);
        product_denominator=cpct1[2][1..(len1-len2)];
    ,
        len1 = length(cpct2[1]);
        len2 = length(cpct1[1]);
        product_numerator=List(cpct2[1][1..(len1-len2)]);
        product_denominator=cpct2[2][1..(len1-len2)];
    );
    for(i = 1, len2,
        numer = nfeltmul(G, cpct1[1][length(cpct1[1])-len2+i ], cpct2[1][length(cpct2[1])-len2+i ]);
        denom = cpct1[2][length(cpct1[1])-len2+i ]*cpct2[2][length(cpct2[1])-len2+i];
        common = gcd(nfbasistoalg(G,numer), denom);
        listput(~product_numerator, numer/common);
        product_denominator = concat(product_denominator, denom/common);
    );

    \\recompactify(G, product_numerator, product_denominator);
    \\print("Expensive verification in mul_compact. Delete once confirmed");
    \\GP_ASSERT_TRUE(compact_reconstruct(G, product_numerator, product_denominator) ==
    \\nfeltmul(G, compact_reconstruct(G, cpct1[1], cpct1[2]), compact_reconstruct(G, cpct2[1], cpct2[2])));

    return([product_numerator, product_denominator]);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ check if terms of compact representation are small enough and shrink them
\\ if needed
\\ INPUT:
\\ - G a number field
\\ - numerator vector of compact representation
\\ - denominator vector of compact representations
\\ OUTPUT:
\\ - A compact representation where all the terms are small
\\\\\\\\\\\\\\\\\\\\\\\\\
recompactify(G, cpctNum, cpctDenom)={
    GP_ASSERT_EQ(length(cpctNum), length(cpctDenom));
    print("Not written");breakpoint();
}
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Mainly use for testing and debugging.
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
        intermediate = nfeltpow(G, intermediate, 2^(length(vec_alpha)-j));
        alphafinal = nfeltmul(G,alphafinal, intermediate);
        \\print("(alpha/d)^",2^(length(vec_alpha)-j), " ", intermediate, "\nintermediate product: ", alphafinal);
    );
    return(alphafinal);
}

\\ Input: Number field G and a compact representation
\\ compute the complex logarithm vector (no absolute values on the minkowski embedding)
complex_log_from_cpct(G, cpct_rep)={
    my(intermediate, lvec_complex = vector(G.r1+G.r2,k,0) ) ;
    for(j=1, length(cpct_rep[1]),
        lvec_complex*=2;
        print("embedding");
        intermediate = log(nfeltembed(G, nfeltdiv(G, cpct_rep[1][j], cpct_rep[2][j])));
        print(precision(intermediate,10));
        print("done embedding");
        lvec_complex+= intermediate;
      );
      return(double_complex_coordinates(G.r1,lvec_complex));
}

\\ Input: Number field G and a compact representation
\\ compute the logarithm vector of absolute values
\\ doubling of the complex coords
log_from_cpct(G, cpct_rep)={
    my(intermediate, log_vector = vector(G.r1+G.r2,k,0) ) ;
    nf_argcheck(G);
    cpct_rep_argcheck(cpct_rep);
    for(j=1, length(cpct_rep[1]),
        log_vector*=2;
        intermediate = log(abs(nfeltembed(G, nfeltdiv(G, cpct_rep[1][j], cpct_rep[2][j]))));
        log_vector+= intermediate;
      );
      return(double_complex_coordinates(G.r1, log_vector));
}

\\ Given a list of compact representations
\\ Compute the corresponding log-vector matrix
log_lattice_from_compact_set(G, cpct_reps)={
    my(logarithm_lattice, urank = G.r1+G.r2-1);
    logarithm_lattice = Mat(log_from_cpct(G, cpct_reps[1])[1..urank]~);
    for(i=2, length(cpct_reps),
        logarithm_lattice = matconcat([logarithm_lattice,log_from_cpct(G, cpct_reps[i])[1..urank]~ ] );
    );
    return(logarithm_lattice);
}

\\ strictly a test function used to verify ideal generators within Jump and
\\ compact representation algorithm
intermediate_check(G, alphavec, dvec, ideal1)={
    my(reconstruction, intermediate_norm, intermediate, finalnorm =1);
    for(j=1, length(alphavec),

        intermediate = nfeltdiv(G, alphavec[j], dvec[j]);   \\#alpha/d
        intermediate = nfeltnorm(G, intermediate);
        \\intermediate = intermediate^(2^(length(alphavec)-j));
        \\#print("(alpha/d)^",2^(length(alphavec)-j), " ", intermediate);

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
    print("Deprecated. See compact_rep_full_input"); breakpoint();
    return();
}

\\ subalgorithm of compact representation function for computing alpha and d
get_alpha_and_d(G, idealB, beta, inverseflag)={
    GP_ASSERT_TRUE(inverseflag == 0 || inverseflag == 1);
    if(inverseflag == 1, "Please ensure using the inverse is correct here");
    my(ideal_denom, alpha_i);
    ideal_denom = get_ideal_denom(idealB);
    if(!inverseflag,
    alpha_i = nfeltmul(G, ideal_denom, nfbasistoalg(G, beta));,
    alpha_i = nfeltdiv(G, ideal_denom, nfbasistoalg(G, beta)););
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

    \\new_eps = ceil(log(abs(G.disc))/(2*n) + 3/2)+1;
    \\new_eps = 2^(-new_eps);
    if(length(log_beta) == G.r1+G.r2,
        boundary_vector = abs(desired_betalog - log_beta);
        bv2 = vector(unit_rank, j, 4* sqrt(abs(G.disc)));
        bv2 =concat(bv2,(4*abs(G.disc))^((G.r1+G.r2)/2));
        print("search bound: ", precision(boundary_vector,10), "\n", precision(bv2,10));
    ,
        print("WARNING: Using generic boundary vector in cpct rep final step");
        boundary_vector = vector(unit_rank, j, 4* sqrt(abs(G.disc)));
        boundary_vector =concat(boundary_vector,(4*abs(G.disc))^((G.r1+G.r2)/2));
    );

    neighbours_output = COLLECT(G, idealB, boundary_vector, eps);               \\ is a list of distinct neighbors of 1, as column vectors
    if(DEBUG_CPCT >0, print("OUTPUT of COLLECT ", neighbours_output););

    if(length(neighbours_output) == 0, print("No neighbours of 1 satisfy the condition. Error."); return(-1));
    checkvec = vector(length(desired_betalog),j, log(abs( nfeltembed(G,neighbours_output[ctr]) ))[j]);
    if(DEBUG_CPCT >0, print(ctr ," minima ", precision(checkvec,10)););

    checkvec += log_beta[1..unit_rank];
    while(!samevecs(desired_betalog, checkvec, eps),
        ctr += 1;

        if(ctr > length(neighbours_output), print("No neighbours of 1 satisfy the condition. Error."); return(-1));
        print("target: ", precision(desired_betalog, 10), "  ",
            precision(log_beta+log(abs( nfeltembed(G,neighbours_output[ctr]))), 10), "  ",
            precision(log(abs( nfeltembed(G,neighbours_output[ctr]) ))));
        checkvec = vector(length(desired_betalog),j, log(abs( nfeltembed(G,neighbours_output[ctr]) ))[j]);
        if(1, print(ctr ," minima ", precision(checkvec,10)););
        checkvec += log_beta[1..unit_rank];
    );

    idealB = idealdiv(G, idealB, neighbours_output[ctr]);
    beta = nfeltmul(G, beta, neighbours_output[ctr]);
    log_beta = checkvec;

    return([idealB, log_beta, beta]);
}

cpct_rep_final_enum(G, idealB, beta, log_beta, desired_betalog, alphaOK, eps, testFlag = 0)={

    my(DEBUG_CPCT_ENUM = 0,
        unit_rank = G.r1+G.r2-1,
        degree = poldegree(G.pol),
        new_eps);

    if (type(testFlag) == "t_INT" && (testFlag == 1),
        DEBUG_CPCT_ENUM = testFlag;
        print("Debugging CPCT Enumeration step");
    );

    if (DEBUG_CPCT_ENUM, print("final enum target: ", precision(desired_betalog, 10)));

    \\#see thesis Prop 3.4.3. The addition 5 bits of precision is just a precaution
    new_eps = ceil(log(abs(G.disc))/(2*n*log(2)) + 3/2)+1;
    new_eps = 2^(-(new_eps +5 ));

    my(neighbours_output, boundary_vector, ctr=1,
        exp_boundary, latticeB, lll_basis_change_matrix, latticeB_lll, scan_elements);
    GP_ASSERT_EQ(length(log_beta), G.r1+G.r2);

    \\#determine the radius we expect to need to search for the desired element
    \\#there 'should' be an element in the search region, assuming the input
    \\#of cpct rep was the log vector of an element.
    boundary_vector = (desired_betalog - log_beta);

    \\#temporarily change precision
    temp_precision = ceil(idealPrecision(G, idealB, normlp(boundary_vector)));
    oldbitprecision = change_precision(temp_precision);

    \\# scale the ideal lattice and convert to real entries
    exp_bvec = exponentiate_logvec(G.r1+G.r2-1, unsquare_log_embeddings(G, boundary_vector), 1);
    scaled_latticeB = mulvec(G[5][1]*idealB, exp_bvec);
    scaled_latticeB = embed_real(G, scaled_latticeB);

    \\# determine the LLL reduced basis and enumerate
    lll_CoB = qflll(scaled_latticeB);   \\# LLL transformation matrix
    lll_ideal = idealB*lll_CoB;         \\# LLL reduced coeff matrix
    scaled_lll_lattice = scaled_latticeB*lll_CoB; \\# LLL-reduced real matrix
    scan_elements = qfminim(scaled_lll_lattice~*scaled_lll_lattice,degree+new_eps*1.0,,2)[3];

    if (DEBUG_CPCT_ENUM,
        print("number of scanned elements: ", length(scan_elements), "  scanned radius = ", degree+new_eps*1.0);
        for(i=1, length(scan_elements),
            check_beta = scaled_lll_lattice*scan_elements[,i];
            print("elt: ", scan_elements[,i], " ",   precision(norml2(check_beta),10));
        );
    );
    \\# loop through returned elements and hope the droid we are looking for is here
    for(i=1, length(scan_elements),
        check_beta = lll_ideal*scan_elements[,i];    \\# beta wrt integral basis
        checkvec = log_beta + get_normalized_log_vector(G, check_beta); \\log(abs( nfeltembed(G, check_beta)));
        \\print("boundary: ", precision(boundary_vector[1..2],10));
        \\print("log beta" , precision(log_beta[1..2],10));

        if(samevecs(desired_betalog, checkvec, new_eps),
            idealB = idealdiv(G, idealB, check_beta);
            change_precision(oldbitprecision);
            \\# the below should never happen
            if(alphaOK != idealB, print("(1) Vector match, but ideal mismatch in cpct enum"); breakpoint(););
            return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
        );

        \\# if the above doesn't work, try the inverse
        checkvec = log_beta - get_normalized_log_vector(G, check_beta); \\log(abs( nfeltembed(G, check_beta)));
        if(samevecs(desired_betalog, checkvec, new_eps),
            idealB = idealmul(G, idealB, check_beta);
            change_precision(oldbitprecision);
            \\# the below should never happen
            if(alphaOK != idealB, print("(2) Vector match, but ideal mismatch in cpct enum");breakpoint());
            return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
        );
    );
    \\# produces a 'bad return value' that should be caught in the calling function
    print("No elements satisfy the condition. Error.");
    return([-1,-1,-1]);
}

cpct_rep_final_enum2(G, idealB, target, desired_betalog, eps, testFlag = 0)={
    if (type(testFlag) == "t_INT" && (testFlag == 1),
        print("final enum");
    );

    SCREEN(0, "Using final enum v2.");
    my(ctr=1,
        unit_rank = G.r1+G.r2-1, field_deg = poldegree(G.pol),
        square_ideal, u_square,
        exp_boundary, latticeB, lll_basis_change_matrix, latticeB_lll, scan_elements);

    square_ideal = idealsquare(idealB, G);
    u_square = pointwise_vector_mul(target,target)~;

    temp_precision = max(50, ceil(idealPrecision(G, square_ideal, normlp(desired_betalog))));
    oldbitprecision = change_precision(temp_precision);
    print(precision(exponentiate_logvec(G.r1+G.r2-1, abs(desired_betalog), 1),10), "  ", precision(u_square, 10));

    \\\ reddiv but look for specific beta
    numerical_ideal = G[5][1]*square_ideal;                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(~numerical_ideal, ~u_square);                             \\ ideal u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = square_ideal*LLL_change_of_basis;                             \\ LLL basis, coords wrt the integral basis
    real_mat_uY = embed_real(~G, ~LLL_numerical_uY);
    scan_elements = qfminim(real_mat_uY~*real_mat_uY, field_deg,,2)[3];
    print("-- elts scanned: ", length(scan_elements));

    for(i=1, length(scan_elements),
        check_beta = LLLcoeffmat*scan_elements[,i];          \\# beta wrt integral basis
        checkvec = get_normalized_log_vector(G, check_beta);
        \\log(abs( nfeltembed(G, check_beta)));
        if(samevecs(desired_betalog, checkvec, eps),
            idealB = idealdiv(G, idealB, check_beta);
            change_precision(oldbitprecision);
            return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
        );
        checkvec = -get_normalized_log_vector(G, check_beta);
        \\ -log(abs( nfeltembed(G, check_beta)));
        if(samevecs(desired_betalog, checkvec, eps),
            idealB = idealmul(G, idealB, check_beta);
            change_precision(oldbitprecision);
            return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
        );
    );

    print("No elements satisfy the condition. Error."); breakpoint();
    return([-1,-1,-1]);
}
/******************************************************************************/
\\ INPUT:
\\ - alpha should be a row vector of the log embedding of a unit
\\ - alphaOK is the Z basis for principal ideal generated by alpha (O_K if alpha a unit)
\\ - G a number field
\\ - eps some error
\\ - avp a prime for which we want all denominators to be coprime to
\\ - see Thiel's description of compact representation. This implementation is applicable as long as
\\   alphaOK is reduced.
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

  print("This is the old compact representation algorithm! No longer verified in the test cases, use at own risk"); breakpoint();
  \\# MAIN LOOP: following the algorithm of Thiel, under the assumption alpha is a unit.
  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  \\# Note: In this case, gamma = 1, and beta1 = 1 also. This means we can skip iteration 1
  \\# If alphaOK is reduced, we should always be able to pick gamma = 1

  [s_term, kprime] = compute_initial_s(alpha, kbound);                                     \\ This is technically the value -s1, since Log(rho) = -alpha
  exp_s = create_target(G, 2*s_term);
  [idealB, target, log_beta, beta] = reddiv_compact(idealB, exp_s, G, G[5][1],avp);        \\ reddiv_compact computes A_2, target, Log(beta_2), beta2

  s_term = -2*s_term;                                                                      \\ this is s2
  log_rho = log_beta;                                                                      \\ rho2 = beta2

  [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 1);
  d_vec = concat(d_vec,ideal_denom);                                                       \\ Get d_i= d(A) and append to tracking vector
  alpha_vec = concat(alpha_vec, alpha_i);                                                  \\ compute alpha_i = d/beta and append to tracking vector

  if(DEBUG_CPCT >0,
    print("log beta and s: ", precision(log_beta,10), "   ", precision(s_term,10));
    print(" - ROUND ", 2, ": beta close enough?: ", samevecs(log_beta[1..unit_rank], s_term, kbound), ": NORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
    );

  \\# TRIVIAL CASE HANDLING
  if (kprime == 1, return([ alpha_vec,d_vec ]));

  for(i=3, kprime,
      [s_term, desired_betalog] = update_tracking_values(s_term, log_rho[1..unit_rank]);  \\s_i and the target beta
      for(ind =1, length(desired_betalog),
        if(desired_betalog[i] < 0, print("negative input"));
      );
      [idealB, target, log_rho, beta] = double_and_reduce(G, idealB, target,log_rho,avp);

      log_beta= log(abs(G[5][1]*beta));

      [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 1);
      d_vec = concat(d_vec,ideal_denom);
      alpha_vec = concat(alpha_vec, alpha_i);

      if(DEBUG_CPCT,
        print("logbeta and s: ", precision(log_beta,10), "   ", precision(s_term,10));
        print(" - ROUND ", i, ": beta close enough?: ", samevecs(log_beta[1..unit_rank], desired_betalog, kbound), ".\nNORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
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

    [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 1);
    d_vec = concat(d_vec,ideal_denom);
    alpha_vec = concat(alpha_vec, alpha_i);

    if(DEBUG_CPCT >0, print(" - ROUND ", kprime+1, ": NORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB)););
    if(DEBUG_CPCT, print(" -- final rho_i ", precision(2*log_rho[1..length(log_beta)] + log_beta,10)););

    return( [alpha_vec, d_vec] );
} \\ end compact_rep_buchmann


\\ # This is compact_rep_buchmann, but it assumes we have r+1 coordinates
\\ # modified so that alphaOK is used as part of the check to ensure the cpct rep
\\ # is corresponding to the same element.
\\ # Note that the minima is always computed assuming it is in Ok. This function is not able to
\\ # is no longer able to use a different ideal as the starting point. See comment [*]
compact_rep_full_input(G, alpha, alphaOK , eps, avp=1, testFlag = 0)={
    \\checkalpha = [-1.94140625, 2.8828125, 0.46093750, -1.99218750, 5.44531250];
    \\if(samevecs(alpha, checkalpha, 1/1000), print("test alpha"); breakpoint());
    \\print("Inside compact rep ");logsum =0;for(i=1, length(alpha),logsum += alpha[i];if (i > G.r1, logsum += alpha[i]););print("log sum : ", precision(logsum,10));breakpoint();
    my(
        unit_rank = G.r1 +G.r2 -1,
        kprime = 1,                                               \\ determines no. steps the algorithm runs for
        kbound,                                                   \\ holds the boundary of the area W
        idealB = matid(poldegree(G.pol)),
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
    GP_ASSERT_EQ(length(alpha), unit_rank+1);

    if(type(alpha) == "t_COL", alpha = alpha~);
    kbound = log(sqrt( abs(G.disc) ))/log(2)/2;                 \\ defines the boundary of the area W

    new_eps = 2^(- (ceil(log(4*poldegree(G.pol))/log(2))+1 + ceil(0.5*log(poldegree(G.pol))/log(2))) );

    \\# MAIN LOOP: following the algorithm of Thiel, under the assumption alpha is a unit.
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\# Note: In this case, gamma = 1, and beta1 = 1 also. This means we can skip iteration 1
    \\# If alphaOK is reduced, we should always be able to pick gamma = 1
    parF = ((2/Pi)^(G.r2))*sqrt(abs(G.disc));               \\\ used to define how much to shrink alpha
    s_term = unsquare_log_embeddings(G, alpha);
    if(normlp(s_term) < ceil(log(parF/poldegree(G.pol))),
        new_t = 0; \\# kprime is by default 1
    ,
        new_k = floor(log( normlp(s_term)*poldegree(G.pol)/ log(parF)) / log(2) )+1;
        new_t = max(0, new_k);
        kprime = new_t+1;
    );
    s_term = 2^(-new_t)*s_term;
    exp_s =  exponentiate_logvec(unit_rank, s_term);

    \\ # [*] observe here that the initial ideal reduction is always done in matid(poldegree(G.pol))
    [idealB, target, log_beta, beta] = reddiv_compact_invert_u(idealB, exp_s, G, G[5][1],avp);
    \\#print("Initial reddiv difference inf norm: ", ceil(normlp(target)));

  s_term = 2*s_term;                                                                       \\ this is s2
  log_rho = log_beta;                                                                      \\ rho2 = beta2
  \\print(G.sign, "   s_term2: ", precision(s_term,10), " log_rho= ", precision(log_rho,10));

  [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 0);
  d_vec = concat(d_vec,ideal_denom);                                                       \\ Get d_i= d(A) and append to tracking vector
  alpha_vec = concat(alpha_vec, alpha_i);                                                  \\ compute alpha_i = d/beta and append to tracking vector

  if ((idealB == alphaOK) && (normlp(log_rho-alpha) < new_eps),
    if (testFlag == 2, GP_ASSERT_TRUE(idealB == alphaOK);); \\#ensure correct ideal returned
    \\#early return statement for cpct rep
    return([alpha_vec, d_vec]);
  );
  if(DEBUG_CPCT > 0,
    print("log beta and s: ", precision(log_beta,10), "   ", precision(s_term,10));
    print(" - ROUND ", 2, ": beta close enough?: ", samevecs(log_beta, s_term, kbound),
        ": NORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
    );

  for(i=3, kprime,
      [s_term, desired_betalog] = update_tracking_values(s_term, log_rho);  \\s_i and the target beta
      [idealB, target, log_rho, beta] = double_and_reduce_invert(G, idealB, target,log_rho,avp);
      \\print(i ,"  iteration reddiv difference inf norm: ", ceil(normlp(target)), "  ", floor(sqrt(abs(G.disc))) );
      \\print("inf-norm target : ", i, "  ", precision(normlp(target),10), "  ", precision((sqrt(abs(G.disc))*(2/Pi)^(G.r2))^(1/poldegree(G.pol)),10) );
      if(abs(ceil(normlp(target))) > floor(sqrt(abs(G.disc))),
          print("WARNING: CPCT REP - target norm is too large");
      );
      log_beta= get_normalized_log_vector(G, beta);
      \\print("log_beta: ", precision(log_beta,10) );
      \\log_beta2= double_complex_coordinates(G.r1, log(abs(G[5][1]*beta)));
      \\GP_ASSERT_TRUE(samevecs(log_beta, log_beta2));

      [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 0);
      d_vec = concat(d_vec,ideal_denom);
      alpha_vec = concat(alpha_vec, alpha_i);

      if(0,
        print(beta, " ", idealB);
        print("- logbeta ", precision(log_beta,10), "\ntar ", precision(target,10));
        print("- log rho ", precision(log_rho,10), "   ", precision(s_term,10));
        print(" - ROUND ", i, ": beta close enough?: ", samevecs(log_beta, desired_betalog, kbound),
            ".\nNORMCHECK: ", intermediate_check(G, alpha_vec, d_vec, idealB) );
      );
  ); \\# end main for loop
    \\print("Check same: ", precision(log_from_cpct(G, [alpha_vec, d_vec]),10), "  " precision(log_rho,10));
    GP_ASSERT_NEAR(norml2(log_from_cpct(G, [alpha_vec, d_vec])- log_rho),0, 2^(-10) );

    \\\# FINAL PART OF THE ALGORITHM
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    s_term = 2*s_term;                                              \\# this is s_{k'+1} = log rho
    desired_betalog = (alpha) - (2*log_rho);                       \\# this equation is "rho - 2*rho_{k'}"
    \\print("desired betalog ", precision(desired_betalog,10));
    [idealB, target, log_rho, beta] = double_and_reduce_invert(G, idealB, target,log_rho);
    log_beta= log(abs(G[5][1]*beta))~;
    log_beta = double_complex_coordinates(G.r1, log_beta);
    \\# use enumeration to find the right value for beta
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\# this is a debug setting for the tests

    if (type(testFlag) == "t_INT" && (testFlag == 1),
        if(!samevecs(abs(desired_betalog), abs(log_beta),new_eps),
            [idealB, log_beta, beta] = cpct_rep_final_enum(G, idealB, beta, log_beta, desired_betalog, alphaOK, new_eps, testFlag);
            if(idealB == -1, ,\\do nothing in this case, this is an error and we need to print the instance
                log_rho += log_beta;
            );
        );
    );
    if(!samevecs(abs(desired_betalog), abs(log_beta),new_eps),
        [idealB, log_beta, beta] = cpct_rep_final_enum(G, idealB, beta, log_beta, desired_betalog, alphaOK, new_eps, testFlag);
        if(idealB == -1, ,\\do nothing in this case, this is an error and we need to print the instance
            log_rho += log_beta;
        );
    );
    \\#clause to catch when cpct_enum fails
    if(idealB == -1,
        print("Final Enum has Failed, printing instance data for debugging: ");
        print(G.pol, "\n ", alphaOK, "  alpha: ", precision(alpha ,50), "  eps: ", new_eps);
        breakpoint();
    );

    \\# finalize the alpha and denom vectors
    [alpha_i, ideal_denom]=get_alpha_and_d(G, idealB, beta, 0);
    d_vec = concat(d_vec,ideal_denom);
    alpha_vec = concat(alpha_vec, alpha_i);

    if(DEBUG_CPCT, print(" -- final rho_i ", precision(2*log_rho[1..length(log_beta)] + log_beta,10)););
    if (testFlag == 2,
        if(idealB != alphaOK, print("Error!: ", G.pol, "\n", idealB," ", alphaOK, "\n", alpha ));
        GP_ASSERT_TRUE(idealB == alphaOK);
    );

    return( [alpha_vec, d_vec] );
} \\ end compact_rep_enum


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
    GP_ASSERT_EQ(matsize(lglat)[1], G.r1+G.r2-1);
    for (i=1, length(lglat),
        extended_log = concat(lglat[,i], extra_log_coordinate(G.r1, G.r2, lglat[,i]));
        \\print("cpct-loglattice: Unit number: ", i, " ", precision(extended_log,25));
        crnew = compact_rep_full_input(G, extended_log, aOK, eps, avp);
        c_rep_list = concat(c_rep_list, [crnew]);
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
        print("denominator ", dlist);
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

\\ G is a number field
\\ r is the unit rank
\\ tracker is a specific vector used in BSGS where
\\ The first r entries are pairs of (compact representation, exponent)
\\ Sub-function of baby-step giant step to compute the logarithm of a list of compact representations/minima
trackerLogarithm(~G, ~tracker, r)={
    logResult = vector(r+1,j,0);
    for(i=1, length(tracker),
        if(i <= r,
            GP_ASSERT_TRUE(length(tracker[i]) == 2);
            logResult += tracker[i][2]*log_from_cpct(G, tracker[i][1]);
        ,\\else
            logResult+= tracker[i][2]*get_normalized_log_vector(G, tracker[i][1]);
        );
    );
    return(logResult);
};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ # FUNCTIONS USED FOR PMAX algorithm (Rigorous index divisor test)
\\ # Utilize a specific means of storing an element as a logarithm matrix along
\\ # with a matrix indicating their exponent
\\ # the exponent matrix is updated at each step by this function
\\ # since it is done per prime, the exponents can be reduced modulo p
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT:
\\ - eta_exp_mat is a integer matrix of exponents
\\ - updatevec indicates how to adjust the new exponent matrix
\\ - k indicates which eta_i need to be updated
\\ - p is the prime index divisor being tested
\\ OUTPUT:
\\ - eta_exp_mat is an updated exponent matrix
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ # THE FUNCTIONS HERE ARE USED FOR THE HEURISTIC INDEX DIVISOR TEST
\\ # REFERRED TO HERE AS pari_prime_check
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ #This function is used in the in the heuristic index divisor test to embed
\\ #a compact representation into a residue field
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

\\ # Embed the element alpha of G represented by a list of cpct reps and exponent
\\ # into a residue field. alpha = prod_{i=1}^{r} (cpct_list[i]^(expvec[i]))
\\ # INPUT:
\\ # - G is a number field
\\ # - cpct_list is a list of r cpct reps
\\ # - expvec is a length r vector
\\ # - res_field is the residue field (a pari object, see function nfmodprinit)
\\ # OUTPUT:
\\ # -
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ # Functions used for testing and debugging
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\#brief - Ensure that ideal = Ok/cpct
\\#param G is a number field
\\#param ideal is an ideal of G
\\#param cpct is a compact representation
\\#detail ASSERT will fail if the ideal Ok/cpct is not equal to input ideal
verify_generator(G, ideal, cpct)=
{
    if(VERIFY_GENERATORS,
        VERIFY_GENERATORS = 0;
        print("WARNING: Verifying ideal generators, expect performance decrease");
    );
    my(
        ndeg = poldegree(G.pol),
        Ok = matid(ndeg),
        tempIdeal = Ok,
        quotientIdeal = idealdiv(G, Ok, compact_reconstruct(G, cpct[1], cpct[2]))
    );

    for(i=1, length(cpct[1]),
        tempIdeal = idealmul(G, tempIdeal, tempIdeal);
        tempIdeal = idealdiv(G,tempIdeal, nfeltdiv(G, cpct[1][i], cpct[2][i]));
    );
    GP_ASSERT_EQ(tempIdeal, quotientIdeal);
    GP_ASSERT_EQ(ideal, quotientIdeal);
}

\\# Debugging function to verify that ideal is generated the product of all
\\# elements in cpctlist
\\# cpct list takes a special form where the firt r elements are cpct reps,
\\# all later elements are just number field elements.
\\ # variable quotientIdeal is computed in the 'obvious way' by reconstructing
\\ # compact reps and dividing the ideal by each of them
\\ # tempIdeal uses the 'better' that does not reconstruct each cpct rep
verify_generator_with_list(G, ideal, cpctlist)=
{
    if(VERIFY_GENERATORS,
        VERIFY_GENERATORS = 0;
        print("WARNING: Verifying ideal generators, expect performance decrease");
    );
    my(
        ndeg = poldegree(G.pol),    r= G.r1+G.r2-1,
        Ok = matid(ndeg),
        tempIdeal = Ok,
        quotientIdeal = Ok;
    );

    for(j=1, length(cpctlist),
        if(j <= r,
            intermediateIdeal = Ok;
            cpct = cpctlist[j][1];
            for(i=1, length(cpct[1]),
                intermediateIdeal = idealmul(G, intermediateIdeal, intermediateIdeal);
                intermediateIdeal = idealdiv(G,intermediateIdeal, nfeltdiv(G, cpct[1][i], cpct[2][i]));
            );
            for(k=1, cpctlist[j][2],
                tempIdeal = idealmul(G, tempIdeal, intermediateIdeal);
                quotientIdeal = idealdiv(G, quotientIdeal, compact_reconstruct(G, cpct[1], cpct[2]));
            );
            GP_ASSERT_EQ(tempIdeal, quotientIdeal);
        ,
            tempIdeal = idealdiv(G, tempIdeal, cpctlist[j][1]);
            quotientIdeal = idealdiv(G, quotientIdeal,cpctlist[j][1] );
        );
    );
    GP_ASSERT_TRUE(check_ideal_reduced(G, ideal));
    GP_ASSERT_EQ(tempIdeal, quotientIdeal);
    GP_ASSERT_EQ(ideal, quotientIdeal);
}
