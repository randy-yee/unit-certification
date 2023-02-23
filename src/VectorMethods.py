\\ This file contains all the more basic vector functionalities from Ha's code
/*
make_poly
invert_coordinates
mulvec
pointwise_vector_mul
dividevector
pointwise_vector_div
remove_last_coordinate
concat_negative_sum
sumvec
expvec
logvector
vector_approximate
is_trace_zero
check0
samevecs
vec_flip_positive
increment_coordinates
column_lin_comb
change_precision

vec_less_than
*/

sqrt2 = sqrt(2);

read("src/MLLL.py");
read("src/bounds.gp")
read("test/TestingUtility.gp")

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ OUTPUT: a polynomial in x with coefficients determined by L,
\\ L[1] is the coefficient of the largest term, L[last] = constant coefficient
\\ INPUT: in a list of coefficients */
make_poly(coeff_list) = {
  my(output_polynomial = 0);
  for(k =1, length(coeff_list),  output_polynomial += coeff_list[k]*x^(length(coeff_list)-k) );
  return(output_polynomial)
};

/********************************************/
/* VECTOR and MATRIX multiplication methods */
/********************************************/

\\ Inverts the coordinates of a column vector y in the ring R^n.
\\ Note that this can take in a row vector, but it makes more sense for columns
invert_coordinates(~y)=vector({length(y)},{i},{1/y[i]})~;

/* 2. Pointwise multiply each column vector of mat by vector vec. Essentially the matrix version of pointwise_vector_mul */
/* Outputs a matrix with the same dimensions as mat. The rows of the new matrix are v[i]*row[i] */
mulvec(~mat,~vec) = matrix({matsize(mat)[1],matsize(mat)[2]},{i},{j},{mat[i,j]*vec[i]});

/*2b. multiply 2 vectors in the ring R^n (pointwise) */
\\ Note that again, columns make more sense than row vectors, but both work
pointwise_vector_mul(~u1, ~u2) = vector(length(u1),{i},{u1[i]*u2[i]})~;

/*3. multiply a matrix y with the inverse of a vector v.
  Note that invert_coordinates(u) inverts entries coordinate wise, then apply mulvec */
dividevector(y,u) = mulvec(y,invert_coordinates(u));

/* Coordinate-wise multiply vector u1 with the inverse of vector u2 in R^n*/
pointwise_vector_div(~u1,~u2) = pointwise_vector_mul(u1,invert_coordinates(u2));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v
\\ OUTPUT: v with an extra term equal to the negative sum of the entries of v
concat_negative_sum(v) = concat(v,-sumvec(v));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v
\\ OUTPUT: v with last coordinate removed
remove_last_coordinate(v) = {v[1..length(v)-1]};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v
\\ OUTPUT: The sum of the elements of the vector v in R^n
sumvec(~v)= {
    my(entry_sum:small = 0);
    for (i = 1,length(v),
        entry_sum += v[i]
    );
    return(entry_sum);
}
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v
\\ OUTPUT: A vector [exp(-v[1]),exp(-v[2]), ... exp(sum)]
\\ Note that the output has length one more than v
expvec(~v) = {
    my(xn);
    xn = concat(v,-sumvec(v));                \\ compute the negative sum of x, concatenate to xn
    vector(length(xn), i, exp(-xn[i]) );    \\ create a vector [exp(-x[1]),exp(-x[2]), ... exp(sum)]
}
inverse_expvec(v) = {
    my(xn);
    xn = concat(v,-sumvec(v));                  \\ compute the negative sum of x, concatenate to xn
    xn = vector(length(xn), i, exp(xn[i]) );    \\ create a vector [exp(x[1]),exp(x[2]), ... exp(sum)]
    return(xn);
}

/******************************************************************************/
/* Returns an approximation of the vector v */
/******************************************************************************/
vector_approximate(~v,eps)= {
    my( vn1=eps*round(v*1/(eps))*1.00 );
    vn1;
}

/******************************************************************************/
/*12. test that v has sum of the coordinate = 0 --> precision for comparision to 0 - similar to the comparision to 1 in the function check0 */
/******************************************************************************/
is_trace_zero(v,eps)={
    if(sumvec(v)<eps^2,
        1,
    0);       \\ else
}:bool

/******************************************************************************/
/*9. check that the coordinates of vector v in R^n are all abs. val. less than 1 or not (PRECISION for comparisons all coordinates of v with 1)*/
check0(v,eps) = {
    for (ctr = 1, length(v),
        if ( 1-abs(v[ctr]) < eps,
            return(0)
        );
    );
    return(1);
}:bool;

\\ Check whether two vectors are the same, up to some error
\\ INPUT:
\\ - v1 and v2 are both vectors, eps is the acceptable error
\\ OUTPUT:
\\ - 0 if the entries of v1 and v2 are not each within eps of each other, 1 otherwise.
samevecs(~v1,~v2, eps)={
  if(length(v1) != length(v2), return(0));

  for (i = 1, length(v1),
    if( abs(v1[i] - v2[i]) > eps, return(0));
  );
  return(1);
};

\\ Method for ensuring the that first nonzero entry of a vector is positive

vec_flip_positive(v1)={
  for(i=1, length(v1),
      if(v1[i] != 0,
          if(v1[i] < 0, return(-v1), return(v1));
      );
  );
  return(v1);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Checks whether each of the entries of v1 are smaller than the corresponding entry of v2
\\ INPUT:
\\ - v1, v2 are vectors
\\ - axis is an integer. Ignore comparison of that entries for that index.
\\   For application in the expand_minset algorithm. Leave as 0 otherwise
\\ OUTPUT
\\ - 0 if v1 equals or exceeds v2 on any coordinate (excluding axis)
\\ - 1 if v1 is smaller than v2 on each coordinate (except possibly axis)
vec_less_than(~v1,~v2, axis = 0)={
		if(length(v1)!= length(v2), print("v1, v2 not the same length. Error."); return(-1));
		for(i =1, length(v1),
				if(i != axis,
						if(v1[i]-v2[i] >=0 , return(0););
				);
		);
		return(1);
}

\\ Computes the infinity norm of a vector v
\\ INPUT:
\\ - A real vector v
\\ OUTPUT:
\\ - Real number corresponding to the infinity norm of v
infinity_norm(v)={
  my(len=length(v), abv, val = 0);
  for (i=1, len,
    abv = abs(v[i]);
    if (abv > val, val = abv);
  );
  return(val);
}


\\# INPUT:
\\ - A vector of integers A which defines the subdivision for the babystock
\\ - A vector of integers V that defines a particular 'giant step'
\\   for all i, V[i] < A[i]
\\# OUTPUT:
\\ - The vector V 'incremented' by 1, similar to a number in which each digit
\\  is a different base
increment_coordinates(a_vec, ~current_vec)={
    my(place = length(current_vec), carryflag = 0;);
    current_vec[place]+=1;
    if(current_vec[place] >= a_vec[place], carryflag =1);
    while(carryflag == 1,
        current_vec[place] = 0; carryflag = 0;
        if(place == 1, return(current_vec));
        place -= 1;
        current_vec[place]+=1;
        if(current_vec[place] >= a_vec[place], carryflag =1);
    );
    return;
}

column_lin_comb(~lattice, ~coeff_vector)=
{
    GP_ASSERT_EQ(length(coeff_vector) , length(lattice));
    lc_vector = vector(matsize(lattice)[1], i, 0)~;
    GP_ASSERT_TRUE(type(lc_vector) == "t_COL");
    for(i=1, length(lattice),
        lc_vector += (coeff_vector[i]*lattice[,i]);
    );
    return(lc_vector);
}

change_precision(newbitprec)=
{
    oldbitprecision = default(realbitprecision);
    default(realbitprecision, newbitprec);
    return(oldbitprecision);
}


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


gram_schmidt(~real_lattice)=
{
    my(rank= length(real_lattice), ortho_basis=real_lattice);
    for(i =2, rank,
        for(j=1, i-1,
            mu_ij = (real_lattice[,i]~ * ortho_basis[,j])/norml2(ortho_basis[,j]);
            ortho_basis[,i] -= mu_ij*ortho_basis[,j];
        );
    );
    return(ortho_basis);
}

get_enumeration_bounds(degree, lll_lattice)=
{
    my(rank = length(lll_lattice),
        ortho_basis, norm_vector, k_vector
    );

    if(norml2(lll_lattice[,1])<1, return 0);
    ortho_basis = gram_schmidt(lll_lattice);
    norm_vector = vector(rank, i, norml2(ortho_basis[,i]));
    k_vector = vector(rank, i, (3/sqrt(2))^(degree - i)*sqrt(degree)/norm_vector[i] );
    k_vector = floor(k_vector);
    return(k_vector);
}


check_ideal_reduced(G, ideal)=
{
    if(abs(1/matdet(idealhnf(G, ideal)))>sqrt(abs(G.disc)), return(0));

    my(rank = length(ideal),
        k_vector, zero_vec, iteration_vector);
    ideal_lattice = G[5][1]*ideal;
    ideal_lattice = embed_real(G, ideal_lattice);

    lll_basis_change_matrix = qflll(ideal_lattice);
    lll_lat = ideal_lattice*lll_basis_change_matrix;    \\#real lattice
    lll_ideal = ideal*lll_basis_change_matrix;          \\#ideal representation
    ortho_basis = gram_schmidt(lll_lat);                \\#orthogonalized
    k_vector = get_enumeration_bounds(rank, lll_lat);  \\# compute the k vector
    k_vector = vector(rank, i, k_vector[i]+1);
    zero_vec = vector(rank, i , 0);
    iteration_vector = zero_vec;
    increment_coordinates(k_vector, ~iteration_vector);
    one_vec = vector(G.r1+G.r2, i , 1);
    temp_bit_precision = max(10, ceil(log(denominator(ideal))+4+(rank^2*log(4*denominator(ideal)^2))+2));
    mainbitprecision = default(realbitprecision);
    default(realbitprecision, temp_bit_precision);  \\#save and change precision
    while(iteration_vector != zero_vec,
        test_vector = column_lin_comb(~lll_ideal, ~iteration_vector);
        if(vec_less_than(abs(nfeltembed(G, test_vector)), one_vec),
            default(realbitprecision, mainbitprecision);    \\#restore precision
            return(0)
        );
        increment_coordinates(k_vector, ~iteration_vector);
    );
    default(realbitprecision, mainbitprecision);    \\#restore precision
    return(1);  \\# no minima found inside of the normed body of 1
}
