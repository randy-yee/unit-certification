\\ A loose collection of functions operating on vectors
/*
invert_coordinates
mulvec
pointwise_vector_mul
pointwise_vector_div
concat_negative_sum
sumvec
logvector
vector_approximate
is_trace_zero
check0
samevecs
vec_flip_positive
vec_less_than
double_complex_coordinates
infinity_norm
increment_coordinates
column_lin_comb
change_precision
gram_schmidt
get_enumeration_bounds
compute_subgroup
get_scaled_ideal_lattice
check_ideal_reduced
*/

read("src/MLLL.py");
read("src/bounds.gp")
read("test/TestingUtility.gp")

/********************************************/
/* VECTOR and MATRIX multiplication methods */
/********************************************/

\\ Inverts the coordinates of a real vector y
\\ Note that this can take in a t_VEC, but always returns a t_COL
invert_coordinates(~y)=vector({length(y)},{i},{1/y[i]})~;

/* 
 *  BRIEF: Pointwise multiply each row of mat by vector vec.
 *  Corresponds to operation V * M, where V is the diagonal matrix whose entries are given by vec
 *  Outputs a matrix with the same dimensions as mat whose rows are v[i]*row[i] 
*/
mulvec(~mat,~vec) = matrix({matsize(mat)[1],matsize(mat)[2]},{i},{j},{mat[i,j]*vec[i]});


\\ BRIEF: Pointwise multiplication of two vectors/columns 
\\ OUTPUT:
\\ - coordinate-wise product of the two inputs. Always returns a column
pointwise_vector_mul(~u1, ~u2) = vector(length(u1),{i},{u1[i]*u2[i]})~;

/*3. multiply a matrix y with the inverse of a vector v.
  Note that invert_coordinates(u) inverts entries coordinate wise, then apply mulvec */
dividevector(~y,~u) = mulvec(y,invert_coordinates(u));

/* Coordinate-wise multiply vector u1 with the inverse of vector u2 in R^n*/
pointwise_vector_div(~u1,~u2) = pointwise_vector_mul(u1,invert_coordinates(u2));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v of nonzero length
\\ OUTPUT: The sum of the coordinates of the vector v
sumvec(~v)= {
    my(entry_sum:small = 0);
    for (i = 1,length(v),
        entry_sum += v[i]
    );
    return(entry_sum);
}
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ BRIEF: Return (v || c) where c is equal to the negative sum of entries of v
\\ INPUT:
\\ - A vector v of length t
\\ OUTPUT:
\\ - A vector (v || c) of length (t+1)
concat_negative_sum(v) = concat(v,-sumvec(v));


/******************************************************************************/
/* Returns an approximation of the vector v */
/* Not terribly efficient. maybe this can be done using strings? */
/******************************************************************************/
vector_approximate(~v,eps)= {
    my( vn1=eps*round(v*1/(eps))*1.00 );
    vn1;
}
vector_approximate2(~v,eps)= {
    my( vn1=eps*round(v*1/(eps))*1.00 );
    vn1;
}


/******************************************************************************/
/*12. sum coordinates, and see if close enough to 0 */
/******************************************************************************/
is_trace_zero(v,eps)={
    GP_ASSERT_TRUE(eps > 0);
    if(abs(sumvec(v))<eps, return(1), return(0));
}:bool


\\ BRIEF: Checks if the real vector v is in the normed body of 1.
\\ i.e. determines whether the magnitude of the coords of v are less than 1
\\ NOTE: Requires PRECISION -- eps is used in comparing coordinates of v with 1)
\\ DEPRECATE: Appears in 
\\ scanball_map, one_massive_qfminim, overlap_scanball, compact_storage_overlap_scanball
check0(~v,eps) = {
    GP_ASSERT_TRUE(eps > 0);
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
samevecs(~v1,~v2, ~eps)={
  GP_ASSERT_TRUE(eps > 0);
  if(length(v1) != length(v2), return(0));
  for (i = 1, length(v1),
    if( abs(v1[i] - v2[i]) > eps, return(0));
  );
  return(1);
};


\\ BRIEF: Return adjusted vector so that the first nonzero entry is positive
\\ INPUT:
\\ - a vector v
\\ OUTPUT:
\\ - either v or -v, based on the sign of the first nonzero coefficient
vec_flip_positive(v1)={
  for(i=1, length(v1),
      if(v1[i] != 0,
          if(v1[i] < 0, return(-v1), return(v1));
      );
  );
  return(v1);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ BRIEF: Checks if the entries of v1 are smaller than that of of v2
\\ - If the value of axis is set, skip comparison on the specified coordinate.
\\ INPUT:
\\ - VEC_t : v1, v2
\\ - t_INT axis
\\ OUTPUT
\\ - 0 if v1 equals or exceeds v2 on any coordinate (except possibly axis)
\\ - 1 if v1 is smaller than v2 on each coordinate (except possibly axis)
vec_less_than(~v1,~v2, axis = 0)={
    GP_ASSERT_EQ(length(v1), length(v2));
	for(i =1, length(v1),
		if(i != axis,
			if(v1[i]-v2[i] >=0 , return(0););
		);
	);
	return(1);
}

\\ BRIEF:
\\ - Assuming the input is a log vector for some number field element
\\ - doubles the entries corresponding to the complex embeddings.
\\ INPUT:
\\ - r1 is the number of real coordinates
\\ - v is the vector
\\ OUTPUT:
\\ - double the coordinates of v that do not correspond to complex embeddings
double_complex_coordinates(~r1, ~v)=
{
    for(i=r1+1, length(v),
        v[i] *=2
    );
    return(v);
}

\\ wrapper function for normlp(v), because I forgot that function existed..
\\ INPUT:
\\ - A real vector v
\\ OUTPUT:
\\ - Real number corresponding to the infinity norm of v
infinity_norm(~v)={
  return(normlp(v));
}

\\ BRIEF:
\\ - Simulates addition by 1, where current_vec is a vector of digits
\\   and a_vec are the maximums at each position.
\\   ex1) Let a_vec = [4,5,6]. Then increment_coordinates([2,3,4]) = [2,3,5]
\\   ex2) increment_coordinates(a_vec, [2,3,5]) = [2,4,0]
\\   ex3) increment_coordinates(a_vec, [3,4,5]) = [0,0,0]
\\ INPUT:
\\ - a_vec: A vector of integers denoted by A
\\ - current_vec: A vector of integers V such that  for all i, V[i] < A[i]
\\ OUTPUT:
\\ - MODIFY IN PLACE: V -> V' incremented in the way described above
\\ - Note that when called, the argument current_vec MUST be provided with a '~' character in front
increment_coordinates(a_vec, ~current_vec)={
    \\ It is assumed that a_vec[i] > current_vec[i] for all i, but not checked!
    my(place = length(current_vec), carryflag = 0;);
    current_vec[place]+=1;
    if(current_vec[place] >= a_vec[place], carryflag =1);
    while(carryflag == 1,
        current_vec[place] = 0; carryflag = 0; \\# implies rollover; set value to 0
        if(place == 1, return);                \\# if at leftmost place value, return
        place -= 1;                            \\# move place-value left
        current_vec[place]+=1;                 \\# apply carryover +1
        if(current_vec[place] >= a_vec[place], carryflag =1); \\#determine if this value rolls over
    );
    return;
}

\\ BRIEF:
\\ - compute a linear combination of the columns of a matrix. 
\\ INPUT:
\\ - A matrix called lattice of dimension (t,s)
\\ - A coeff_vector of length s
\\ OUTPUT:
\\ - a column vector of size t
column_lin_comb(~lattice, ~coeff_vector)=
{
    GP_ASSERT_EQ(length(coeff_vector) , length(lattice));
    my(lc_vector = vector(matsize(lattice)[1], i, 0)~);
    GP_ASSERT_TRUE(type(lc_vector) == "t_COL");
    for(i=1, length(lattice),
        lc_vector += (coeff_vector[i]*lattice[,i]);
    );
    return(lc_vector);
}

\\ BRIEF:
\\ - helper function to change the working precision of GP
\\ INPUT:
\\ - an integer value newbitprec
\\ OUTPUT:
\\ - an integer value equal to the old bit precision
change_precision(newbitprec)=
{
    oldbitprecision = default(realbitprecision);
    default(realbitprecision, newbitprec);
    return(oldbitprecision);
}

\\ BRIEF:
\\ - perform gram schmidt on the input matrix over the reals
\\ INPUT:
\\ - real_lattice is a matrix
\\ OUTPUT:
\\ - a matrix whose columns should be orthogonal
gram_schmidt(~real_lattice)=
{
    my(rank= length(real_lattice), ortho_basis=real_lattice);
    for(i = 2, rank,
        for(j=1, i-1,
            mu_ij = (real_lattice[,i]~ * ortho_basis[,j])/norml2(ortho_basis[,j]);
            ortho_basis[,i] -= mu_ij*ortho_basis[,j];
        );
    );
    return(ortho_basis);
}

\\ BRIEF:
\\ - determine a vector indicating the maximum possible coefficients of a
\\ - a shortest vector; used to determine which vectors need to be checked (DEBUGSCALING- norml2)
\\ INPUT:
\\ - degree of underlying number field,
\\ - an matrix representation of a lattice (LLL-reduced)
\\ OUTPUT:
\\ - a vector whose size matches the number of columns of the input matrix
get_enumeration_bounds(degree, ~lll_lattice)=
{
    my(rank = length(lll_lattice),
        ortho_basis, norm_vector, k_vector
    );

    if(norml2(lll_lattice[,1])<1, return 0);
    ortho_basis = gram_schmidt(lll_lattice);
    norm_vector = vector(rank, i, sqrt(norml2(ortho_basis[,i])));
    k_vector = vector(rank, i, (3/sqrt(2))^(degree - i)*sqrt(degree)/norm_vector[i] );
    k_vector = floor(k_vector);
    return(k_vector);
}

\\ BRIEF:
\\\ A function used for testing. Given a matrix corresponding to Lambda_K
\\\ construct a new matrix corresponding to a sublattice of Lambda_K
\\\ Helpful for testing the ability to find index divisors
\\ INPUT:
\\ - G is a number field
\\ - unimat is the matrix corresponding to log lattice of G
\\ - modpair/extype are used to modify the lattice (see code)
\\ OUTPUT:
\\ - a matrix corresponding to a sublattice of the input lattice
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\# IDEAL METHODS
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ - K is a number field
\\ - ideal_I is a coefficient matrix for an ideal of K
\\ - scale_vector is a length r+s+1 vector.
\\ OUTPUT - the ideal lattice corresponding to ideal_I; with ith row scaled
\\          by the ith entry of scale_vector, converted to a real matrix
get_scaled_ideal_lattice(K, ~ideal_I, scale_vector=0)=
{
    my(result);
    result = K[5][1]*ideal_I;
    if(scale_vector == 0, return(embed_real(K, result)););
    result = mulvec(result, scale_vector);
    return(embed_real(K, result));
}


\\ BRIEF:
\\ - Check if the input ideal is a reduced ideal of number field G
\\ INPUT:
\\ - G a number field
\\ - ideal is a coefficient matrix of an ideal wrt integral basis of G
\\ OUTPUT:
\\ - boolean indicating whether the ideal is reduced or not
check_ideal_reduced(~G, ~ideal)=
{
    \\# if the norm is too small, then its not reduced
    if(abs(1/idealnorm(G, ideal))>sqrt(abs(G.disc)), return(0));
    my(rank = length(ideal),
        k_vector, zero_vec, iteration_vector,
        complex_basis, real_basis, 
        basis_change,
        real_basis_lll,complex_lattice_lll,
        check_elements, one_vec, test_vector_real
        );

    complex_basis = G[5][1]*ideal;
    real_basis = embed_real(G, complex_basis);

    basis_change = qflll(real_basis);
    real_basis_lll = real_basis*basis_change;    \\#real lattice (LLL)

    \\# first element is in normed body of 1, so not reduced
    if(norml2(real_basis_lll[,1]) < 1, return(0));

    \\#complex lattice (LLL)
    complex_lattice_lll = complex_basis*basis_change;

    \\# scan radius with a small error factor, capture potential elements
    \\# in the the normed body of 1
    check_elements = qfminim(real_basis_lll~*real_basis_lll,poldegree(G.pol)+0.065,,2);
    \\#if no elements, then the normed body only contains the 0 vector
    if(check_elements[1] == 0, return(1));
    one_vec = vector(G.r1+G.r2, i, 1);
    for(i=1, length(check_elements[3]),
        \\#check if any elements are in the normed body of 1
        \\(DEBUGSCALING - verify abs vs norm)
        test_vector_real = abs(complex_lattice_lll*check_elements[3][,i]);
        if(vec_less_than(test_vector_real, one_vec),
            return(0)
        );
    );
    return(1);

}
