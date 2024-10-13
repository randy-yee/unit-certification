\\ This file contains all the more basic vector functionalities from Ha's code
/*
make_poly
invert_coordinates
mulvec
pointwise_vector_mul
dividevector
pointwise_vector_div
concat_negative_sum
sumvec
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
\\ INPUT:
\\ - in a list L containing t coefficients ordered largest exponent to smallest
\\ OUTPUT:
\\ - a polynomial f in x of degree (t-1) with coefficients determined by L

make_poly(coeffs) = {
  my(poly_f = 0);
  for(k =1, length(coeffs),  poly_f += coeffs[k]*x^(length(coeffs)-k) );
  return(poly_f)
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
dividevector(~y,~u) = mulvec(y,invert_coordinates(u));

/* Coordinate-wise multiply vector u1 with the inverse of vector u2 in R^n*/
pointwise_vector_div(~u1,~u2) = pointwise_vector_mul(u1,invert_coordinates(u2));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ BRIEF: Return (v || c) where c is equal to the negative sum of entries of v
\\ INPUT:
\\ - A vector v of length t
\\ OUTPUT:
\\ - A vector (v || c) of length (t+1)
concat_negative_sum(v) = concat(v,-sumvec(v));

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v of length t
\\ OUTPUT: The sum of the coordinates of the vector v
sumvec(~v)= {
    my(entry_sum:small = 0);
    for (i = 1,length(v),
        entry_sum += v[i]
    );
    return(entry_sum);
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
check0(~v,eps) = {
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
  if(length(v1) != length(v2), return(0));

  for (i = 1, length(v1),
    if( abs(v1[i] - v2[i]) > eps, return(0));
  );
  return(1);
};

\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
\\ BRIEF: Checks if each of the entries of v1 are smaller than that of of v2
\\ - axis allows us to skip comparison on the specified coordinate, which is
\\   needed in the expand_minset function. Leave as 0 otherwise
\\ INPUT:
\\ - vectors v1, v2
\\ - integer value axis
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

double_complex_coordinates(~r1, ~v)=
{
    for(i=r1+1, length(v),
        v[i] *=2
    );
    return(v);
}

\\ Computes the infinity norm of a vector v
\\ INPUT:
\\ - A real vector v
\\ OUTPUT:
\\ - Real number corresponding to the infinity norm of v
infinity_norm(~v)={
  my(len=length(v), abv, val = 0);
  for (i=1, len,
    abv = abs(v[i]);
    if (abv > val, val = abv);
  );
  return(val);
}

\\ BRIEF:
\\ - Simulates addition by 1, where current_vec is a vector of digits
\\   and a_vec is the roll-over value.
\\   ex) a_vec = [4,5,6], Then increment_coordinates(a_vec, [2,3,5]) = [2,4,0]
\\   ex) increment_coordinates(a_vec, [3,4,5]) = [0,0,0]
\\ INPUT:
\\ - A vector of integers A
\\ - A vector of integers V such that  for all i, V[i] < A[i]
\\ OUTPUT:
\\ - MODIFY IN PLACE: V -> V' incremented in the way described above
increment_coordinates(a_vec, ~current_vec)={
    \\# Note that it is assumed that a_vec[i] > current_vec[i] for all i
    \\# but this is not checked!
    my(place = length(current_vec), carryflag = 0;);
    current_vec[place]+=1;
    if(current_vec[place] >= a_vec[place], carryflag =1);
    while(carryflag == 1,
        current_vec[place] = 0; carryflag = 0; \\# implies rollover; set value to 0
        if(place == 1, return);                \\# if this was the leftmost place value, return
        place -= 1;                            \\# move one place-value left
        current_vec[place]+=1;                 \\# apply the carryover +1
        if(current_vec[place] >= a_vec[place], carryflag =1); \\#determine if this value rolls over and repeat
    );
    return;
}

\\ BRIEF:
\\ - compute a linear combination of the columns of a matrix
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
\\ - changes the working precision of GP
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
\\ - perform gram schmidt on the input matrix over the reals (DEBUGSCALING- norml2)
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
\\          by the ith entry of scale_vector
get_scaled_ideal_lattice(K, ideal_I, scale_vector)=
{
    my(result);
    result = K[5][1]*ideal_I;
    result = mulvec(result, scale_vector);
    return(embed_real(K, result));

}


\\ BRIEF:
\\ -
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
        k_vector, zero_vec, iteration_vector);
    ideal_lattice = G[5][1]*ideal;
    ideal_lattice = embed_real(G, ideal_lattice);

    lll_basis_change_matrix = qflll(ideal_lattice);
    lll_lat = ideal_lattice*lll_basis_change_matrix;    \\#real lattice
    lll_ideal = ideal*lll_basis_change_matrix;          \\#ideal representation
    if(norml2(lll_lat[,1]) < 1, return(0));
    ortho_basis = gram_schmidt(lll_lat);                \\#orthogonalized
    k_vector = get_enumeration_bounds(rank, lll_lat);  \\# compute the k vector

    check_elements = qfminim(lll_lat~*lll_lat,poldegree(G.pol)+0.1,,2);
    \\#if no elements, then the normed body only contains the 0 vector
    if(check_elements[1] == 0, return(1));
    one_vec = vector(G.r1+G.r2, i , 1);
    for(i=1, length(check_elements[3]),
        \\#if there were elements in the scan region, check if they are in the normed body of 1
        \\(DEBUGSCALING- abs vs norm)
        test_vector_real = abs(G[5][1]*lll_ideal*check_elements[3][,i]);
        if(vec_less_than(test_vector_real, one_vec),
            return(0)
        );
    );
    return(1);

}

get_scaled_M(K)=
{
    my(
        scaled_M = K[5][1],
        num_embeddings = K.r1 + K.r2,
        fld_degree = length(K.zk)
    );
    for(i = K.r1+1, num_embeddings,
        scaled_M[i,] = 2*scaled_M[i,];
    );
    return scaled_M;
}
/*
\\ BRIEF:
\\ - DEPRECATED
\\ INPUT:
\\ - G a number field
\\ - ideal is a coefficient matrix of an ideal wrt integral basis of G
\\ OUTPUT:
\\ - boolean indicating whether the ideal is reduced or not
check_ideal_reduced_old(G, ideal)=
{
    \\ if num elts is 0, reduced, else check each of them if they are in the normed body of 1
    k_vector = vector(rank, i, k_vector[i]+1);
    print(k_vector);
    zero_vec = vector(rank, i , 0);
    iteration_vector = zero_vec;
    increment_coordinates(k_vector, ~iteration_vector);
    one_vec = vector(G.r1+G.r2, i , 1);
    temp_bit_precision = max(10, ceil(log(denominator(ideal))+4+(rank^2*log(4*denominator(ideal)^2))+2));
    mainbitprecision = default(realbitprecision);
    default(realbitprecision, temp_bit_precision);  \\#save and change precision

    complex_ideal_lattice = G[5][1]*lll_ideal;
    while(iteration_vector != zero_vec,
        \\test_vector = column_lin_comb(~lll_ideal, ~iteration_vector);
        test_vector_real = abs(complex_ideal_lattice*iteration_vector~);
        \\print(precision(abs(nfeltembed(G, test_vector)),10), "  \n", precision(test_vector_real,10), "\n");
        if(vec_less_than(test_vector_real, one_vec),
            default(realbitprecision, mainbitprecision);    \\#restore precision
            return(0)
        );
        increment_coordinates(k_vector, ~iteration_vector);
    );
    default(realbitprecision, mainbitprecision);    \\#restore precision
    return(1);  \\# no minima found inside of the normed body of 1
}


*/
