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

vec_less_than
*/

sqrt2 = sqrt(2);

read("src/MLLL.py");
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

\\ just a method for printing out matrices in a slightly nicer fashion
printmat(M, digits=10) = {
    my(mat_dimensions = matsize(M),
    output_string = "";);
    for(i = 1, mat_dimensions[1],
        output_string="";
        for(j = 1, mat_dimensions[2],
            output_string = concat(concat(output_string, precision(M[i,j],digits)), "     ");
        );
        print(output_string);
    );
};

/********************************************/
/* VECTOR and MATRIX multiplication methods */
/********************************************/

\\ Inverts the coordinates of a column vector y in the ring R^n.
\\ Note that this can take in a row vector, but it makes more sense for columns
invert_coordinates(y)=vector({length(y)},{i},{1/y[i]})~;

/* 2. Pointwise multiply each column vector of mat by vector vec. Essentially the matrix version of pointwise_vector_mul */
/* Outputs a matrix with the same dimensions as mat. The rows of the new matrix are v[i]*row[i] */
mulvec(mat,vec) = matrix({matsize(mat)[1],matsize(mat)[2]},{i},{j},{mat[i,j]*vec[i]});

/*2b. multiply 2 vectors in the ring R^n (pointwise) */
\\ Note that again, columns make more sense than row vectors, but both work
pointwise_vector_mul(u1,u2) = vector(length(u1),{i},{u1[i]*u2[i]})~;

/*3. multiply a matrix y with the inverse of a vector v.
  Note that invert_coordinates(u) inverts entries coordinate wise, then apply mulvec */
dividevector(y,u) = mulvec(y,invert_coordinates(u));

/* Coordinate-wise multiply vector u1 with the inverse of vector u2 in R^n*/
pointwise_vector_div(u1,u2) = pointwise_vector_mul(u1,invert_coordinates(u2));


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
sumvec(v)= {
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
expvec(v) = {
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v
\\ OUTPUT: Vector with entries of [log(abs(v[i]))] as a column
\\ Construct the log vector for x (coordinatewise). If length of x is r, then the logvector has length r-1 */
logvector(v) = vector(length(v), {i}, {log(abs(v[i]))} )~;

/******************************************************************************/
/* Returns an approximation of the vector v */
/******************************************************************************/
vector_approximate(v,eps)= {
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
samevecs(v1,v2, eps)={
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
vec_less_than(v1,v2, axis = 0)={
		if(length(v1)!= length(v2), print("v1, v2 not the same length. Error."); return(-1));
		for(i =1, length(v1),
				if(i != axis,

						if(v1[i]-v2[i] >=0 , return(0););
				);
		);
		return(1);
}

\\ Given two vectors v1 and v2, determines if all entries are within w of each other.
\\ INPUT:
\\ - v1 and v2 are real vectors,
\\ - w is a positive real number.
\\ OUTPUT:
\\ - Return 1 if they are close, 0 otherwise.
check_closeness(v1, v2, w)={
    for(i=1, length(v1),
        if(abs(v1[i]-v2[i]) > w, return(0));
    );
    return(1);
} \\ end check_closeness


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
