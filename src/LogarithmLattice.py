/*

testequalOK
minkowski_absval
valuationvec
absoluteval_nvec
unsquare_log_embeddings
units_to_matrix
log_determinant


embed_real
undo_real
get_real_vec
checkred
get_ideal_denom
ideal_contains1
limitminvector

is_vec_in_lattice
get_real_vec
process_complex_loglattice

get_log_lattice_bnf
get_unscaled_determinant

is_in_axis_box
verify_lattice_containment
*/


/******************************************************************************/
/* Tests if a frac. ideal y is equal to OK by computing the HNF and comparing the identity matrix. */
/* @params y is a matrix representing a frac ideal
\\           G is the number field output of nfinit
*/
/******************************************************************************/
testequalOK(y,K)={
    my(n=length(y));
    my(yhnf);

    yhnf=idealhnf(K,y);
    if(yhnf==matid(n),
        1,
    0);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A vector v, r1 the number of real embeddings
\\ OUTPUT: Takes the abs on all entries 1...r1, norm on entries r1+1 ... r1+r2
\\ The point of the is function is that if you use nfeltembed on a number field
\\ element the complex embeddings are not squared, which means using log(abs())
\\ will not use the regulator. Useful for debugging
minkowski_absval(v, r1)={
    for(i=1, length(v),
        if(i <=r1, v[i] = abs(v[i]),
        v[i] = norm(v[i]));
    );
    return(v);
}


/******************************************************************************/
\\ INPUT: G is a number field, elt is a number field elt that is either a polynomial or a column
\\ Should set the variable column = 1 if the input is in column form, but
\\ Pari's nfalgtobasis function knows to not do anything if you forget to do this
\\ OUTPUT: the dimension r1+r2 vector (|x_1|, ... |x_r1|, |x_{r1+1}|^2, ... |x_r2|^2)
/* IMPORTANT to note the squares on the complex embeddings.

/******************************************************************************/
valuationvec(~G,~elt, column = 0)={
	my(embedvec, column_elt);
    \\ process input so that it is a column vector
    if((column == 0) || (column ==1 && type(elt)== "t_INT"), column_elt = nfalgtobasis(G,elt);, column_elt = elt );
    embedvec = (G[5][1]*column_elt)~;  \\ convert to vector
    embedvec = minkowski_absval(embedvec, G.r1);
	return(embedvec);
}


\\ INPUT:
\\ - G a number field
\\ - v a coefficient vector for some element in G
\\ - the default param side=2 is used when you want complex embeddings next to each other.
\\ Change this to another value if you want it so that the order has distinct embeddings first,
\\ then the conjugates in the same order. ie. [s1, s2, s_{r1+1}, ... s_{r1+r2}, conj(s_{r1+1}), conj(s_{r1+2}), ... conj(s_{r1+r2})]
\\ OUTPUT:
\\ A vector of length n consisting of the absolute value of each embedding.
absoluteval_nvec(~G, ~v, side = 2)={
    my(embed_v, outvec);
    embed_v = abs(G[5][1]*v)~;                                                  \\ vector of length r1+r2 corresponding to valuations

    if(side == 2,
        \\ complex conjugates will be side by side
        outvec = embed_v[1..G.r1];
        for(i = G.r1+1, G.r1 + G.r2, outvec = concat(outvec, [embed_v[i], embed_v[i]]));
        return(outvec);
    ,\\else
        embed_v = concat(embed_v, embed_v[(G.r1+1) .. (G.r1+G.r2)]);
        return(embed_v);
    );
}


\\ When you're dealing with log vectors of the form (log|x_1|, ... log|x_r1|, log|x_{r1+1}|^2, ... log|x_(r2-1)|^2)
\\ this eliminates the squares by dividing by 2.
\\ INPUT:
\\ - G a number field
\\ - logvec, the log embedding of an element, of the form described above.
\\ OUTPUT:
\\ - a log vector of the form  (log|x_1|, ... log|x_r1|, log|x_{r1+1}|, ... log|x_(r2-1)|). Note the removal of the squares on the complex entries.
unsquare_log_embeddings(G, logvec)={
    my(newlog);
    newlog = vector(length(logvec), i, if(i <= G.r1, logvec[i], logvec[i]/2));
    return(newlog);
}


\\ Reads in the fundamental units of a number field (bnf.fu)
\\ returns those units as a coefficient vector in terms of the integral basis.
units_to_matrix(nf, f_units)={
    my(unitmat = nfalgtobasis(nf, f_units[1]));
    if(length(f_units) > 1,
        for(i =2, length(f_units),
            unitmat = matconcat( [unitmat, nfalgtobasis(nf, f_units[i])] );
        );
    );
    unitmat = Mat(unitmat);
    return(unitmat);
};

\\ If the field regulator is known, can be used to see if a sublattice of the unit group is actually the unit group itself, or determine the index
\\ INPUT:
\\ - A Number field G
\\ - a matrix representing a lattice in G
\\ OUTPUT:
\\ - Gives the (equivalent of the) regulator of the matrix.

log_determinant(G, unitmat)={
    my(unit_embedding, square_matrix);
    unit_embedding = G[5][1]*unitmat;
    square_matrix = matrix( length(unit_embedding), length(unit_embedding), i,j, unit_embedding[i,j] ) ;
    square_matrix = log(abs(square_matrix));
    for(i = 1, length(square_matrix), if(i > G.r1, square_matrix[i,]*=2; ));

    return(abs(matdet(square_matrix)));
}



\\ Used to power an element in a numberfield where the coefficients are reduced mod a prime p (i.e. in the ring Zk / pZk )
\\ INPUT:
\\ - A Number field G
\\ - a prime p
\\ - an element elt
\\ - an exponent expo
\\ OUTPUT:
\\ - Gives the (equivalent of the) regulator of the matrix.

nfeltpow_modp(G, p, elt, expo )={
    my(bin, base, result);
    bin = binary(expo);
    base = elt%p;
    result = base;

    for(i=2, length(bin),
        if(bin[i] == 1,
            result = nfeltmul( G, nfeltpow(G,result, 2)%p, base ) % p;,
        \\else
            result = nfeltpow(G,result, 2)%p;
        );
    );
    return(result);
}

\\ converts the embedding vector of v (which could be complex), into a real vector
\\ INPUT:
\\ - G a number field
\\ - v a dimension G.r1+G.r2 vector
\\ OUTPUT:
\\ - the real vector of v having dimension R^n, n is the degree of the number field
get_real_vec(~G, ~v)={
  my(vec1, realvec=[]);
  sqrt2 = sqrt(2);
  vec1 = v;
  realvec = vec1[1..G.r1];
  for(i= G.r1+1, G.r1+G.r2,                                            \\ loop over the complex embeddings
      realvec = concat(realvec, sqrt2*real(vec1[i]));
      realvec = concat(realvec, sqrt2*imag(vec1[i]));
  );
  return(realvec);
};

\\ INPUT:
\\ - G a number field
\\ - M a (complex) matrix of dimension (r x n), where r = r1 +r2
\\ OUTPUT:
\\ - A matrix M' which has dimension n x n, obtained by splitting complex rows into real and imag parts, scaling them by sqrt2
embed_real(~G,~M)={
    my(
      outmatrix,                                                  \\ holder for the output matrix
      tempvec,                                                    \\ intermediate vector holder
      column_num,                                                 \\ holds the number of columns, should be equal to the numfield degree
      b_mat                                                       \\ holds the IB matrix
    );

    if(G.r1 == length(G.zk), return(M) );                         \\ if totally real num field, do nothing

    b_mat = M;
    column_num = matsize(b_mat)[2];                                   \\ initiate column number
    outmatrix = b_mat[1..G.r1,];
    \\outmatrix = matrix(G.r1,column_num, i, j, b_mat[i,j]);          \\ copy the real part of the matrix

    for(i= G.r1+1, G.r1+G.r2,                                         \\ loop over the complex embeddings
        tempvec = sqrt(2)*real(b_mat[i,]);
        outmatrix = matconcat([outmatrix; tempvec]);                  \\ concat to outmatrix

        tempvec = sqrt(2)*imag(b_mat[i,]);
        outmatrix = matconcat([outmatrix; tempvec]);
    );

    return(outmatrix);

}; \\ end get_R_basis


\\ INPUT:
\\ - K a number field
\\ - lglat is the output of bnf[3], which is the complex log lattice
\\ OUTPUT:
\\ This function takes in a logarithm lattice (output from bnfinit argument [3])
\\ takes the real part, and divides the rows for complex embeddings in half
\\ The point of this function is for when the bnf log lattice is computed ahead of time
process_complex_loglattice(K, lglat)={
    my(r, complexlogunit, lambda1, LambdaK);
    r = K.r1 +K.r2 -1;
    lambda1 = real(lglat);                                                       \\ equivalent to getting the log of the abs vals
    LambdaK = lambda1[1..r,];
    for(i =1, length(LambdaK),
        if(i > K.r1, LambdaK[i,] = LambdaK[i,]/2)
    );
    \\LambdaK=LambdaK*qflll(LambdaK);
    return(LambdaK);
}
\\ INPUT:
\\ - a pari bnf type (output of bnfinit)
\\ OUTPUT:
\\ - the corresponding real log lattice. Dimensions (G.r1 + G.r1) x r. Note that there is NOT a factor of 2 on the complex components, as one would obtain from real(bnf[3])
\\ (pari must take the norm of the complex values, followed by the complex log rather than the abs value)
get_log_lattice_bnf(bnf1)={
    my(r, complexlogunit, lambda1, LambdaK);
    r = bnf1.r1 +bnf1.r2 -1;
    complexlogunit=bnf1[3];                                                     \\ in pari, bnf[3] is the complex log embeddings of the independent units
    lambda1 = real(complexlogunit);                                             \\ equivalent to getting the log of the abs vals

    LambdaK = lambda1[1..r,];
    for(i =1, length(LambdaK),
        if(i > bnf1.r1, LambdaK[i,] = LambdaK[i,]/2)
    );
    \\LambdaK=LambdaK*qflll(LambdaK);
    return(LambdaK);
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ INPUT: A number field and a log lattice
\\ OUTPUT: The equivalent of the regulator for the lattice. This is different from just using matdet
\\ as it multiplies by 2 for each row corresponding to a complex embedding
unscaled_determinant(num_field, log_lat)={
    mult = max(0, num_field.r2-1);
    return( abs(matdet(log_lat))*(2^mult) );
}



\\trunc should be an integer, truncates the real up to trunc decimal digits
truncatereal(real1, trunc)={
    my(temp);
    temp = round(real1*10^trunc);
    temp = temp/(10^trunc);
    return(temp);
}


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Determines the least positive integer d such that d*idealB is an integral matrix.
\\ Input:
\\ - A matrix with rational coefficients representing an ideal B
\\ OUTPUT:
\\ - least positive integer d such that d*idealB is an integral matrix.
get_ideal_denom(idealB)={
  my(ideal_denom = 1, [row_num, col_num] = matsize(idealB));
  for(j=1, row_num,
    for(k=1, col_num,
      ideal_denom = lcm (ideal_denom, denominator(idealB[j,k]) );      \\ gets the lcm of the denominators of idealB, this is the value d(A)
    )
  );
  return(ideal_denom);
} \\ close get_ideal_denom

\\ Determines if the ideal y contains the element 1
ideal_contains1(y)={
  my(onevec, one_solution, nn=length(y));
  onevec = concat(1,vector(nn-1))~;             \\ create the vector [1, ... 0]
  one_solution = y^(-1)*onevec;                 \\ solve for vector x in the equation y*x = onevec
  for (i = 1, nn,
    if(type(one_solution[i]) != "t_INT",return(0));
  );
  return(1);                                                \\ if the function gets to here, that means there was an integer lin comb that gives 1
}

\\ Check that an ideal y is reduced: y in Q^{nxn}, USE APPROXIMATIONS m1, y1= m1*y \in R^{nxn}, qflll(y1), function in 8 and 9
\\ Does this by checking by ruling out ideals with too big of a norm, then checks that 1 is contained and is a minimum (via cubescan).
\\ Input: y is a matrix representing an ideal,
\\       G the number field
\\       eps is the error
/******************************************************************************/
checkred(y,G,eps)={
    if(abs(1/matdet(y))<=sqrt(abs(G.disc)), return(0));
    if(ideal_contains1(y)==1 && is_minimum(y, nfalgtobasis(G,1), G, eps),
        return(1),
        return(0)                                                   \\else return 0
    );
}
/******************************************************************************/
limitminvector(yy,bound)={
    my(gramy=yy~*yy);
    lis=qfminim(gramy,bound,,2)[3];
    lis;
}
checkred_old(y,G,eps)={
    if(abs(1/matdet(y))>sqrt(abs(G.disc)), return(0));
        if(ideal_contains1(y)==1,
            my(y1, lisy, n = poldegree(G.pol), gramy);
            y1=G[5][1]*y;
            y1 = embed_real(G, y1);
            y1 = y1*qflll(y1);                                                  \\ return LLL reduced y1, coeffs in terms of integral basis
            gramy = y1~*y1;
            lisy=qfminim(gramy,n,,2)[3];                                        \\ This function uses qfminim with radius n

            for (i = 1, length(lisy),
                if ( check0(y1*lisy[,i],eps) != 0, return(0))
            );
            return(1),

        0);
}

compute_sublattice(lglat, FILENAME ,extype = 0)={
    my(coord, index, pow, modpair);
    coord = random(length(lglat))+1;
    pow = random(12) +1;
    while(pow == 1, pow = random(12) +1);
    lglat[,coord] = pow*lglat[,coord];
    modpair = [coord, pow];
    write(FILENAME,"Column ", coord, " scaled by ", pow );
    if(extype == 0 || length(lglat) == 1,
        return([lglat, modpair]);
    );

    my(coord2 = random(length(lglat))+1);
    while(coord == coord2,
        coord2 = random(length(lglat))+1;
    );
    pow = random(12) +1;

    modpair = concat(modpair, [coord2, pow]);
    if(extype ==1,
        lglat[,coord2] = pow*lglat[,coord2];
        write(FILENAME,"Column ", coord2, " scaled by ", pow );
        return([lglatm, modpair]);
    );

    if(extype == 2,
        lglat[,coord] += pow*lglat[,coord2];
        write(FILENAME,"Column ", coord, " increased by ", pow, " * column ", coord2);
        return([lglat,modpair]);
    )
}

\\ given a column vector v and a lattice L, check if v is contained in the lattice
\\ eps is an error value
is_vec_in_lattice(~v,~L,eps)={
    my(v_solution, round_solution);
    GP_ASSERT_EQ(matsize(L)[1],matsize(L)[2]);
    v_solution=L^(-1)*v;                            \\ solves L*x = v for x.
    round_solution=round(v_solution);               \\ rounds the entries, this strategy is checking if x is an integer vector.
    \\
    if(norml2(v_solution-round_solution)<eps,   \\ If v is in L, then cov should be equal to covin up to some error
        1,
    return(0));                           \\ else
}

\\ element_logvec is the r+s-1 length log vector
\\ box_corners are the defining elements of an axis-aligned box
\\ RETURN:
\\ - 1 if element_logvec is contained in the box, 0 otherwise
is_in_axis_box(element_logvec, box_corners)={
    for(i =1, length(element_logvec),
        if( element_logvec[i] < box_corners[1][i] || element_logvec[i] > box_corners[2][i],
            return(0);
        )
    );
    return(1);
}

verify_lattice_containment(~G, ~new_vec)=
{
print("New vector found: Initial reg = ", precision(matdet(lattice_lambda),10));
print("WARNING: Expensive verification. Check if new found vector is truly a unit");
bnflattice = get_log_lattice_bnf(bnfinit(G));
print(precision(bnflattice^(-1)*new_vec~,10));

breakpoint();
}
