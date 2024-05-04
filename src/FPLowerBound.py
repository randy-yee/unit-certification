read("src/VectorMethods.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ LOWER BOUND METHODS (POHST-FIEKER)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
/*
Function List:
my_hermiteconstant
m_star
get_pohst_j
get_T2_value
is_non_torsion_unit
sort_by_T2
get_maximal_independent_unitset

lower_regbound

get_index_bound2
*/


\\ code for computing hermite constants for Pohst Fieker lower bound method.
\\ adapted from Pari source code
my_hermiteconstant(n)={
  /* Original C code for reference
  GEN h,h1;
  pari_sp av;
  switch(n)
  {
    case 1: return gen_1;
    case 2: return mkfrac(utoipos(4), utoipos(3));
    case 3: return gen_2;
    case 4: return utoipos(4);
    case 5: return utoipos(8);
    case 6: return mkfrac(utoipos(64), utoipos(3));
    case 7: return utoipos(64);
    case 8: return utoipos(256);
  }
  av = avma;
  h  = powru(divur(2,mppi(DEFAULTPREC)), n);
  h1 = sqrr(ggamma(sstoQ(n+4,2),DEFAULTPREC));
  return gerepileuptoleaf(av, mulrr(h,h1));
  */

  if(n == 1, return(1));
  if(n == 2, return(4/3));
  if(n == 3, return(2));
  if(n == 4, return(4));
  if(n == 5, return(8));
  if(n == 6, return(64/3));
  if(n == 7, return(64));
  if(n == 8, return(256));
  if(n > 8,
      my(h,h1);
      h = (2/Pi)^n;
      h1 = (gamma( (n+4)/2 ))^2;
      return(h1*h);
  );
};

\\ SUBROUTINE: Computes the value M* in the P-F lower bound
\\ G is a number field, const_C and jv are specific constants.
\\ Returnts the value M* from the Pohst-Fieker.
m_star(G, const_C, jv)={
  my(val1, val3, dim = length(G.zk));
  val1 = (const_C-jv)/(dim-jv);
  val3 = acosh(val1);
  return( ((dim-jv)/4)*(val3^2));
};

get_pohst_j(G)={
    return(length(G.zk)-2);
    if(G.r2 == 0,
        return(length(G.zk)-2); \\ PF suggest this is 0. But we have counterexamples
    ,\\ else
      return(min(2*G.r2, length(G.zk)-2));                                      \\ choice of j described in in F-P, p2772
    );
}

get_T2_value(G, element)={
    return( norml2(log(absoluteval_nvec(G, element) )) );
}

is_non_torsion_unit(G, candidate_elt, zerovec, eps)={
    my(zerovec = vector(G.r1+G.r2, i, 0));
    if(abs(nfeltnorm(G, candidate_elt)) != 1, return(0),
        if(!samevecs( log(abs(valuationvec(G, candidate_elt, column = 1))), zerovec, eps ), return(1), return(0) );
    );
    \\return ( abs(nfeltnorm(G, candidate_elt)) == 1 && !samevecs( log(abs(valuationvec(G, candidate_elt, column = 1))), zerovec, eps ));
}

sort_by_T2(G, unit_set)={
    my(sorted_units = [], sorted_SK_units = Mat());
    for(j =1, length(unit_set),
      sorted_units = concat(sorted_units, [[ get_T2_value(G, unit_set[j]), unit_set[j] ]]);
      \\if(DEBUG_LREG , print(precision(sorted_units[j],10), "   Norm: ", nfeltnorm(G, sorted_units[j][2]), "  powerbasis coeffs: ", G[8]^-1*sorted_units[j][2]););
    );
    sorted_units = vecsort(sorted_units);

    for(i=1, length(sorted_units),
      sorted_SK_units = matconcat( [sorted_SK_units, sorted_units[i][2]] );
    );
    return(sorted_SK_units);
}

get_maximal_independent_unitset(G, SK_units,eps)={
    my(unitmat, lunits, indextracker = [1]);
    unitmat = Mat(SK_units[,1]);
    lunits = log(abs(G[5][1]*unitmat));

    for(i=2, length(SK_units),

        \\ METHOD 1: USES APPROXIMATE MATRIX SOLVING TO DETERMINE LINEAR DEPENDENCE
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\ NOTE: Use determinant of GRAM MATRIX
        lunits = concat(lunits, log(abs(G[5][1]*SK_units[,i])));
        if( abs(matdet(lunits~*lunits)) > eps,
            indextracker = concat([indextracker, i]);
            unitmat = matconcat([unitmat, SK_units[,i]]);,
        \\else
            lunits = lunits[,1..length(lunits)-1];
        );
    );

    if(0,
        print("Independent set of units in SK: " );
        for(i=1, length(unitmat), print(unitmat[,i], "   Log embedding: ", precision(lunits[,i],10), "\n Numericals: ", precision(abs(G[5][1]*unitmat[,i]),10) ));
        print("Rank of the maximal set in SK: ", matrank(lunits), "\n");
        if(length(unitmat) ==0, write("empties.txt", G.pol));
    );

    SK_urank = length(unitmat);                                             \\ this is the value k, the size of the maximal ind. set in S_K
    return([unitmat, lunits, SK_urank]);
}

\\ INPUT:
\\ Follows the method of Fieker Pohst to get a lower bound on the Regulator
\\ - G is a number field from nfinit
\\ - K is some number larger than (1+ sqrt(2))*n
\\ - input_j is an integer value between 0 and n-2.
\\      Defaults to -1, which compute the value j as suggested by Pohst-Fieker.
\\      Otherwise use the provided value
\\ - k_limit places a maximum size on the radius to search for units,
\\      as sometimes the default value takes too long or fails
\\ OUTPUT:
\\ - A lower bound on the regulator
\\
\\ NOTE: THE DEFAULT CALCULATION OF J DOES NOT SEEM TO ALWAYS WORK,
\\ and does not match Magma for some signatures. Magma is doing something different
\\ than what is described in Pohst-Fieker which fixes the issue.
\\ Adjustments can be made using optional input variable input_j, but care needs to be taken.
\\ further investigation needed

lower_regbound(G, K, eps, input_j = -1, k_limit = 0)={
    \\if(DEBUG_LREG, print("----------------------"); print("Lowerbound: K = ", precision(K,10)); );

    my(field_deg = length(G.zk), unit_rank = G.r1+G.r2-1,
        j_val,
        K_star,
        real_IB = G[5][2],                      \\ This is a matrix G such that (Gv)^T * Gv = T2(v) \\embed_real(G, G[5][1])
        SK = [],
        SK_units = Set(),
        zerovec = vector(G.r1+G.r2, i, 0),
        inverse_element,
        SK_urank = 0,
        finalbound = 1,unitmat, lunits
    );

    if(input_j == -1,
        j_val=get_pohst_j(G);
    , \\else
        j_val = input_j;
    );

    \\if(K < (1+sqrt(2))*length(G.zk), print("Error K is too small. "); return(0));       \\ check if K is big enough
    K_star = m_star(G, K, j_val);                                                         \\ define the value K*

    \\ COMPUTE THE SET SK
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(K < k_limit || k_limit == 0,
        tt = getabstime();
        SK = qfminim(real_IB~*real_IB, K+ eps, , 2 )[3];                                    \\ get the elements whose norm is smaller than K
        tt2 = getabstime();
    , \\
        print("Pohst K is too large, searching up to the limit bound ", floor(k_limit));
        K = k_limit;
        tt = getabstime();
        SK = qfminim(real_IB~*real_IB, K+ eps, , 2 )[3];
        tt2 = getabstime();
    );
    \\ COMPUTE SK_units, THE SET OF UNITS WHICH LIE IN S_K ALONG WITH THEIR INVERSES
    for (i=1, length(SK),
        \\if ((abs(nfeltnorm(G, SK[,i])) == 1 && !samevecs( log(abs(valuationvec(G, SK[,i], column = 1))), zerovec, eps )) ==1, print("unitfound"));
        if( is_non_torsion_unit(G, SK[,i], zerovec, eps),
            \\if(DEBUG_LREG, print( "L2 Norm of TK element ", precision(norml2(((absoluteval_nvec(G, SK[,i])))), 10 ) ););
            inverse_element = vec_flip_positive(nfeltpow(G, SK[,i], -1));
            SK_units = setunion( SK_units, Set([vec_flip_positive(SK[,i]), inverse_element]) );
        );
    );
    \\tt3 = getabstime(); print("Lower Bound search time:  ",tt3-tt2 );
    SK_units = sort_by_T2(G, SK_units);         \\ Reorder SK_units as matrix columns with increasing T2(eps) value

    \\ determine a maximal independent set of units. Uses a Greedy approach to the problem
    if(length(SK_units) != 0,
        [unitmat, lunits, SK_urank] = get_maximal_independent_unitset(G, SK_units,eps);
    );

    \\ BEGIN COMPUTATION OF THE LOWER REGULATOR BOUND
    new_j = j_val;

    \\print("Fixed j value is :  ", new_j);
    for(i =1, SK_urank,                                                     \\ this is the value k, the size of the maximal ind. set in S_K

        minlog = lunits[,i]; minlog = concat(minlog, minlog[G.r1+1 .. G.r1+G.r2]);
        minlog = minlog~ * minlog;
        Mi_tilde = min(K_star, minlog );

        \\ counts number of conjugates with abs value near 1 and adjust j if possible
        \\ one_conj = 0;
        \\ numvec = abs(G[5][1]*unitmat[,i]);
        \\ for(ctr =1, length(numvec), if( abs(numvec[ctr] -1) <eps, one_conj +=1 ));
        \\ new_j = max(one_conj, new_j);
        \\ print("M_tilde_",i, " = ", precision(Mi_tilde,8), "   current j= ", new_j );

        finalbound *= Mi_tilde;
    );
    if(SK_urank < unit_rank,
      finalbound *= m_star(G, K, new_j)^(unit_rank - SK_urank);
    );
    finalbound *= (2^G.r2)/ (length(G.zk)*my_hermiteconstant(unit_rank));
    finalbound = sqrt(finalbound);

    return(finalbound);
};


\\ nf is a number field
\\ lattice_lambda is a real matrix representing a sublattice of the unit lattice
\\ eps is the error
\\ jval -1 and search limit are arguments for calculating the lower regulator bound
\\ OUTPUT is a bound in the index of lattice_lambda as a subgroup of the unit lattice
get_index_bound2(nf,lattice_lambda, eps, jval = -1, search_limit = 0, outstring = "")={
    my(constK, lrb, degree = poldegree(nf.pol));

    \\# nf.zk is specified as LLL reduced in the pari documentation
    constK = nf[5][2]*nfalgtobasis(nf, nf.zk[length(nf.zk)]); constK = constK~*constK;
    constK = 2*constK;
    numElementEstimate = (1+2*constK)^(degree);
    if(outstring != "",
        if (search_limit == 0,
            write(outstring, "K: ", constK);,
            write(outstring, "effective K ", floor(min(constK, search_limit)), " Limit: ", search_limit);
        );
    );
    lrb = lower_regbound(nf, constK, eps, jval, search_limit);
    \\write("bound2.txt", "LOGLAT: ", type(lattice_lambda),precision(lattice_lambda,10),"\nDET: ", precision(abs(matdet(lattice_lambda)),10), "  \nLowerbound: ", lrb );
    ceil(unscaled_determinant(nf, lattice_lambda)/lrb );
}
