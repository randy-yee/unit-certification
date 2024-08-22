/*
Implementation of the method of Pohst/Arenz for obtaining a set of
fundamental units from a full set of independent units.
Also known as the p-max algorithm

Notation reminder:
The matrix nf[5][1] is the matrix whose columns are the embeddings (up to complex conjugate) of integral basis elements
Sometimes we need to embed this into R^n by splitting the complex embeddings into real and imaginary parts, and scaling by sqrt(2)



Function List:
log_determinant
get_R_basis

\\ CALLS ROUTINES FROM THE file "pohst-lower-bound.py", which is all functions
\\ specifically related to building the lower_regbound function
lower_regbound

\\ The remaining function are all related to implementing the standard version of
\\ Pohst's fundamental unit from independent units method.
get_pohst_prime
get_search_ideal
get_new_eta
update_eta_set
update_eta_set_torsion
has_pth_root_normal
has_pth_root_inversion
pohst_check_normal
bruteforce_proot
LLL_reduce_units
example_gen


create_beta
get_log_lattice
#########################################
#########################################
SAGE WRAPPERS
get_integral_basis
get_units
get_disc
get_index_bound
get_column

*/

read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/FPLowerBound.py");

print("Loaded PmaxNormal functions");



\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ DEBUG flags, set to 1 if debugging a particular function.
DEBUG_MAIN = 0;
DEBUG_IDEAL_SEARCH = 0;
DEBUG_UPDATE = 0;
DEBUG_POHST = 0;
DEBUG_PROOT = 0;
DEBUG_PROOT_I =0;
DEBUG_NEW_ETA = 0;
DEBUG_LREG = 0;
DEBUG_LLLreduce = 0;

ERROR_FILE = "data/pohst-errors.txt";


\\ If the field regulator is known, can be used to see if a sublattice of the unit group is actually the unit group itself, or determine the index
\\ INPUT:
\\ - G a number field
\\ - a matrix representing a lattice in G
\\ OUTPUT:
\\ - Gives the (equivalent of the) regulator of the matrix.
log_determinant(G, unitmat)={
    my(square_matrix);
    square_matrix = G[5][1][1..length(unitmat),];
    square_matrix = square_matrix*unitmat;
    square_matrix = log(abs(square_matrix));
    for(i = 1, length(square_matrix), if(i > G.r1, square_matrix[i,]*=2; ));

    print("Field sig: ",G.r1, " ", G.r2);
    print(precision(abs(matdet(square_matrix)),10), "  " precision(abs(matdet(log(abs(get_scaled_M(G)[1..length(unitmat),]*unitmat)))),10));
    breakpoint();
    return(abs(matdet(square_matrix)));
};


\\ converts the numerical integral basis matrix of a number field with signature (r1,r2) into an appropriate
\\ matrix in R^n for the purpose of using in qfminim. n = n = r1 +2*r2
\\ Takes the r2 complex embedding values and treats their real and imag parts separately
\\ Need to scale by sqrt2 so as to maintain a consistent norm.
\\ INPUT:
\\ - G is a number field obtained from nfinit of degree n
\\ OUTPUT
\\ - An n by n square matrix as described above.
get_R_basis(G)={
    my(
      outmatrix,                                                 \\ holds output matrix
      tempvec,                                                   \\ temporary vector holder
      b_mat = G[5][1],                                           \\ holds integral basis embedding matrix
      column_num = matsize(b_mat)[2]
    );

    if(G.r1 == length(G.zk), return(b_mat) );                           \\ if totally real then do nothing
    outmatrix = b_mat[1..G.r1,];
    outmatrix = matrix(G.r1,column_num, i, j, b_mat[i,j]);              \\ copy real part of the matrix

    for(i= G.r1+1, G.r1+G.r2,                                         \\ loop over complex embeddings
        tempvec = vector(column_num, j, sqrt2*real(b_mat[i,j]));      \\ create scaled vector from real parts
        outmatrix = matconcat([outmatrix; tempvec]);                  \\ concat to outmatrix

        tempvec = vector(column_num, j, sqrt2*imag(b_mat[i,j]));      \\ create scaled vector from imag parts
        outmatrix = matconcat([outmatrix; tempvec]);
    );
    return(outmatrix);
}; \\ end get_R_basis



\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ FUNCTIONS DIRECTLY RELATED TO THE POHST METHOD OF GETTING FUNDAMENTAL UNITS
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ obtain a prime fitting specifications of Pohst/Arenz
\\ INPUT:
\\ - disc is the discriminant of a number field,
\\ - p_num is a rational prime
\\ - low is a lower bound for the resulting prime q to be found
\\ OUTPUT:
\\ - A prime q satisfying the conditions that q does not divide disc,
\\ and p_num divides (q-1)

get_pohst_prime(disc, p_num, low =1)={
    my(prime_q = nextprime(low));
    while( ((disc %prime_q) ==0) || (((prime_q-1) % p_num) != 0)  ,
        prime_q = nextprime(prime_q+1);
    );
    return(prime_q);
}; \\ end get_pohst_prime

get_pohst_prime2(disc, p_num, low =1, denom_lcm=1)={
    my(
        q_minus_one = ceil(low/p_num)*p_num,
        prime_q = q_minus_one+1
    );
    while( isprime(prime_q)==0 || ((disc %prime_q) ==0) || gcd(denom_lcm, prime_q) !=1,
        q_minus_one += p_num;
        prime_q = q_minus_one+1;
    );
    return(prime_q);
}; \\ end get_pohst_prime

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ INPUT:
\\ - G an number field,
\\ - p a prime number,
\\ - pohst_eta (See Pohst or Arenz, this is eta_k) is a number field elt (either polynomial or vector)
\\      it just needs to be a valid input for nfmodpr
\\ OUTPUT:
\\ - An ideal Q as illustrated in Pohst/Arenz
get_search_ideal(G,p, pohst_eta)={
    my(check_poly,
      field_deg = length(G.zk),
      prime_q,                                \\ holds the pohst prime
      ideal_factors,                          \\ holds ideal factorization of the prime p
      factor_candidate,                       \\ holds prime ideal
      residue_field,
      embedded_eta,
      ideal_candidate_found = 0,              \\ flag to indicate an ideal has been found
      lowbound =1,                             \\ the lower bound for get_pohst_prime
      prime_ctr = 1
    );

    \\ start of while loop which loops over primes satisfying the pohst condition until a suitable ideal is found
    while(!ideal_candidate_found,

        prime_q = get_pohst_prime(G.disc, p, lowbound);                                       \\ get the pohst prime number to be factored.
        \\print("pohst-primes: ", prime_q, "   ", get_pohst_prime2(G.disc, p, lowbound));breakpoint();
        ideal_factors = idealprimedec(G,prime_q);                                              \\ Decompose prime_q as ideals in G (vector of prid structures)

        if(DEBUG_IDEAL_SEARCH,
            print("get_ideal: Found Pohst prime q: ", prime_q); print("get_ideal: Number of Factors of q: ", length(ideal_factors));
        );

        \\ loop over all the prime ideals and check if polynomial t^p - pohst_eta = 0 is irreducible in res field
        for(i = 1, length(ideal_factors),

            factor_candidate = ideal_factors[i];                                               \\ holds the ith ideal of the decomposition, call it P

            if(DEBUG_IDEAL_SEARCH,
                res_units = idealstar(G, factor_candidate, flag=1);
                print("\n", "Prime ideal: ", factor_candidate, "Size of the residue field -1: ", res_units.no);
                \\print("bid structure: ", res_units); print("(res_units.no), " ", res_units.gen, " ", res_units.mod);
            ) ;

            residue_field = nfmodprinit(G,factor_candidate);                                    \\ get residue field OK/P
            embedded_eta = nfmodpr(G, pohst_eta, residue_field);

            \\ irreducibility check of Pohst and/or Arenz for polynomial T^p-eta (in OK/P)
            if (embedded_eta^( (factor_candidate[1]^factor_candidate.f-1)/p) != 1,
                \\print("Primes checked: ", prime_ctr);
                return(factor_candidate);                                           \\ Should always exit the function here unless the input is invalid
            );
            prime_ctr+=1;
        );                                                                      \\ end of for_loop

        lowbound = prime_q++;
    );                                                                          \\ end of while_loop

    \\ ERROR CHECK
    print( "Error: Did not find any prime ideal with the irreducibility property.");
    return(-1);
}; \\ end function get_search ideal

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Recall that prime ideals are represented as 5-component vectors: (p, a ,e,f,b)
\\ a is a coefficient vector such that p Z_k + a Z_k generates the ideal
\\ e is the ramification index, f is the residual index

\\ INPUT:
\\        Two algebraic numbers eta_i ,eta_j,
\\        A number field G
\\        A prime p,
\\        searchideal a modpr structure (derived from a prime ideal Q), i.e a prime ideal
\\ OUTPUT:
\\ - outputs the new_eta = (eta_i^t)*eta_j, with new_eta a pth power in the residue field of searchideal.
get_new_eta(G, p, k, searchideal, eta_i, eta_j, cpct_input = 0)={
    if( DEBUG_NEW_ETA, " ---- get_new_eta: Getting a new eta");
    my(
      res_field, resfield_gen,
      eta_i_mod, eta_j_mod,
      dlog_i, dlog_j, new_exp
    );

    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: The size of the res.field multiplicative group: ", searchideal.p^searchideal.f););

    res_field = nfmodprinit(G, searchideal);                                    \\ initialize the residue field

    if(cpct_input == 0,
        if(type(eta_j) == "t_FFELT",
            eta_i_mod = eta_i;
            eta_j_mod = eta_j;
        ,   \\else
            eta_i_mod =  nfmodpr(G, eta_i, res_field);                                  \\ embed the eta_i, eta_j into residue field
            eta_j_mod =  nfmodpr(G, eta_j, res_field);
        );
    ,   \\ else
        eta_i_mod = cpct_rep_modp(G, eta_i, res_field);
        eta_j_mod = cpct_rep_modp(G, eta_j, res_field);
    );
        resfield_gen = ffprimroot(ffgen(eta_i_mod));                                \\ gets a multiplicative generator

    if(DEBUG_NEW_ETA, print(" ---- get_new_eta: Is ", resfield_gen , " really a generator? ", resfield_gen.p - resfield_gen^((searchideal.p^searchideal.f -1)/2) ););

    dlog_i =  fflog(eta_i_mod, resfield_gen);                                   \\ get the DLOG wrt the mult. generator
    dlog_j = fflog(eta_j_mod, resfield_gen);

    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: DLOGS of eta_i, eta_j  ", dlog_i, "  ,   ", dlog_j ););

    new_exp = lift(-dlog_j*Mod(dlog_i,p)^(-1));
    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: Computed eta_",k," power: ",new_exp ););

    return(new_exp);                                                            \\JL7_UPDATE
};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ This is a subroutine for the pohst search function.
\\ INPUTS:
\\    A number field G
\\    A prime p
\\    eta is a matrix of independent units
\\    k is the index of the 'first' element, which the later elements are modified by
\\    eta_k is the modifying element
\\    eta_exp_mat keeps track of the power of eta_k which we modify by
\\    new_units is the matrix of the current independent set of units.
\\ OUTPUT:
\\ - The new set of eta, which keeps the first k elements the same, and the rest are modified.
update_eta_set(G, k, eta_k, p, eta_exp_mat, new_units)={
    my(eta_i, searchideal, res_field, LLLnewmat = Mat(), exponent_tracker=[], exp_i);

    searchideal = get_search_ideal(G,p, eta_k);                                 \\ obtain the search ideal
    res_field = nfmodprinit(G, searchideal);
    eta_k = nfmodpr(G, eta_k, res_field);

    for(i=k+1, length(new_units),
        eta_i = nfmodpr(G, 1, res_field);
        for(j=1, length(eta_exp_mat), eta_i *= (nfmodpr(G,new_units[,j], res_field)^eta_exp_mat[j,i]) );

        if(DEBUG_UPDATE, print("Update_eta: Updating eta_",i ););
        exp_i = get_new_eta(G,p, k, searchideal, eta_k, eta_i);                 \\ find the power of eta_k such that eta_k^t * eta_i is a pth power
        exponent_tracker = concat(exponent_tracker, exp_i);                     \\ JL7_UPDATE

    );
    if(DEBUG_UPDATE, print("exponent tracker ", exponent_tracker););            \\ JL7_UPDATE
    return([exponent_tracker]);
};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ This is a subroutine for the pohst search function. Updates the units by the torsion element
\\ INPUTS:
\\    A number field G
\\    A prime p
\\    eta is a matrix of independent units
\\    k is the 'first' element, which the later elements are modified by
\\ OUTPUT:
\\ - The new set of eta, which keeps the first k elements the same, and the rest are modified.
update_eta_set_torsion(G, t_gen, p, new_units)={
    my(eta_i, searchideal, LLLnewmat = Mat(), exponent_tracker=[], exp_i);

    t_gen = nfbasistoalg(G, t_gen);                                             \\ turn t_gen to poly representation
    searchideal = get_search_ideal(G,p, t_gen);                                 \\ obtain the search ideal

    for(i=1, length(new_units),
        eta_i = new_units[,i];
        eta_i = nfbasistoalg(G, eta_i);                                         \\ turn eta_i to poly representation

        if(DEBUG_UPDATE, print("Update_eta: Updating eta_",i ););
        exp_i = get_new_eta(G,p, 0, searchideal, t_gen, eta_i);                 \\ find the power of t_gen such that t_gen^t * eta_i is a pth power

        exponent_tracker = concat(exponent_tracker, exp_i);                     \\ JL7_UPDATE
    );
    if(DEBUG_UPDATE, print("Update_eta: exponent tracker ", exponent_tracker););            \\ JL7_UPDATE

    return([exponent_tracker]);
};
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ This is a subroutine for the pohst_check.
\\ INPUT: Number field G,
\\          p a rational prime
\\          unit_coeffs is a column vector representing an independent unit eta
\\          eps is a global variable indicating the error
\\ OUTPUT: Return a pair consisting of:
\\        a flag indicating whether a solution has been found,
\\        A solution or 0

has_pth_root_normal(G, real_IB, p, unit_coeffs,eps)={
    my(
      basis_numerical,
      unit_numerical,
      scaled_basis,
      quad_form,
      minima_list,
      normcheck,
      checkflag, candidate, candidate1, change_of_basis
    );
    basis_numerical = real_IB;

    unit_numerical=absoluteval_nvec(G, unit_coeffs, side = 2)^(-1/p);           \\ Ha's idea, get the norms before turning it into a real matrix

    scaled_basis = matdiagonal(unit_numerical)*basis_numerical;                 \\ Ha's idea, get the norms before turning it into a real matrix
    change_of_basis = qflll(scaled_basis);
    LLLscaled_basis = scaled_basis*change_of_basis;
    quadform = (LLLscaled_basis~)*LLLscaled_basis;                                    \\ get the quadratic form (as a symmetric matrix)

    if(DEBUG_MAIN,print("--- p_root: Entering qfminim search on length: ", length(real_IB)););
    minima_list = qfminim(quadform, length(real_IB)+eps, 100, 2);                     \\ length(real_IB)) + eps is the target norm + some error
    if(DEBUG_MAIN, print("--- p_root: Number of elements found: ", minima_list[1]); );

    if (minima_list[1] != 0,                                                    \\ check if the list returned was empty

        \\ check all elements found for one whose power yields the original unit
        for(i=1, length(minima_list[3]),

            candidate1 = change_of_basis*minima_list[3][,i] ;
            \\print("candidate: ", precision(candidate1,10));
            if(abs(nfeltnorm(G, candidate1))!=1, ,
                candidate = nfeltpow(G,candidate1, p);
                checkflag = 1;

                if(DEBUG_PROOT, print("--- p_root: p-power quotient yields root of unity: ", precision(log(abs(nfeltembed(G, nfeltdiv(G, candidate, unit_coeffs)))),10) ););

                for(j=1, length(candidate),
                    if( abs(abs(candidate[j])- abs(unit_coeffs[j])) >= eps, checkflag = 0;);
                );
                if( checkflag,
                    if(DEBUG_PROOT, print("--- p_root: FOUND VALUE raised to p = ", nfeltpow(G,minima_list[3][,1], p ) ););
                    return([1, candidate1 ]);
                );
            );
        );
        if(DEBUG_MAIN, print("Found values all failed the check"););
        return([0, 0]);

        );
    return([0,0]);

};


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ SUBROUTINE: of has_pth_root_inversion
\\ This function allows you to compare the log embeddings of elements even when
\\ one of the entries is *incorrectly* 0 (due to precision errors)
\\ INPUT:
\\ - G a number field
\\ - factor is an integer by which to scale the first element
\\ - el1 el2 are elts of G as coefficient vectors
\\ OUTPUT:
\\ returns 1 if factor * logembedding(el1) == logembedding(el2), 0 otherwise
compare_logs(G, fac, el1, el2,eps)={
    my(zeroflag = 0, el1_num = abs(G[5][1]*el1), el2_num = abs(G[5][1]*el2));

    for(i=1, length(el1_num),
        if(el2_num[i] == 0, zeroflag = 1;
        , \\else
            if( abs( log(el2_num[i])-fac*log(el1_num[i]) ) > eps, return(0); );
        );
    );
    if(zeroflag > 0, print("WARNING: Log comparison may be inaccurate"));
    return(1);
}

get_dual_basis(G)={
    my(dualbasis, complexdupes);
    if(G.r2 != 0,
        complexdupes = G[5][1][G.r1+1.. G.r2+G.r1,];                                \\ copy complex rows, append conjugate rows to bottom
        dualbasis = matconcat( [ G[5][1], conj(complexdupes)]~);
    ,\\else
        dualbasis = G[5][1];
    );
    dualbasis = abs(dualbasis^(-1));                                            \\ obtain the numerical dual basis matrix
    return(dualbasis);
}

\\ compute a possible solution with coefficients modulo the main_prime
get_candidate_p_root(G, main_prime,unit_coeffs,p_inverse)={
    local(candidate);
    candidate = nfeltpow_modp(G, main_prime, unit_coeffs, p_inverse);           \\ compute a candidate solution
    if (type(candidate) == "t_INT", candidate = concat([candidate], vector(length(G.zk)-1,i, 0))~; );
    for(i =1, length(candidate),
        if(candidate[i] > main_prime/2, candidate[i] = candidate[i] - main_prime;);
    );
    return(candidate);
}

get_maximum_T(G, unit_coeffs, p)={
    my(dualbasis, pth_root, T, Tmax);
    dualbasis = get_dual_basis(G);
    pth_root = absoluteval_nvec(G, unit_coeffs, side = 1)^(1/p);                \\ given the unit eta, obtain vector |eta_i|^1/p
    T = dualbasis*(pth_root~);                                                  \\ get vector of T_i
    Tmax = T[1];
    for(i=2, length(T), Tmax = max(Tmax, T[i]) );                               \\ Determine the max value
    if(DEBUG_PROOT_I, "max Ti = ",print(precision(Tmax,10)););
    return(Tmax);
}
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ This is a subroutine for the pohst_check.
\\ INPUT: Number field G,
\\          p a rational prime
\\          unit_coeffs is a column vector representing an independent unit eta
\\          eps is a global variable indicating the error
\\ OUTPUT: Return a pair consisting of:
\\        a flag indicating whether a solution has been found,
\\        A solution or 0
has_pth_root_inversion(G, p, unit_coeffs,eps)={
    if(DEBUG_PROOT_I, print("----  has_pth_root_inv BEGIN ----"); );
    my(Tmax, main_prime,
        prime_conditions,                                                       \\ flag for getting suitable prime
        hp_primes,                                                              \\ to hold set of primes dividing the h_P for a candidate P
        p_inverse,
        hp,
        componentsize,
        candidate
    );

    \\ STEP 1: Compute the dual basis to the integral basis and the values Ti
    Tmax = get_maximum_T(G, unit_coeffs, p);

    \\ STEP 2: Get the big prime P
    main_prime = nextprime(ceil(2*Tmax));                                       \\ chose a prime big enough in relation to the T_i
    prime_conditions = 0;                                                       \\ set the condition flag to false
    hp = 1;

    while(!prime_conditions,
        if(G.disc % main_prime != 0,                                            \\ if P does not divide discriminant
            hp_primes = Set(); hp =1;
            p_factors = idealprimedec(G,main_prime );                           \\ factors P in G

            for(i =1, length(p_factors),
                \\print("prime and residue degree ", p_factors[i][1], "^", p_factors[i][4], "-1 = ", p_factors[i][1]^p_factors[i][4] -1);
                componentsize = p_factors[i][1]^p_factors[i][4] -1;
                hp *= componentsize^p_factors[i][3];                            \\ p_factors[i][1]^p_factors[i][4]
                hp_primes = setunion( hp_primes, Set( factor(componentsize)[,1]));
            );

            prime_conditions = !setsearch(hp_primes, p);                        \\ set to 0 if mainprime is in hp_primes, 1 if not
        );
        if(!prime_conditions, main_prime = nextprime(main_prime+1); );
    );
    \\print("mainprime P: ", main_prime, "\n", hp_primes);

    \\ STEP 3: determine the inverse of p mod P
    CRTmods = [];
    for(i=1, length(hp_primes), CRTmods = concat(CRTmods, valuation(hp, hp_primes[i])));

    if (DEBUG_PROOT_I,
        p_inverse = 1; for(j = 1, length(hp_primes), p_inverse *= hp_primes[j]^CRTmods[j]);   \\ reconstruct hp from the primes and valuations
        \\print(" hp = reconstructed? ", p_inverse == hp);
    );

    for ( j =1, length(CRTmods),
        CRTmods[j] = hp_primes[j]^CRTmods[j];
        hp_primes[j] = Mod(p, CRTmods[j])^-1;
    );

    p_inverse = lift(chinese(hp_primes));
    if(DEBUG_PROOT_I,  print("inverse check:  ", p_inverse == lift(Mod(p, hp)^(-1) ) ); );

    \\ STEP 4: Determine if there is a pth root
    candidate = get_candidate_p_root(G, main_prime,unit_coeffs,p_inverse);

    if(DEBUG_PROOT_I ,
      print("Candidate root: ", candidate);
      print("Expected power: ", unit_coeffs);
      print("Candidate raised to p ", nfeltpow(G, candidate, p));
      print( precision(p*log(valuationvec(G, candidate, column = 1)) ,20), "\n", precision(log(valuationvec(G, unit_coeffs, column = 1)) ,20)  );
      print("log vector difference: " precision( p* log(valuationvec(G, candidate,column =1)) - log(valuationvec(G, unit_coeffs, column = 1)) ,20) );
    );
    if(DEBUG_PROOT_I, print("----  has_pth_root_inv END ----"););

    \\write(ERROR_FILE, "\n ", compare_logs(G, p, candidate, unit_coeffs,eps),"  " ,samevecs(p*log(valuationvec(G, candidate, column = 1)), log(valuationvec(G, unit_coeffs, column = 1)), eps) );
    if( compare_logs(G, p, candidate, unit_coeffs,eps) ,
      print("pth_root_inv found a solution!");
      return([1,candidate]),
    \\else
      return([0,0]);
    );

}; \\ close the alternate pth root method



\\ For now assume that the unit lattice is provided as coefficients in terms of the integral basis
\\ B is the bound of primes to check for, should be computed using some function based on the regulator.
pohst_check_normal(G, L, B,eps)={
    print("---- Running pmax_normal ----");
    my(
        index_holder=1,                                                         \\ keeps track of what the index currently is believed to be
        eta_matrix,
        new_units = L,                                                          \\ copy L into the eta matrix
        k,
        solutionflag,
        solution,
        index_bound = B,                                                        \\ upper bound on index
        betavec,
        eta_exp_mat,                                                            \\ JL7_UPDATE eta_coeff, update variables added
        updatevec,
        [torsion, torsiongen] = nfrootsof1(G),
        p = 2,
        real_IB = get_R_basis(G);                                               \\ holds the integral basis real-embedding matrix
    );
    torsiongen2 = nfalgtobasis(G,torsiongen);


    \\ Loop over primes which may divide the index
    while(p < = floor(index_bound),

        \\ construct the input to pari_prime_check
        betavec = [];
        if(torsion %p == 0, betavec = concat(betavec, [torsiongen2]));
        for(i=1, length(new_units), betavec = concat(betavec, [new_units[,i]]) );

        if(pari_prime_check(G, betavec, p, 1,0) == 0,

            \\ handling the torsion:
            if(gcd(p, torsion) != 1,
                if(DEBUG_POHST, print("Update torsion: ", p, "   ", torsion););
                [updatevec] = update_eta_set_torsion(G, torsiongen, p, new_units);
                if(DEBUG_POHST, print("torsion exponent vec: ", updatevec););

                for(i =1, length(new_units), new_units[,i] = nfeltmul(G, nfeltpow(G, torsiongen, updatevec[i]), new_units[,i]) );
            );

            \\eta_matrix = new_units;
            eta_exp_mat = matid(length(new_units));                          \\ JL7_UPDATE
            if(DEBUG_MAIN,
            print("--- MAIN LOOP DATA ----","\n--- P = ", p);
            for(j=1, length(new_units), print("--- eta_", j , " = ", new_units[,j], " Norm = ", nfeltnorm(G, new_units[,j]) ) );
            print("--- Current lattice log determinant: ", precision(log_determinant(G, new_units), 10));
            print("--- MAIN LOOP DATA END ----");
            );
            \\eta_matrix = new_units;                                                 \\ needed to reset the eta_matrix so that it doesn't get too big.

            k = 1;                                                                  \\ Each time a new prime is considered, reset k to 1
            solutionflag = 1;                                                       \\ flag indicating if a solution is found
            solution = 0;                                                           \\ variable to hold solution

            while(solutionflag,

                if(DEBUG_MAIN, print("\n--- For P= ", p, " Checking eta_", k, " for pth roots..."););

                \\ compute test_eta_k prior to input to pth root, there's no way to avoid this
                test_eta_k = 1;
                for(i = 1, length(new_units), test_eta_k = nfeltmul(G, test_eta_k, nfeltpow(G, new_units[,i], eta_exp_mat[i,k])));

                if(p < length(G.zk)+1,
                    [solutionflag, solution] = has_pth_root_normal(G, real_IB, p, test_eta_k, eps);,
                \\else
                    [solutionflag, solution] = has_pth_root_inversion(G, p, test_eta_k, eps);
                );

                if(DEBUG_MAIN, print("--- ( P= ", p, " eta_", k, ") Found a solution?: ", solutionflag , "\n"););

                \\ when a solution is found, execute the following code
                if(solutionflag,
                    if(DEBUG_MAIN, print("--- pth root found for p = ", p ); print("found solution: ", solution); );
                    \\print("Current units: ", new_units);

                    \\ UPDATE new_units by replacing the kth entry and performing LLL reduction on the real matrix
                    new_units[,k] = solution;
                    \\print("after replacement: ", precision(log(abs(K[5][1]*new_units)),10));

                    new_units = LLL_reduce_units(G, new_units);
                    print("Regulator of new_units after replacement+LLL: ", precision(log_determinant(G, new_units), 10));


                    \\ REACCOUNT for the torsion after performing LLL
                    if(gcd(p, torsion) != 1,
                        \\print("Update torsion - post-LLL: ", p, "   ", torsion);
                        [updatevec] = update_eta_set_torsion(G, torsiongen, p, new_units);
                        for(i =1, length(new_units), new_units[,i] = nfeltmul(G, nfeltpow(G, torsiongen, updatevec[i]), new_units[,i]) );
                    );

                    \\ UPDATE INDEXBOUND AND CHECK TERMINATION CONDITION
                    index_bound = ceil(index_bound/p);
                    index_holder = index_holder*p;
                    k = 1;
                    if(index_bound ==1, break);

                );  \\ closes solution found code
                \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                if(solutionflag == 0,                                               \\ in the case no solution to the equation, we have to update the eta matrix
                    if(k == length(L), solutionflag = 0; break );                   \\ when k == length(L), we are finished with the current prime.
                                                                                    \\ break out of the while loop and go to the next prime p
                    if(DEBUG_MAIN,
                        print("eta update inputs: ");
                        print(k, "  ", test_eta_k, "\n\n", eta_exp_mat, "\n\n", new_units );
                      );
                    [updatevec] = update_eta_set(G,k,test_eta_k,p,eta_exp_mat ,new_units);
                    for(j = k+1, length(eta_exp_mat),
                        eta_exp_mat[,j] = eta_exp_mat[,j] + updatevec[j-k]*eta_exp_mat[,k];         \\ add the correct multiple of the kth column
                    );

                    eta_exp_mat = eta_exp_mat %p ;                                                  \\ reduce coeffs mod p
                    \\ ensure exponents lie in the range [-p/2, p/2]
                    for(bb =1, matsize(eta_exp_mat)[1],
                        for(aa =1, matsize(eta_exp_mat)[2],
                          if( eta_exp_mat[aa,bb] > p/2, eta_exp_mat[aa,bb] = eta_exp_mat[aa,bb] -p );
                    ));

                    k++;

                    solutionflag = 1;
                );  \\ close the no solution found code
                \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            );,                                                                 \\ close while loop and start the 'else' of the pari_prime_check 'if'
            if(DEBUG_POHST,print("pari check succeeds"));
        );                                                                      \\ close the pari_prime_check 'if statement'
        p = nextprime(p+1);
    );

    print("Finished!");
    return(new_units);     \\ Note, the modified matrix does indeed generate the right lattice -- See Arenz pg 165

};



\\ Used to verify the results of the pohst-check-normal via a brute force strategy.
\\ INPUT:
\\ - G a number field
\\ - p a prime
\\ - unitmat a matrix representing a set of r independent units
\\ OUTPUT:
\\ - Runs through every exponent combination from 0..p of the elements in unitmat. Uses the has_pth_root_normal to determine if there is a pth root.
\\
bruteforce_proot(G, p , unitmat, eps)={
  print("Brute force check: ");
  my(exponentvec = vector(length(unitmat), i, 0));                              \\ initialize exponent vec to 0 vector


  for(i=1, p^length(unitmat)-1,                                                 \\ there are p^r-1 exponent vectors to consider, excluding 0
    exponentvec[1] += 1;
    for(j=1, length(exponentvec)-1,                                             \\ the loop updates the exponent vector by incrementing by 1 as a base p number
        if(exponentvec[j] == p,
            exponentvec[j+1] += 1;
            exponentvec[j] = 0;
        );

    );
    \\print(exponentvec);
    test_element = nfalgtobasis(G,1);                                           \\ get the element 1

    \\ this for loop computes the element with the corresponding exponents on each element of unitmat
    for (j=1, length(unitmat),
        \\print("Components: ", nfbasistoalg(G, test_element), "  ", nfbasistoalg(G,unitmat[,j]), "  ", exponentvec[j]);
        test_element = nfeltmul(G, test_element, nfeltpow(G, unitmat[,j], exponentvec[j]));
        \\print(test_element, "    ", nfbasistoalg(G, test_element) "\n");
    );
    real_IB = G[5][2];
    print(test_element, "   ", exponentvec);
    print(has_pth_root_normal(G, real_IB, p, test_element ,eps));               \\ checks if pth root exists.

  );
};



\\ INPUT:
\\ - G a number field
\\ - a matrix of independent units, as coefficients of the integral basis
\\ OUTPUT:
\\ - A matrix of units whose corresponding log lattice is LLL reduced
LLL_reduce_units(G,uLat)={
    my(loglat, LLL_cb_matrix, tempvec, newLat);

    if(length(uLat)==1,
      LLL_cb_matrix = Mat([1]),
    \\ else
      loglat = log(abs(G[5][1]*uLat));                                             \\ get the log lattice of uLat
      LLL_cb_matrix = qflll(loglat);
      \\print(precision(loglat,10));print("LLL transform mat: ", LLL_cb_matrix);
    );

    newLat = matrix(matsize(uLat)[1], matsize(LLL_cb_matrix)[2],i,j, 0 );

    \\ outer for loop for computing the ith unit of the LLL reduced units
    for(i =1, length(LLL_cb_matrix),
        tempvec = nfalgtobasis(G, 1);                                           \\ initialize the new unit as 1

        \\ inner loop
        for(j=1, length(uLat),
            if(DEBUG_LLLreduce, print(i,"  ", j, "  ", LLL_cb_matrix[j,i]  ););
            \\print("powering: ", uLat[,j], "   ", LLL_cb_matrix[j,i], " = ", nfeltpow(G, uLat[,j], LLL_cb_matrix[j,i]));
            tempvec = nfeltmul(G, tempvec, nfeltpow(G, uLat[,j], LLL_cb_matrix[j,i]) );
            if(tempvec ==1, tempvec = nfalgtobasis(G,1));
        );

        newLat[,i] = tempvec;
    );

    return(newLat);
};



\\ INPUT:
\\ K nfinit of a polynomial f
\\ K1 bnfinit of a polynomial f
\\ rsize an integer indicating the max scaling factor of units
\\ COMBINATION is a flag that indicates if you want one unit to be a power product of the others, or to simply power each unit individually

\\ OUTPUT
\\ coefficient matrix of a sublattice of the fundamental units, determined by the random exponent vector
example_gen(K,K1, rsize, COMBINATION)={
  print("------ EXAMPLE DATA -------");
  my(unitrank, randexp, power_units, tempunit);
  unitrank = length(K1.fu);
  randexp = vector(unitrank,i, random(rsize)+1);
  \\randexp = vector(unitrank,i, rsize);
  power_units = [];
  COMBINATION = 0;

  if (unitrank >1 && COMBINATION,
      tempunit = 1;
      for (j=1, length(K1.fu), tempunit = tempunit* K1.fu[j]^(randexp[j]));     \\ This makes the first unit a power product of all the fundamental units
        power_units = concat(power_units, tempunit);
      for (j=2, length(K1.fu), power_units = concat(power_units, K1.fu[j]));,   \\ the remaining units all stay the same.
  \\else
      power_units = vector(length(K1.fu), i, K1.fu[i]^randexp[i] );             \\ otherwise just power each unit to the ith coord of the randexp vector
  );

  punit_matrix = units_to_matrix(K, power_units);

  print("Random exponent vector: ", randexp);
  print("Generated sublattice:"); for(i=1, length(punit_matrix), print(punit_matrix[i,]) );
  print("Sublattice *regulator*: ", precision( log_determinant(K,punit_matrix), 10));
  print("------ END EXAMPLE DATA ------- \n");
  return(punit_matrix);
};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ SAGE WRAPPER CLASSES

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
get_integral_basis(nf)={return(nf.zk)}
get_degree(nf)={return(poldegree(nf.pol))}
get_units(bnf)={return(bnf.fu)}
get_disc(nf)={return(nf.disc)}
get_pol(nf)={return(nf.pol)}
get_index_bound(bnf,eps)={
    my(constK);
    constK = bnf.nf[5][2]*nfalgtobasis(bnf, bnf.zk[length(bnf.zk)]); constK = constK~*constK;
    constK = 2*constK;
    ceil(log_determinant(bnf.nf, units_to_matrix(bnf.nf, bnf.fu))/lower_regbound(bnf.nf, constK, eps) );
}


get_column(gpmat, val)={
    return(gpmat[,val]);
}
replace_column(M1, k, col1)={
    M1[,k] = col1; return(M1);
}
divide_vec(vec, den)={
    return(vec/den);
}
compute_K(nf)={
    my(constK);
    constK = nf[5][2]*nfalgtobasis(nf, nf.zk[length(nf.zk)]); constK = constK~*constK;
    constK = 2*constK;
    return(constK);
}
index_bound_from_lower(latticelambda, lbound)={
    return(ceil(abs(matdet(latticelambda))/lbound));
}

\\ INPUT:
\\ - torsion is the size of the torsion group
\\ - torsiongen2 is a generator for the torsion group
\\ - p is a prime
\\ - new_units is the set of independent units
\\ - M a matrix of independent units, as coefficients of the integral basis (r x n)
\\ OUTPUT:
\\ - A vector beta for use in the pari_prime_check function
create_beta(torsion, torsiongen2, p, new_units)={
    my(betavec = []);
    if(torsion %p == 0, betavec = concat(betavec, [torsiongen2]));
    for(i=1, length(new_units), betavec = concat(betavec, [new_units[,i]]) );
    return(betavec);
}

\\ INPUT:
\\ - nf1 a number field
\\ - M a matrix of independent units, as coefficients of the integral basis (r x n)
\\ OUTPUT:
\\ - The corresponding rxr log lattice matrix, where r is the unit rank

\\ NOTE: Implemented for testing. In practice it may be beneficial to add an LLL step at the end
get_log_lattice(nf1, M)={
    r = nf1.r1 +nf1.r2 -1;
    numMat = log(abs(nf1[5][1]*M));
    numMat = numMat[1..r,];
    \\ numMat = numMat*qflll(numMat);
    return(numMat);
}
