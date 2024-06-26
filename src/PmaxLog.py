
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ IMPORTS

print("Loaded Log-Pohst Functions");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ IMPORTS
read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/CompactRepresentation.py");
read("src/Neighbours.py")
read("src/bounds.gp")
read("src/PmaxNormal.py");

/*
PMAX LOG FUNCTION LIST:
The functions below are used to implement the logarithm version of
Arenz's method / Pohst's On Computing Fundamental Units
in combination with the pari/gp prime_check method
For obtaining a set of fundamental units from a full set
of independent units.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
POHST ALGORITHM (LOG VERSION)


get_scaled_basis
check_in_unitlattice
complex_check_in_unitlattice
get_search_ideal_cpct
get_new_eta_cpct
update_eta_set_log
lcm_denominators
log_pohst_pari

*/
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Global variables
sqrt2 = sqrt(2);



\\ DEBUG FLAGS
DEBUG = 0;
DEBUG_CPCT = 0;
DEBUG_BSGS = 0;
DEBUG_MLLL = 0;
DEBUG_NEW_ETA = 0;
DEBUG_MEMBERSHIP = 0;
DEBUG_LPOHST = 0;
TIMING = 0;



\\ START THE FUNCTION IN THE LOG-POHST FILE

/******************************************************************************/
/******************************************************************************/


\\ subalgorithm of check_in_unitlattice. Embeds an ideal of G into R^n,
\\ scales the lattice by v and returns the LLL reduced form
\\ INPUTS:
\\ - G is a number field,
\\ - v the embedding vector of an element of G (length is G.r1+G.r2, no abs values taken yet)
\\ - videal is the coefficient matrix of an ideal
get_scaled_basis(G, v, videal)={
    my(
        abs_v = abs(v),
        real_ideal_basis = G[5][1]*videal,
        change_of_basis,
        real_vec = abs_v[1..G.r1];
    );
    if(G.r2 != 0,                                                               \\ if there are complex embeddings, convert to real matrix (nxn) and real n-vector
        real_ideal_basis = embed_real(G, real_ideal_basis);
        for(i=G.r1+1, G.r1+G.r2,
            real_vec = concat(real_vec, [abs_v[i],abs_v[i] ]);
        );
    );
    real_ideal_basis = matdiagonal(real_vec)*real_ideal_basis;                  \\ scale the rows
    change_of_basis = qflll(real_ideal_basis);
    real_ideal_basis = real_ideal_basis*change_of_basis;                        \\ perform LLL

    return([real_ideal_basis,change_of_basis] );                                                \\ return scaled basis
}

\\ This is used for checking whether an element is part of the log lattice Lambda_K
\\ given that we only have the log image of a potential element.
\\ For us, this will be (1/p)eta, where eta is a log image obtained via the pohst algorithm.
\\ as such, we also know that the missing coordinates of v becaues they must sum to 0
\\ INPUT:
\\ - G an number field,
\\ - v a log vector
\\ OUTPUT:
\\ - 1 or 0 indicating True or False
check_in_unitlattice(G, v, eps)={
    \\# force v to be a row
    if(type(v) == "t_COL", v = v~);
    if(norml2(v)<eps^2, return(1););  \\# v is zero, automatically return 1

    my(
        m1=G[5][1], n=poldegree(G.pol), r=G.r1+G.r2-1,
        radius_S, jump_output, log_mu, new_y, log_mu_p, babystock, babystock_log, scaled_basis, quadform, change_of_basis
    );

    radius_S = 1/4*sqrt(r)*log(abs(G.disc));                                \\# This is the search radius S
    \\#jump algorithm determines nearby divisor to v (Alg 8 - Step 2)
    jump_output = giantstep(matid(n),v,G,n,eps);                            \\# obtain a minima mu that is 'close' to v using JUMP

    \\debug_compare(jump_output, giantstep_high_precision(matid(n),v,G,n,eps));

    extra_coord = 0;
    for(i=1, length(v),
        if (i<=G.r1,
            extra_coord-=v[i];,
            extra_coord-=2*v[i];
        );
    );
    if (G.r2 != 0, extra_coord/=2;);
    v2= concat(v, extra_coord);

    log_mu = jump_output[3];                                                \\# this is a log vector of mu,  Psi(mu) in the Fonteine-Jacobson paper
    new_y = jump_output[1];                                                 \\# this is the reduced ideal (1/mu)*y

    \\# Alg 8 - Step 3
    log_mu_p = v - log_mu[1..r];
    \\log_mu_p = v - log_mu;
    log_mu_p2 = v2 - log_mu;

    print("log_mu_p  ", precision(log_mu_p,10), "   ", precision(log_mu_p2,10));
    if(DEBUG_MEMBERSHIP,
    print("\n Original membership check: \n \
        - Checking if v is in unit lattice, v: ", precision(v_extra,10),
        "\n- psi(mu)=log_mu: ", precision(log_mu,10), \
        "\n- Norm(v - psi(mu)) = ", precision(norml2(log_mu_p),10) );
    );

    \\# return 0 if too far away, return 1 if we get 0 since implies exact match
    if(norml2(log_mu_p2) > radius_S^2, return(0));
    if(norml2(log_mu_p2) < eps,
        print("  Element of unit lattice found");
        return(1);
    );

    \\# step 7, which exponentiates the value, then computes matrix product
    \\# which represents lattice in which we will enumerate
    exp_log_mu_p2 = real(exp(log_mu_p2)); print("exp_log_mu_p2 ", precision(exp_log_mu_p2, 10));
    exp_log_mu_p = create_target(G, log_mu_p); \\#print(" l2norm of log_mu' = ", precision(norml2(exp_log_mu_p),10));
    print("exp_log_mu_p", precision(exp_log_mu_p, 10));
    [scaled_basis, change_of_basis] = get_scaled_basis(G, exp_log_mu_p2, new_y);

    quadform = scaled_basis~*scaled_basis;

    babystock = qfminim(quadform, n+eps, , 2)[3];
    if(length(babystock) == 0,
        \\return(0)
    ,\\else
        for(i=1, length(babystock),
            candidate = new_y*change_of_basis*babystock[,i];
            if(nfeltnorm(G, candidate) != 1, ,
                candidate = log(abs(nfeltembed(G, candidate)));
                print("Compare elements: ", precision(candidate, 10), "   ", precision(log_mu_p,10));
            );
        );
    );
    for(i=1, length(babystock),
        babystock_log = log(abs(valuationvec(G,new_y*change_of_basis*babystock[,i], column=1)));
        babystock_log = unsquare_log_embeddings(G, babystock_log);
        \\print( "  UNITLATTICE CHECK babystock log", precision(babystock_log, 10),
        \\            "\n log_mu_p ",precision(log_mu_p,10),"\n difference: ", precision(abs(log_mu_p) - abs(babystock_log[1..r]),10));
        if(norml2(abs(log_mu_p) - abs(babystock_log[1..r])) < eps && idealdiv(G, matid(n), babystock[i]) == new_y,
            print("  UNIT_LATTICE_CHECK SUCCEEDS");
            return(1);
        );
    );

    return(0);
}

exponentiated_embedding(G, logvec)={
    elmp = vector(poldegree(G.pol), i, 0);
    for(i=1, G.r1, elmp[i]=exp(real(logvec[i])));
    for(i=1, G.r2,
        val = exp(logvec[i]);
        elmp[G.r1+2*i-1] = val;
        elmp[G.r1+2*i] = conj(val);
    );
    return(elmp);
}

make_embedding_matrix_square(G)=
{
    res = G[5][1][1.. G.r1,];
    for(i=1, G.r2,
    res = matconcat([res~, G[5][1][G.r1+i,]~])~;
    res = matconcat([res~, conj(G[5][1][G.r1+i,])~])~;
    );
    return(res);
}

\\ This is used for checking whether an element is part of the log lattice Lambda_K
\\ given that we only have the complex log image of a potential element.
\\ INPUT:
\\ - G an number field,
\\ - v a complex log vector
\\ - eta (See Pohst or Arenz, this is eta_k) as a polynomial
\\ OUTPUT:
\\ - 1 or 0 indicating True or False
complex_check_in_unitlattice(G, v, eps)={
    if(type(v) == "t_COL", v = v~);
    if(norml2(v)<eps^2,
        print("close to zero");
        return(1),                          \\# v is zero, automatically return 1
    \\else
        M_inv = make_embedding_matrix_square(G)^(-1);

        my(
            m1=G[5][1], n=poldegree(G.pol), r=G.r1+G.r2-1,
            radius_S, jump_output, log_mu, new_y, log_mu_p, babystock,
            babystock_log, scaled_basis, quadform, change_of_basis
        );

        \\ use jump to get a nearby divisor
        radius_S = 1/4*sqrt(r)*log(abs(G.disc));                                            \\ This is the search radius S
        jump_output = giantstep_high_precision(matid(n),real(v[1..r]),G,n,eps,1);           \\ obtain a minima mu that is 'close' to v using JUMP

        log_mu = jump_output[3];                                                \\ this is a log vector of mu,  Psi(mu) in the Fonteine-Jacobson paper
        new_y = jump_output[1];                                                 \\ this is the reduced ideal (1/mu)*y
        log_mu_p = v - log_mu;                                                  \\ this is log(mu') as in the FJ paper.
        print("2 v       ", precision(v,10));
        print("2 log_mu  ", precision(log_mu,10));
        print("2 log_mu_p  ", precision(log_mu_p,10));
        if(DEBUG_MEMBERSHIP,
        print("\n Complex membership check: \n",
            "- Checking if v is in unit lattice, v: ", precision(v,10),
            "\n- psi(mu)=log_mu: ", precision(log_mu,10), \
            "\n- Norm(v - psi(mu)) = ", precision(norml2(log_mu_p),10) );
        );

        if(abs(norml2(real(log_mu_p)) ) <= eps, print("return element close to 0"); breakpoint(); return(1) );
        elmp = exponentiated_embedding(G, log_mu_p);
        print("exp embed ", precision(exponentiated_embedding(G, log_mu_p),10));
        coeffs = M_inv*elmp~;

        \\print(precision(coeffs,10));

        \\ check if the vector is all integers
        for(i=1, length(coeffs),
            if( abs(imag(coeffs[i]) ) >= eps, return(0));
            if( abs(real(coeffs[i])-round(coeffs[i])) >= eps, return(0););
        );
        if(DEBUG_MEMBERSHIP,print("coeffs are all integers"););
        \\
        \\ Also need to check that the ideal 1/mu' I = O_K
        check_ideal = idealdiv(G, new_y, round(coeffs));
        if(DEBUG_MEMBERSHIP, print("Verifying (1/mu')I = O_K"),check_ideal == matid(n));
        if (check_ideal == matid(n),
            print("Ideal is identity"); return(1),
            return(0));

    );
}


\\ simple function which given a list of cpct reps, computes the lcm of all of the denominators
\\ INPUT: A list of compact representations
\\ OUTPUT: LCM of the denominators of those compact representations
lcm_denominators(cpct_reps)={
    least = 1;
    for(i=1, length(cpct_reps),
        least = lcm(least, lcm(cpct_reps[i][2]));
    );
    return(least);
}


\\ INPUT:
\\ - G an number field,
\\ - p a prime number,
\\ - k the index of the unit eta_k (See Pohst or Arenz). k =0 implies we are adjusting with torsion
\\ - cpctreps is the cpct reps of the current original set of independent units (no updates)
\\ - expmat indicates how to get the updated units {eta_i}
\\ - logmat corresponds to the logarithm embeddings of the elements of cpctreps
\\ OUTPUT:
\\ - An ideal Q as illustrated in Pohst/Arenz
\\ - along with a new compact representation (if needed)
get_search_ideal_cpct(G,p, k, cpctreps, expmat, logmat, torsion_coeffs=[])={
    my(check_poly,
      field_deg = length(G.zk),               \\ degree of the number field
      prime_q,                                \\ holder for the pohst prime
      prime_decomposition,                    \\ holds ideal factorization of the prime p
      factor_candidate,                       \\ variable for holding a prime ideal
      res_field,                              \\ holds the residue field generated by factor_candidate
      embedded_eta,                           \\ used to check irreducibility mod the factor_candidate
      expvec,
      ideal_candidate_found = 0,              \\ flag to indicate an ideal has been found
      lowbound =1,                            \\ the lower bound for get_pohst_prime
      denom_lcm = lcm_denominators(cpctreps),
      prime_ctr = 1
    );

    while(!ideal_candidate_found,
        prime_q = get_pohst_prime2(G.disc, p, lowbound, denom_lcm);
        prime_decomposition = idealprimedec(G,prime_q);                         \\ Decompose prime_q as ideals in G (as vector of prid structures)

        if(DEBUG_IDEAL_SEARCH,print(" - Found Pohst prime q: ", prime_q, "\n - Number of Factors of q: ", length(prime_decomposition)););
        \\ loop over all the prime ideals, look for one where the polynomial t^p - pohst_eta = 0 is irreducible
        for(i = 1, length(prime_decomposition),

            factor_candidate = prime_decomposition[i];                          \\ holds the ith ideal of the decomposition, call it P

            res_field = nfmodprinit(G,factor_candidate);                        \\ converts prime ideal P into residue field OK/P
            cpctreps = cpct_denom_switch(G, logmat, cpctreps, eps, prime_q);    \\ ensure the component representations are coprime modulo prime_q

            \\ handles the torsion case
            if(k == 0,
                embedded_eta = nfmodpr(G, nfalgtobasis(G,nfrootsof1(G)[2]),res_field);
            ,\\else
                expvec = expmat[,k];
                embedded_eta = cpct_modp_from_factors(G, cpctreps, expvec, res_field);

                if(length(torsion_coeffs) != 0,
                    embedded_eta *= nfmodpr(G, nfalgtobasis(G,nfrootsof1(G)[2]),res_field)^torsion_coeffs[k];
                );
            );


            \\ irreducibility check of Pohst/Arenz
            if (embedded_eta^( (factor_candidate[1]^factor_candidate.f-1)/p) != 1,
                \\write("findingQ.log", "p = ",p, ". Ideals checked = ",prime_ctr, " Q = ", (prime_q-1)/p , "*", p, "+1. Ideal norm ", idealnorm(G,factor_candidate)  );
                return([factor_candidate,cpctreps]);                                           \\ Should always exit the function here unless the input is invalid
            );
            prime_ctr+=1;
        );
        lowbound = prime_q++;
    );
    print( "ERROR: Did not find any prime ideal with the irreducibility property.");
    return(0);
}; \\ end function get_search ideal


\\ INPUT:
\\ - G an number field,
\\ - p a prime number,
\\ - k the index of the modifying element (I don't think we actually use this)
\\ - eta (See Pohst or Arenz, this is eta_k) as a polynomial
\\ - the optional flag cpct_input should be set to 1 if the elements are provided as compact representations
\\ OUTPUT:
\\ - outputs the new_eta = (eta_k^t)*eta_j, with new_eta a pth power in the residue field of searchideal.
get_new_eta_cpct(G, p, k, j, searchideal, cpctreps, expmat,torsion_coeffs=[])={

    my(
      res_field, resfield_gen,
      eta_k_mod, eta_j_mod,
      dlog_k, dlog_j, new_exp
    );

    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: The size of the res.field multiplicative group: ", searchideal.p^searchideal.f););

    res_field = nfmodprinit(G, searchideal);                                    \\ initialize the residue field
    if(k==0,
        eta_k_mod = nfmodpr(G, nfalgtobasis(G, nfrootsof1(G)[2]),res_field);
    ,
        eta_k_mod = cpct_modp_from_factors(G, cpctreps, expmat[,k], res_field);
    );
    eta_j_mod = cpct_modp_from_factors(G, cpctreps, expmat[,j], res_field);

    if(length(torsion_coeffs)!=0,
        embedded_torsion = nfmodpr(G, nfalgtobasis(G,nfrootsof1(G)[2]),res_field);
        eta_k_mod *= (embedded_torsion)^torsion_coeffs[k];
        eta_j_mod *= (embedded_torsion)^torsion_coeffs[j];
    );
    resfield_gen = ffprimroot(ffgen(eta_k_mod));                                \\ gets a multiplicative generator

    if(DEBUG_NEW_ETA, print("eta_k_j ",eta_k_mod ););

    dlog_k = fflog(eta_k_mod, resfield_gen);                                    \\ get the DLOG wrt the mult. generator
    dlog_j = fflog(eta_j_mod, resfield_gen);

    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: DLOGS of eta_k, eta_j  ", dlog_k, "  ,   ", dlog_j ););

    new_exp = lift(-dlog_j*Mod(dlog_k,p)^(-1));
    if( DEBUG_NEW_ETA, print(" ---- get_new_eta: Computed eta_",k," power: ",new_exp ););

    return(new_exp);
};

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ This is a subroutine for the pohst search function.
\\ INPUTS:
\\    G - a number field
\\    p - a prime
\\    k is the index of the 'first' element, which the later elements are modified by
\\    cpct_reps - list of the compact representations corresponding to loglat
\\    expmat - indicates how to compute the logs of the intermediate eta values from loglat
\\    loglat - logarithm matrix corresponding to current system of units (not the eta values)
\\ OUTPUT:
\\ - The new set of eta, which keeps the first k elements the same, and the rest are modified.
update_eta_set_log(G,p,k, cpct_reps, expmat, loglat, torsion_coeffs = [])={
    my(eta_k, eta_i, cebotarev_ideal, exponent_tracker=[], exponent_i, temporary_cpctreps);

    [cebotarev_ideal,temporary_cpctreps] = get_search_ideal_cpct(G,p, k, cpct_reps, expmat, loglat, torsion_coeffs);

    for(i=k+1, length(cpct_reps),
        exponent_i = get_new_eta_cpct(G, p , k, i, cebotarev_ideal, temporary_cpctreps, expmat,torsion_coeffs);
        exponent_tracker = concat(exponent_tracker, exponent_i);
    );
    return(exponent_tracker);
}


torsion_update(~torsionCoeffs, ~updateVector, k)=
{
    for(s = k+1, length(torsionCoeffs),
        torsionCoeffs[s] = (torsionCoeffs[s]+updateVector[s-k]*torsionCoeffs[k]);
    );
    return(torsionCoeffs);
}


log_pohst_pari(G,L,unitvector_cpct, B, eps, OUTFILE1 = "log_pohst_output")={
    print("LPohst: start. p = ", 2);

    my(new_units = L, index_holder = 1, index_bound = B, solution, solutionflag = 0, p = 2);
    my(eta_exp_mat, compact_lcm, test_eta_k, updatevec);
    my( torsion_coeffs, [torsion, torsiongen] = nfrootsof1(G) );
    betavec = unitvector_cpct;

    compact_lcm = lcm_denominators(unitvector_cpct);
    my(time_pthRoot = 0, t_pthRootBefore, t_pthRootAfter,
        time_found = 0, t_foundBefore, t_foundAfter);
    initialTime = getabstime();
    while(p <= index_bound,

        \\ include torsion in pari check when p is not coprime
        if(torsion %p == 0,
            betavec = concat(betavec, [[ List( [nfalgtobasis(G, torsiongen)] ), [1]  ]] );
        );

        lpohst_ploop_t1 = 0;
        \\if(1,                                                                 \\ UNCOMMENT LINE AND COMMENT OUT THE LINE BELOW IF YOU WANT TO SKIP THE FAST PRIME CHECK
        if(pari_prime_check(G, betavec, p, compact_lcm, 1) == 0 ,

            print("\npari_prime_check detects possible index divisor ", p);
            eta_exp_mat = matid(length(new_units));

            torsion_coeffs = [];
            if(gcd(p,torsion)!=1,
                \\ update eta set using the torsion unit, store the coeffs in a vector
                torsion_coeffs = update_eta_set_log(G,p,0,unitvector_cpct, eta_exp_mat, new_units);
            );

            k = 1;                                                              \\# each time a new prime is considered, reset k to 1
            solutionflag = 1;                                                   \\# flag indicates if a solution has been found
            solution = 0;                                                       \\# variable to hold either 0 or a found solutions

            while(solutionflag == 1,

                t_pthRootBefore = getabstime();
                \\\# step 1 of Algorithm 8, the pth root test
                test_eta_k = new_units*eta_exp_mat[,k];
                test_eta_k = test_eta_k/p;

                /*# commented out because method can't be proven
                eta_k_complex_log = vector(G.r1+G.r2, t, 0);
                for(i =1, length(unitvector_cpct),
                    eta_k_complex_log += eta_exp_mat[i,k]*complex_log_from_cpct(G, unitvector_cpct[i]);
                );

                eta_k_complex_log = eta_k_complex_log/p;
                GP_ASSERT_VEC_NEAR(real(eta_k_complex_log)[1..G.r1+G.r2-1]~,test_eta_k,eps);
                print("complex eta_k ",precision(eta_k_complex_log,10));
                */
                lattice_check_t1 = getabstime();

                \\\# here, test_eta_k is known on all coordinates because
                \\\# it is a lin comb. of log vectors of units, so the sum of coordinates is just 0
                solutionflag = check_in_unitlattice(G, test_eta_k~, eps);

                \\solutionflag1=complex_check_in_unitlattice(G, eta_k_complex_log, eps); # commented out because method can't be proven
                \\print(solutionflag1, "  ", solutionflag);
                t_pthRootAfter = getabstime();
                time_pthRoot+= (t_pthRootAfter - t_pthRootBefore);

                \\# commented out because method can't be proven
                \\if (solutionflag != complex_check_in_unitlattice(G, eta_k_complex_log/p, eps), print("flags not matching");breakpoint());

                \\print("DEBUGGING: lattice check time: ", getabstime() - lattice_check_t1);


                if(solutionflag ==1,
                    print("LPohst: Found Solutions --");
                    t_foundBefore = getabstime();
                    \\# These two lines use mlll instead of column replacement
                    \\new_units = new_units*eta_exp_mat;
                    \\new_units = my_mlll(matconcat([new_units, test_eta_k]),eps);

                    \\ These lines replace the kth unit, which is different than described in Pohst
                    new_units = replace_column(new_units, k, test_eta_k);
                    new_units = new_units*qflll(new_units);

                    cpct_t1 = getabstime();
                    unitvector_cpct = cpct_from_loglattice(G, new_units, eps);
                    cpct_t2 = getabstime();
                    print("recomputing compact reps time: ", (cpct_t2 -cpct_t1 ) );
                    \\for(i=1, length(unitvector_cpct),
                    \\    v1 = log(nfeltembed(G,compact_reconstruct(G, unitvector_cpct[i][1],unitvector_cpct[i][2]) ));
                    \\    v2 = complex_log_from_cpct(G, unitvector_cpct[i]) ;
                    \\    print(precision(norml2(v1-v2),10));
                    \\);

                    \\ if the cpct reps have changed, then need to update betavec and the lcm
                    betavec = unitvector_cpct;
                    compact_lcm = lcm_denominators(unitvector_cpct);

                    eta_exp_mat = matid(length(new_units));
                    index_bound = ceil(index_bound/p);
                    index_holder = index_holder*p;
                    if(gcd(p,torsion)!=1,
                        \\ update eta set using the torsion unit, store the coeffs in a vector
                        torsion_coeffs = update_eta_set_log(G,p,0,unitvector_cpct, eta_exp_mat, new_units);
                    );
                    k = 1;
                    t_foundAfter = getabstime();
                    time_found += (t_foundAfter-t_foundBefore);
                    if (index_bound ==1, print("Index is now 1. Ending LPohst"); return(new_units));
                , \\else                                                        # this is the case when no sol is found
                    if(k == length(L),solutionflag = 0; break;);

                    updatevec = update_eta_set_log(G,p,k,unitvector_cpct, eta_exp_mat, new_units, torsion_coeffs);
                    eta_exp_mat = update_expmat(eta_exp_mat, updatevec, k , p );

                    torsion_coeffs = torsion_update(~torsion_coeffs, ~updatevec, k);

                    \\ Update k and reset solutionflag to 1

                    k+=1;
                    solutionflag = 1;
                );
            );

        , \\else:
            if(p%100000 ==1, print("LPohst: pari_check succeeds for, ",p, ". Now checking p = ", nextprime(p+1)));
        );
        if(torsion %p == 0,
            betavec = unitvector_cpct;
            compact_lcm = lcm_denominators(unitvector_cpct);                    \\ used as the 'bad' input to pari_prime_check, ignores non-coprime primes during the equation finding step
        );
        \\if(p>500,
        \\    print("PthRoot: ", time_pthRoot, "\nUpdate: ", time_update ); breakpoint());

        if (getabstime()-initialTime > (12*60*60*1000), write(OUTFILE1, G.pol, " time exceeds 12 hours");break;);

        p = nextprime(p+1);
    );
    return(new_units);

}
