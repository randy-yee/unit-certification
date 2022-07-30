read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
FIRST = 1;
\\ This is the generalized version of Ha's scanball algorithm
\\ subroutine of Schoof's scan algorithm. It checks a small ball B_D around the divisor D = (y,u) and returns the minima within
\\ INPUT:
\\ - y coefficient matrix of an ideal
\\ - u a positive real vector of length r= r1 +r2 -1
\\ - G the number field
\\ - psimu is a log vector used to keep track of the distance in the log lattice.
\\ - web is the max distance between "web" points
\\ - eps is the usual error
\\ OUTPUT:
\\ - Adds minima to the Map object bmap
/******************************************************************************/
scanball_map(~G, ~bmap, y, u, psimu, web, eps, ~repeated_minima)={

    my(
        n = poldegree(G.pol),
        x, scan_bound,
        vecholder, gram_mat,
        scan_elements,
        LLL_reduced_yu
    );
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\\ #Use the y,u to define a lattice to scan for elements
    x = G[5][1]*y;                                                              \\ numerical representation of y (complex)
    x = mulvec(x,u);                                                            \\ compute y*u
    x = embed_real(G,x);
    LLL_reduced_yu = x*qflll(x);                                                \\ lll reduce y*u
    vecholder = LLL_reduced_yu[,1];                                             \\ short vector, 1st element of LLL basis
    scan_bound = sqrt(n)*exp(2*web)*sqrt(norml2(vecholder));                    \\ See schoof alg 10.7, e^(2*var_eps)*sqrt(n)*sqrt(norml2(col))
    gram_mat=LLL_reduced_yu~*LLL_reduced_yu;                                    \\ get the gram matrix

    scan_elements = qfminim(gram_mat,scan_bound^2,,2)[3];
    scan_elements = y*scan_elements;                                            \\ get scanned elements wrt integral basis
    my(
        norm_deltaK = ceil(((2/Pi)^(G.r2))*abs(G.disc)^(1/2)*idealnorm(G,y)),
        eltnorm = 0,
        new_y,
        real_y,
        new_yLLL,
        psi_value,
        vec_numerical
    );
    for(ii=1, length(scan_elements),

        \\\ #Easy necessary condition for minimum'''
        \\\ #norm of a minimum should satisfy 1 < N(s_elt) < N(y)*delta_K
        eltnorm = abs(nfeltnorm(G,scan_elements[,ii] ));
        if(eltnorm>=1 && eltnorm<=norm_deltaK,

            \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            \\\# check if ideal is reduced, which implies we have a minimum
            \\\# note that the norm condition is just a simple
            new_y = idealdiv(G,y,scan_elements[,ii]);
            real_y = embed_real(G, G[5][1]*new_y);                              \\ get the ideal (1/w_i)*y

            if (FIRST, print("WARNING: confirm this ideal norm check is accurate"); FIRST =0);
            new_yLLL = real_y*qflll(real_y);
            if(norml2(new_yLLL) > 1-eps,
                if(checkred_old(new_y,G,eps)==1,
                    vec_numerical = (G[5][1]*scan_elements[,ii])~;
                    psi_value = vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps);
                    if(mapisdefined(bmap, new_y, &existing_entry),
                        repeatflag = is_repeat_babystock(existing_entry, psi_value, eps);
                        if(repeatflag==0,
                            listput( ~existing_entry, psi_value);
                            mapput(bmap, new_y, existing_entry );
                        , \\else
                            repeated_minima+=1;
                        );
                    ,\\else
                        mapput(bmap, new_y, List([psi_value]));
                    );
                );
            );
        );
    );
}

\\ distingished from scanball_map as instead of one psimu, a list of them
\\ is provided. In this way, when the ideal y is repeated, we can reduce overall
\\ number of scans
overlap_scanball(~G, ~bmap, ~y, ~u, ~log_distance_list, web, eps, ~repeated_minima)={

    my(
        n = poldegree(G.pol),
        x, scan_bound,
        vecholder, gram_mat,
        scan_elements,
        LLL_reduced_yu
    );
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\\ #Use the y,u to define a lattice to scan for elements
    x = G[5][1]*y;                                                              \\ numerical representation of y (complex)
    x = mulvec(x,u);                                                            \\ compute y*u
    x = embed_real(G,x);
    LLL_reduced_yu = x*qflll(x);                                                \\ lll reduce y*u
    vecholder = LLL_reduced_yu[,1];                                             \\ short vector, 1st element of LLL basis
    scan_bound = sqrt(n)*exp(2*web)*sqrt(norml2(vecholder));                    \\ See schoof alg 10.7, e^(2*var_eps)*sqrt(n)*sqrt(norml2(col))
    gram_mat=LLL_reduced_yu~*LLL_reduced_yu;                                    \\ get the gram matrix

    scan_elements = qfminim(gram_mat,scan_bound^2,,2)[3];
    scan_elements = y*scan_elements;                                            \\ get scanned elements wrt integral basis
    my(
        norm_deltaK = ceil(((2/Pi)^(G.r2))*abs(G.disc)^(1/2)*idealnorm(G,y)),
        eltnorm = 0,
        new_y,
        real_y,
        new_yLLL,
        psi_value,
        vec_numerical
    );
    print("elements scanned: ", length(scan_elements));
    for(ii=1, length(scan_elements),

        \\\ #Easy necessary condition for minimum'''
        \\\ #norm of a minimum should satisfy 1 < N(s_elt) < N(y)*delta_K
        eltnorm = abs(nfeltnorm(G,scan_elements[,ii] ));
        if(eltnorm>=1 && eltnorm<=norm_deltaK,

            \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            \\\# check if ideal is reduced, which implies we have a minimum
            \\\# note that the norm condition is just a simple
            new_y = idealdiv(G,y,scan_elements[,ii]);
            real_y = embed_real(G, G[5][1]*new_y);                              \\ get the ideal (1/w_i)*y

            if (FIRST, print("WARNING: confirm this ideal norm check is accurate"); FIRST =0);
            new_yLLL = real_y*qflll(real_y);
            if(norml2(new_yLLL) > 1-eps,
                if(checkred_old(new_y,G,eps)==1,
                    vec_numerical = (G[5][1]*scan_elements[,ii])~;
                    for(j = 1, length(log_distance_list),
                        psi_value = vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+log_distance_list[j][1..G.r1+G.r2-1],eps);
                        if(mapisdefined(bmap, new_y, &existing_entry),
                            repeatflag = is_repeat_babystock(existing_entry, psi_value, eps);
                            if(repeatflag==0,
                                listput( ~existing_entry, psi_value);
                                mapput(bmap, new_y, existing_entry );
                            , \\else
                                repeated_minima+=1;
                            );
                        ,\\else
                            mapput(bmap, new_y, List([psi_value]));
                        );
                    );
                );
            );
        );
    );
}
/******************************************************************************/
/*28. Find a new vector v \in \Lambda_K\Lambda' where (O_K, v) \in \ep_B: norml2, fun. in 13 (is_vec_in_lattice),  fun. in 11 (updatelamda)*/
\\ Check the set ball_minima for a new vector in Lambda not already in Lambda'
\\ INPUT:
\\ - ball_minima a list of minima, each of the form [ideal, log vector] = [1/mu*O_K, vector = log(mu)]
\\ - L is the lattice we are checking (Lambda')
\\ - G the number field
\\ - eps the error

/******************************************************************************/
check_units_bstock(~bmap,~L,~G,eps)={

    my(
        ideal_identity = matid(n),
        n:small = poldegree(G.pol),
        new_counter:small,
        candidate, eps_sqr = eps^2,
        minlist
    );
    new_counter = 0;
    if (mapisdefined(bmap, ideal_identity, &minlist),
        for(i=1,length(minlist),                                       \\ if the ideal (1/mu)*O_k = O_k then check further
            candidate=minlist[i][2];                                                \\ candidate = log(mu)
            if(norml2(candidate)>eps_sqr&&is_vec_in_lattice(candidate~,L,eps_sqr)==0,           \\ if nonzero and v is not already in L, then
                new_counter+=1;
                print("Babystock unit found, " precision(L,10), "  ", precision(candidate,10));
                L = my_mlll(matconcat([L, candidate~]), eps);
            )
        );

    );
    return([L,new_counter]);        \\ modifed to return new_counter as well
}
