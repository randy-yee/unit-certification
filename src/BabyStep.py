read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
FIRST = 1;
VERIFY_GENERATORS = 1;
adjusttime1 = 0;
localbitprecision = 30;

\\ INPUT:
\\ - G a number field
\\ - giant_divisor is a 3-element list consisting of
\\      an ideal, a vector of size G.r1+G.r1 corresponding to valuation of an element
\\      a compact representation
\\ - expected position is the coordinates in R^r of your giant lattice element
\\ - distance_ok is the limit of the acceptable distance from the divisor to the target position
\\ - tracker is a list of compact representation corresponding to the current element
\\ - trackerLog is the logarithm value
superctr = 0;
supertimer = 0;
adjust_giant_step_cpct(~G, ~giant_divisor, ~tracker, trackerLog, ~expected_position, s_radius, eps, storage = "LOG")={
    GP_ASSERT_TRUE(eps > 0);
    my(
        r = G.r1 + G.r2 -1,
        divisor_distance,
        localbitprecision = 30  \\ only used to get a rough gauge on distance
    );
    a_time1 = getabstime();

    if (storage == "COMPACT",
        trackerLog += trackerLogarithm(G, ~tracker, r);
    );
    superctr +=1;
    supertimer += (getabstime()-a_time1 );
    \\if(superctr %1000 == 0, print("time spent in adjust: ",supertimer));
    mainbitprecision = default(realbitprecision);
    default(realbitprecision, localbitprecision);
    divisor_distance = expected_position - trackerLog~;
    default(realbitprecision, mainbitprecision);
    divisor_distance = bitprecision(divisor_distance, mainbitprecision);

    a_time2 = getabstime();
    adjusttime1 += (a_time2 - a_time1);

    my(newFactor, logNewFactor, adjustment_divisor, reduced_product, new_distance,
                new_divisor);
    if(sqrt(norml2(divisor_distance)) < s_radius,
        return([giant_divisor, trackerLog]);
    , \\else
        for(i=1, 2,

            adjustment_divisor = get_nearby_rdivisor(G, matid(poldegree(G.pol)), divisor_distance, i%2);
            if (sqrt(norml2(adjustment_divisor[3])) < eps,
                return(giant_divisor);
            ,
                new_divisor = [idealmul(G, giant_divisor[1], adjustment_divisor[1]),
                    pointwise_vector_mul(giant_divisor[2],adjustment_divisor[2] )~ ];

                reduced_product = reddiv_compact(new_divisor[1], new_divisor[2],G, G[5][1] );
                newFactor = nfeltmul(G, adjustment_divisor[4], reduced_product[4]);

                default(realbitprecision, localbitprecision);
                logNewFactor = log(abs(nfeltembed(G,newFactor) ));
                new_distance = norml2(expected_position - logNewFactor~);
                default(realbitprecision, mainbitprecision);
                new_distance = bitprecision(new_distance, mainbitprecision);

                if(new_distance < norml2(divisor_distance)+eps,
                    trackerLog += logNewFactor;
                    listput(~tracker, [newFactor, 1]);
                    \\print("baby verify: adjustment return");verify_generator_with_list(G, reduced_product[1], tracker);
                    return([reduced_product, trackerLog]);
                );
            );
        );
        \\ adjustment fails, so we should use jump to compute something close
        gstep_divisor = jump_compact(matid(length(G.zk)), expected_position, G, length(G.zk), eps);
        trackerLog = log_from_cpct(G, gstep_divisor[3]);
        print("WARNING: recomputing current divisor with jump.");
        return([gstep_divisor, trackerLog]);
    );
};

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
                    psi_value = log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu;
                    if(mapisdefined(bmap, new_y, &existing_entry),
                        repeatflag = is_repeat_babystock(existing_entry, psi_value, eps);
                        if(repeatflag==0,
                            listput( ~existing_entry, psi_value);
                            mapput(~bmap, new_y, existing_entry );
                        , \\else
                            repeated_minima+=1;
                        );
                    ,\\else
                        mapput(~bmap, new_y, List([psi_value]));
                    );
                );
            );
        );
    );
}





\\ distingished from scanball_map as instead of one psimu, a list of them
\\ is provided. In this way, when the ideal y is repeated, we can reduce overall
\\ number of scans
\\ bmap is passed by reference, and any new minima are added to it
overlap_scanball(~G, ~bmap, ~y, ~u, ~log_distance_list, ball_distance, eps, ~repeated_minima)={
    SCAN_STATS = 0;
    if(SCAN_STATS, print("ball dist, eps", precision(ball_distance,10),"  " eps););
    my(
        n = poldegree(G.pol),
        x, scan_bound,
        vecholder, gram_mat,
        scan_elements,
        LLL_reduced_yu
    );
    scanstart = getabstime();
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\\ #Use y,u to define a lattice to scan for elements
    x = G[5][1]*y;                                                              \\ numerical representation of y (complex)
    x = mulvec(x,u);                                                            \\ compute y*u
    x = embed_real(G,x);
    LLL_reduced_yu = x*qflll(x);                                                \\ lll reduce y*u
    vecholder = LLL_reduced_yu[,1];                                             \\ short vector, 1st element of LLL basis
    \\# ensure that this is actually the shortest vector
    \\# or just determine the length of the shortest vector
    scan_bound = sqrt(n)*exp(2*ball_distance)*sqrt(norml2(vecholder));          \\ See schoof alg 10.7, e^(2*var_eps)*sqrt(n)*sqrt(norml2(col))
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
    if(SCAN_STATS, print("enumerated elts: ",  length(scan_elements), "  scan bound: ", precision(scan_bound,10)););
    total_check_reduced_time = 0;
    total_add_element_time = 0;
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
            new_yLLL = real_y*qflll(real_y);
            if(1,
            /*
                if (check_ideal_reduced(G, new_y) != checkred_old(new_y,G,eps),
                    print("ideal check mismatch:");
                    print(G.pol, "  ", new_y, "  ", eps);

                );
                \\print("comparing reduced ideal checks. delete when resolved");
                \\if(checkred_old(new_y,G,eps)==1,
            */
                checkstart = getabstime();
                if(check_ideal_reduced(G, new_y),
                    check_end = getabstime();
                    total_check_reduced_time += (check_end -checkstart);
                    vec_numerical = (G[5][1]*scan_elements[,ii])~;
                    \\ # Start at 2 because j=1 holds the nu value ( see babystock_scan_jump )
                    add_start = getabstime();

                    for(j = 2, length(log_distance_list),
                        \\print("WARNING: Determine what the right level of precision to round these is");
                        psi_value = vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2]))+log_distance_list[j][1..G.r1+G.r2],eps^2);
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
                    add_end = getabstime();
                    total_add_element_time+=(add_end-add_start);
                );
            );
        );

    );
    scan_end = getabstime();
    if(SCAN_STATS, print("done enumerating", "  ", total_check_reduced_time, "  ", total_add_element_time, " Total: ", (scan_end-scanstart)););
}


\\ distingished from scanball_map as instead of one psimu, a list of them
\\ is provided. In this way, when the ideal y is repeated, we can reduce overall
\\ number of scans
\\ bmap is passed by reference, and any new minima are added to it
compact_storage_overlap_scanball(~G, ~bmap, ~y, ~u, ~log_distance_list, ball_distance, eps, ~repeated_minima, cpct_list, ~cpct_bstock )={
    my(
        n = poldegree(G.pol),
        x, scan_bound,
        vecholder, gram_mat,
        scan_elements,
        LLL_reduced_yu
    );
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\\ #Use y,u to define a lattice to scan for elements
    x = G[5][1]*y;                                                              \\ numerical representation of y (complex)
    x = mulvec(x,u);                                                            \\ compute y*u
    x = embed_real(G,x);
    LLL_reduced_yu = x*qflll(x);                                                \\ lll reduce y*u
    vecholder = LLL_reduced_yu[,1];                                             \\ short vector, 1st element of LLL basis
    \\# ensure that this is actually the shortest vector
    \\# or just determine the length of the shortest vector
    scan_bound = sqrt(n)*exp(2*ball_distance)*sqrt(norml2(vecholder));          \\ See schoof alg 10.7, e^(2*var_eps)*sqrt(n)*sqrt(norml2(col))
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
            if(1,
                if(checkred_old(new_y,G,eps)==1,
                    vec_numerical = (G[5][1]*scan_elements[,ii])~;
                    \\ # Start at 2 because j=1 holds the nu value ( see babystock_scan_jump )
                    for(j = 2, length(log_distance_list),
                        psi_value = log(abs(vec_numerical[1..G.r1+G.r2]))+log_distance_list[j][1..G.r1+G.r2];

                        \\print("psi (low  prec): ", precision(psi_value,20), "  ");
                        cpctList = cpct_list[j-1];
                        listput(~cpctList[2], [scan_elements[,ii], 1]);
                        \\print("psi (high prec): ", precision(log_from_cpct(G, cpctList[1])+trackerLogarithm(G, cpctList[2], G.r1+G.r2-1)[1..G.r1+G.r2],20));
                        if(mapisdefined(bmap, new_y, &existing_entry),
                            repeatflag = is_repeat_babystock(existing_entry, psi_value, eps);
                            if(repeatflag==0,
                                listput( ~existing_entry, psi_value);
                                mapput(~bmap, new_y, existing_entry );
                                mapput(~cpct_bstock, new_y, [cpctList]);
                            , \\else
                                repeated_minima+=1;
                            );
                        ,\\else
                            mapput(~bmap, new_y, List([psi_value]));
                            mapput(~cpct_bstock, new_y, List([cpctList]));
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
check_units_bstock(~bmap, ~L, ~G, eps)={
    my(
        n:small = poldegree(G.pol),
        ideal_identity = matid(n),
        new_counter:small,
        candidate, eps_sqr = eps^2,
        minlist
    );
    new_counter = 0;
    if (mapisdefined(bmap, ideal_identity, &minlist),
        for(i=1,length(minlist),                                                               \\ if the ideal (1/mu)*O_k = O_k then check further
            candidate=minlist[i];                                                              \\ candidate = log(mu)
            if(norml2(candidate)>eps_sqr&&is_vec_in_lattice(candidate[1..r]~,L,eps_sqr)==0,          \\ if nonzero and v is not already in L, then
                new_counter+=1;
                print("Babystock unit found, " precision(L,10), "  ", precision(candidate,10));
                L = my_mlll(matconcat([L, candidate~]), eps);
            )
        );

    );
    return([L,new_counter]);        \\ modifed to return new_counter as well
}

\\ check babystock for units, when it is stored as compact rep lists
check_compact_bstock(~cpct_bstock, ~L, ~G, eps)={

    my(
        n:small = poldegree(G.pol),
        ideal_identity = matid(n),
        new_counter:small,
        candidate, eps_sqr = eps^2,
        minlist,
        r = G.r1+G.r2 -1
    );
    new_counter = 0;
    if (mapisdefined(cpct_bstock, ideal_identity, &minlist),
        for(i=1,length(minlist),                                                               \\ if the ideal (1/mu)*O_k = O_k then check further
            candidate=log_from_cpct(G, minlist[i][1])+trackerLogarithm(G, minlist[i][2], r);                                                              \\ candidate = log(mu)
            if(norml2(candidate)>eps_sqr&&is_vec_in_lattice(candidate[1..r]~,L,eps_sqr)==0,          \\ if nonzero and v is not already in L, then
                new_counter+=1;
                print("Babystock unit found, " precision(L,10), "  ", precision(candidate,10));
                L = my_mlll(matconcat([L, candidate~]), eps);
            )
        );

    );
    return([L,new_counter]);        \\ modifed to return new_counter as well
}

babystockPrecision(G, sublatticeBasis)={
    my(X2,sumv = sublatticeBasis[,1]);
    for(j=2, length(sublatticeBasis), sumv+=sublatticeBasis[,j]);
    X2 = prec_baby(poldegree(G.pol), log(abs(G.disc)), normlp(sumv));
    return(ceil(X2));
}

initialize_babystock_edges(~giant_legs, scan_shrink_factor, r)=
{
    my(denoms, babystock_t);
    \\\# Note, these denoms are chosen based on the length r vector, not the r+1 vector
    denoms = vector( r, i, sqrt(norml2(giant_legs[,i][1..r])));

    print("norm of giant legs: ", precision(denoms,10), "  ", precision(scan_shrink_factor,10));
    denoms/= scan_shrink_factor;
    print("denominators after scaling: ", precision(denoms,10));

    for(i=1, length(denoms),
        if (denoms[i]-floor(denoms[i]) < 0.00001,
            denoms[i] = floor(denoms[i] ),
            denoms[i] = ceil(denoms[i]);
        );
        if (denoms[i] ==0, denoms[i] = 1);
    );
    babystock_t = giant_legs;
    for(i=1, length(babystock_t), babystock_t[,i]=babystock_t[,i]/denoms[i];);
    return([babystock_t, denoms]);
}

\\ #Subalgorithm of bsgs. computes baby steps incrementally using ideal
\\ #multiplication to obtain adjacent elements.
\\ #Compact Represenation Version
incremental_baby_steps(y, ~lattice_lambda, ~giant_legs,\
                        ~baby_hashmap, G, scan_radius, eps, outFileInfo=[]) =
{
    print("baby steps using increments (Log version)");
    my(timeout, OUTFILE_BS);
    if(length(outFileInfo) == 2,
        timeout = outFileInfo[2];
        OUTFILE_BS = outFileInfo[1];
    ,
        timeout = 0;
        OUTFILE_BS = 0;
    );

    my(
        field_deg = poldegree(G.pol),   r = G.r1+G.r2-1,
        zero_vec = vector(r, i, 0),     identity = matid(field_deg),
        web_coords = zero_vec,          place_marker = r,
        directions = vector(r, i, 1),
        denoms,
        scan_shrink_factor,
        babystock_t,
        expected_position = vector(r+1, i, 0)~,
        start_time
    );

    GP_ASSERT_TRUE(eps > 0);
    GP_ASSERT_EQ(r, length(giant_legs));

    REQ_BS = babystockPrecision(G, giant_legs);
    print("prec: ", default(realbitprecision), "  ", REQ_BS);
    default(realbitprecision, REQ_BS);


    start_time = getabstime();


    [babystock_t, denoms] = initialize_babystock_edges(~giant_legs, scan_radius, r);
    increment_coordinates(denoms, web_coords);

    my(
        direction_elements,
        inverse_direction_elements,
        scanIdeals = Map(),
        distanceList = List(),
        scanIdealsMatrix,
        repeat_counter = 0,
        repeated_minima = 0,
        logdist,
        compactTracking = List(),
        trackingLog = vector(r+1, i, 0),
        idealCompactGenerator = Map(),
        base_value = [[1], [1]],
        s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    );

    \\ #vectors of length r which are used to compute an adjacent element in
    \\ #a particular direction in log space
    \\ #direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, babystock_t, field_deg, eps);

    \\ # stores the current element as a list of compact representations.
    \\ # the first r are the precomputed elements along with an integer indicating a current multiple
    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );
    emptyCompactList = compactTracking;
    \\ #initalization of the current element
    baby_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Additional timing variables and file writes");
    my(baby_t1, baby_tn, baby_tmid, ctr);
    baby_t1 = getabstime();
    baby_tmid = baby_t1;

    \\# increments are done modulo the product of avec, returns to zero when
    \\# all elements have been visited
    while(web_coords != zero_vec,
        \\print(web_coords, "  BD: ", precision(trackingLog,10));
        if(directions[place_marker] == 1,
            baby_divisor = [idealmul(G, baby_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],direction_elements[place_marker][2] )~];
            expected_position += babystock_t[,place_marker];
            compactTracking[place_marker][2] += 1;
            trackingLog += log_from_cpct(G, direction_elements[place_marker][3]);
        ,\\else
            baby_divisor = [idealmul(G, baby_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],inverse_direction_elements[place_marker][2] )~];
            expected_position -= babystock_t[,place_marker];
            compactTracking[place_marker][2] -= 1;
            trackingLog -= log_from_cpct(G, direction_elements[place_marker][3]);
        );

        wasUpdated = get_next_giant_divisor_cpct(G, ~baby_divisor, ~compactTracking);
        if (wasUpdated,
            trackingLog += log(abs(nfeltembed(G, compactTracking[length(compactTracking)][1])));
        );
        [baby_divisor, tempLog] = adjust_giant_step_cpct(~G, ~baby_divisor,~compactTracking, ~trackingLog, ~expected_position, s_radius, eps);

        \\# debug function to ensure that the divisor actually truly generates the ideal
        \\print("verifying element after adjust");verify_generator_with_list(G, baby_divisor[1], compactTracking);
        logdist = trackingLog;

        \\# identify all ideal to be scanned plus the corresponding u
        if (mapisdefined(scanIdeals, baby_divisor[1], &distanceList),
            repeat_counter+=1;
            if(!is_repeat_babystock(distanceList, logdist, eps),
                \\print("new babystock element: ", precision(logdist, 10));
                listput(~distanceList, logdist);
                mapput(~scanIdeals, baby_divisor[1], distanceList);

                templist = mapget(idealCompactGenerator, baby_divisor[1]);
                listput(~templist, compactTracking);
                mapput(~idealCompactGenerator,baby_divisor[1],templist );
            );
        ,
            \\ # Important! scanIdeals pushes the nu value of the first
            \\ # occurrence. This is accounted for in overlap_scanball
            mapput(~scanIdeals, baby_divisor[1] ,List([baby_divisor[2], logdist]));
            \\print("new elements: ",  baby_divisor[1], "  ",precision(logdist,10));
            mapput(~idealCompactGenerator, baby_divisor[1], List([compactTracking] ));
        );

        \\# increase the tracking variable and update directions
        place_marker = increment_with_place_marker(~denoms, ~web_coords);
        updateDirections(~directions, ~place_marker);

        \\\ regenerate the logs of the babysteps
        if ((ctr%500 == 0) && ctr > 0,
            compactTracking = emptyCompactList;
            base_value = compact_rep_full_input(G, trackingLog, baby_divisor[1], eps, 1, 2);

            GP_ASSERT_NEAR(normlp(log_from_cpct( G, base_value)-trackingLog), 0, 0.00001);
            trackingLog = log_from_cpct( G, base_value);
        );
        ctr++;
        if((ctr % 2000) == 0,
            print(ctr, "  ", web_coords);
            baby_tn = getabstime();
            baby_tmid = baby_tn;
            if((timeout > 0)&&(baby_tn - start_time > timeout),
                write(OUTFILE_BS, "babystock computation ", (baby_tn - start_time)/60000.0, " mins. Exceeds timeout.");return([lattice_lambda, []]);
            );
        );
    );
    \\default(realbitprecision, mainbitprecision);
    \\# go through each unique ideal and enumerate once. The results are multipled
    \\# by each associated u to account for all the distinct minima
    print("Babystock: enumerating");
    if (length(scanIdeals) < 1, return([lattice_lambda, []]); );

    scanIdealsMatrix = Mat(scanIdeals)[,1];
    print("precision before scanball ", default(realbitprecision));
    for(i = 1, length(scanIdealsMatrix),
        distanceList = mapget(scanIdeals, scanIdealsMatrix[i]);

        \\cpct_list = mapget(idealCompactGenerator, scanIdealsMatrix[i]);
        \\GP_ASSERT_EQ(length(distanceList)-1, length(cpct_list));
        if (length(distanceList) > 2, print("Collision found within babystocks"););
        nu = distanceList[1];
        overlap_scanball(~G, ~baby_hashmap, ~scanIdealsMatrix[i], ~nu, ~distanceList, scan_radius, eps, ~repeated_minima);
    );
    my(denom_product = 1);
    for(i=1, r, denom_product*=denoms[i]);
    print("scan time: ", getabstime() - baby_t1, "  ", denom_product);
    [lattice_lambda, newctr] = check_units_bstock(~baby_hashmap,~lattice_lambda,~G,eps);
    if(newctr != 0,
        print("Found a unit in babystock. new reg =", precision(abs(matdet(lattice_lambda)),10) );
        print("Babysteps stored ", length(Mat(baby_hashmap)~), "  Number of times ideals were repeats: ", repeat_counter);
        return([lattice_lambda, [1] ]);
    );
    return([lattice_lambda, []]);
}

\\ #Subalgorithm of bsgs. computes baby steps incrementally using ideal multiplication
\\ #variable for storage option
incremental_baby_steps_compact(y, ~lattice_lambda, ~giant_legs,\
                        ~baby_hashmap, G, scanballRadius, eps, storage, ~cpct_hashmap, outFileInfo=[]) =
{
    if(storage != "COMPACT", print("invalid storage format for this babystep routine. Compact only"); break;);
    my(timeout, OUTFILE_BS);
    if(length(outFileInfo) == 2,
        timeout = outFileInfo[2];
        OUTFILE_BS = outFileInfo[1];
    ,
        timeout = 0;
        OUTFILE_BS = 0;
    );
    print("Storing baby steps via ideal increments");
    GP_ASSERT_TRUE(eps > 0);
    my(
        field_deg = poldegree(G.pol),   r = G.r1+G.r2-1,
        zero_vec = vector(r, i, 0),     identity = matid(field_deg),
        web_coords = zero_vec,          place_marker = r,
        directions = vector(r, i, 1),
        denoms,
        scan_shrink_factor,
        babystock_t,
        expected_position = vector(r+1, i, 0)~,
        start_time
    );
    GP_ASSERT_EQ(r, length(giant_legs));
    REQ_BS = babystockPrecision(G, giant_legs);

    default(realbitprecision, REQ_BS);  \\\ change precision, switches back in giant step algs
    start_time = getabstime();

    \\\ # compute vectors which determine what the babystock look like
    \\\ # the vector denoms informs how many babysteps are done in each direction to cover the babystock
    [babystock_t, denoms] = initialize_babystock_edges(~giant_legs, scanballRadius, r);
    increment_coordinates(denoms, web_coords);

    my(
        \\ #vectors of length r which are used to compute an adjacent element in
        \\ #a particular direction in the logarithm lattice
        \\ #direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
        direction_elements,
        inverse_direction_elements,
        scanIdeals = Map(),
        distanceList = List(),
        scanIdealsMatrix,
        repeat_counter = 0,
        repeated_minima = 0,
        logdist,
        compactTracking = List(),
        emptyCompactList,
        trackingLog = vector(r+1, i, 0),
        idealCompactGenerator = Map(),
        s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc)),
        base_value = [[1], [1]];
    );
    \\ # multiplicands to move in each 'direction', forward and backwards
    \\ # each one is in compact representation form
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, babystock_t, field_deg, eps);

    \\ # initialize the current baby step element and the logarithm
    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );
    emptyCompactList = compactTracking;

    \\ # initialize the first baby-step divisor
    baby_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN("Additional timing variables and file write in jump. Beginning baby step loop ", denoms);
    my(baby_t1, baby_tn, baby_tmid, ctr);
    baby_t1 = getabstime();
    baby_tmid = baby_t1;

    t_divisor_update = 0;
    t_divisor_reduce = 0;
    t_ideal_collection = 0;
    t_enumeration = 0;

    \\##########################################################################
    \\# Begin scanning over the babystock region
    \\# increments are done modulo the product of avec, returns to zero when
    \\# all elements have been visited
    while(web_coords != zero_vec,

        /*
        \\print(web_coords , "  BDC: ", precision(logdist,10));
        \\print("verifying loop start. reset commented out");verify_generator_with_list(G, baby_divisor[1], compactTracking);
        if(place_marker == 1 && ((web_coords[1]%2) == 0),
            for(i = 1, length(directions), GP_ASSERT_TRUE(directions[place_marker] == 1));
            compactTracking = List();
            for( i=1, length(direction_elements),
                listput(~compactTracking, [direction_elements[i][3], 0]);
            );
            compactTracking[1][2] = web_coords[1];
            logdist = web_coords[1]*log_from_cpct(G, direction_elements[place_marker][3]);
            print("BD: ",  precision(baby_divisor,10));
            target = direction_elements[place_marker][2];
            for(l = 1, web_coords[1], target = pointwise_vector_mul(target, direction_elements[place_marker][2]));
            print("reset vals: ", idealpow(G, direction_elements[place_marker][1],web_coords[1]), "  ", precision(target,10));
        );
        */

        time_divisor_update1 = getabstime();
        \\\ # Take a step with ideal and update relevant variables
        GP_ASSERT_TRUE(storage == "COMPACT");
        if(directions[place_marker] == 1,
            baby_divisor = [idealmul(G, baby_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],direction_elements[place_marker][2] )~];

            expected_position += babystock_t[,place_marker];
            compactTracking[place_marker][2] += 1;

        ,\\else
            baby_divisor = [idealmul(G, baby_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],inverse_direction_elements[place_marker][2] )~];

            expected_position -= babystock_t[,place_marker];
            compactTracking[place_marker][2] -= 1;
        );

        time_divisor_update2 = getabstime();
        t_divisor_update += (time_divisor_update2 - time_divisor_update1);
        time_reduce1 = getabstime();

        \\\ # reduce and update variables
        wasUpdated = get_next_giant_divisor_cpct(G, ~baby_divisor, ~compactTracking);
        if(storage == "LOG" && wasUpdated,
            trackingLog += log(abs(nfeltembed(G, compactTracking[length(compactTracking)][1])));
        );
        if(storage == "COMPACT",
            trackingLog = log_from_cpct(G, base_value)
        );
        \\\ # ensure baby-step divisor is within the expected range of the web point
        [baby_divisor, tempLog] = adjust_giant_step_cpct(~G, ~baby_divisor,~compactTracking, ~trackingLog, ~expected_position, s_radius, eps, storage);

        time_reduce2 = getabstime();
        t_divisor_reduce += (time_reduce2 - time_reduce1);
        time_ideals1 = getabstime();

        if(storage == "COMPACT",
            logdist = tempLog;
            original_precision = default(realbitprecision);
            \\default(realbitprecision, 50);
            \\GP_ASSERT_NEAR(norml2(trackerLogarithm(G, ~compactTracking, r)- tempLog+trackingLog), 0, 2^(-40));
            \\default(realbitprecision,original_precision);
        ,
            logdist = trackingLog;
        );

        \\# identify ideals to be scanned plus the corresponding u
        if (mapisdefined(scanIdeals, baby_divisor[1], &distanceList),
            repeat_counter+=1;
            if(!is_repeat_babystock(distanceList, logdist, eps),
                listput(~distanceList, logdist);
                mapput(~scanIdeals, baby_divisor[1], distanceList);

                if (storage == "COMPACT",
                    templist = mapget(idealCompactGenerator, baby_divisor[1]);
                    listput(~templist, [base_value, compactTracking]);
                    GP_ASSERT_TRUE(length(templist[length(templist)]) == 2 );
                    mapput(~idealCompactGenerator,baby_divisor[1],templist );
                );
            );
        ,
            \\ # Important! scanIdeals pushes the nu value of the first
            \\ # occurrence. This is accounted for in overlap_scanball
            mapput(~scanIdeals, baby_divisor[1] ,List([baby_divisor[2], logdist]));
            if (storage == "COMPACT",
                mapput(~idealCompactGenerator, baby_divisor[1], List([[base_value,compactTracking]] ));
                GP_ASSERT_EQ(length(mapget(idealCompactGenerator,baby_divisor[1])),1);
            );
        );

        \\# increase the tracking variable and update directions
        place_marker = increment_with_place_marker(~denoms, ~web_coords);
        updateDirections(~directions, ~place_marker);
        ctr++;
        if((ctr % 100) == 0,
            \\GP_ASSERT_NEAR(norml2(trackerLogarithm(G, ~compactTracking, r)- tempLog+trackingLog), 0, 2^(-40));
            GP_ASSERT_NEAR(norml2(trackerLogarithm(G, ~compactTracking, r)- tempLog+trackingLog), 0, 2^(-40));
            base_value = compact_rep_full_input(G, tempLog, matid(field_deg), eps);
            compactTracking = emptyCompactList;
            \\print("tempLog  ", precision(tempLog, 10));

            print(ctr, "  ", web_coords);
            print("Total Timing: Update: ", t_divisor_update, "  Reduction: ", t_divisor_reduce, "  ideal collect: ", t_ideal_collection,
                " Enum: ", t_enumeration) ;
            baby_tn = getabstime();
            baby_tmid = baby_tn;
            if((timeout > 0)&&(baby_tn - start_time > timeout),
                write(OUTFILE_BS, "babystock computation ", (baby_tn - start_time)/60000.0, " mins. Exceeds timeout.");return([lattice_lambda, []]);
            );
        );
        time_ideals2 = getabstime();
        t_ideal_collection += (time_ideals2 - time_ideals1);
    );

    time_enumerate1 = getabstime();
    \\# go through each unique ideal and enumerate once. The results are multipled
    \\# by each associated u to account for all the distinct minima
    print("Babystock: enumerating");
    if (length(scanIdeals) < 1, return([lattice_lambda, []]); );

    scanIdealsMatrix = Mat(scanIdeals)[,1];
    for(i = 1, length(scanIdealsMatrix),
        distanceList = mapget(scanIdeals, scanIdealsMatrix[i]);
        if (storage == "COMPACT",
            cpct_list = mapget(idealCompactGenerator, scanIdealsMatrix[i]);
            GP_ASSERT_EQ(length(distanceList)-1, length(cpct_list));
        );
        if (length(distanceList) > 2, print("Warning collision found within babystocks"););
        nu = distanceList[1];
        if(storage == "COMPACT",
            compact_storage_overlap_scanball(~G, ~baby_hashmap, ~scanIdealsMatrix[i], ~nu, ~distanceList, scan_shrink_factor, eps, ~repeated_minima, cpct_list, ~cpct_hashmap);
        ,
            overlap_scanball(~G, ~baby_hashmap, ~scanIdealsMatrix[i], ~nu, ~distanceList, scan_shrink_factor, eps, ~repeated_minima);
        );
    );
    my(denom_product = 1); for(i=1, r, denom_product*=denoms[i]);
    if(storage == "COMPACT",
        print("hashmap size: ", length(cpct_hashmap)); GP_ASSERT_EQ(length(baby_hashmap), length(cpct_hashmap));
    );
    print("scan time: ", getabstime() - baby_t1, "  ", denom_product);
    if(storage == "COMPACT",
        [lattice_lambda, newctr] = check_compact_bstock(~cpct_hashmap, ~lattice_lambda, ~G, eps);
    ,
        [lattice_lambda, newctr] = check_units_bstock(~baby_hashmap, ~lattice_lambda, ~G, eps);
    );

    if(newctr != 0,
        print("Found a unit in babystock. new reg =", precision(abs(matdet(lattice_lambda)),10) );
        print("Babysteps stored ", length(Mat(baby_hashmap)~), "  Number of times ideals were repeats: ", repeat_counter);
        return([lattice_lambda, [1] ]);
    );
    time_enumerate2 = getabstime();

    t_enumeration = time_enumerate2 - time_enumerate1;

    print("Total Timing: Update: ", t_divisor_update, "  Reduction: ", t_divisor_reduce, "  ideal collect: ", t_ideal_collection,
        " Enum: ", t_enumeration) ;

    return([lattice_lambda, []]);
}
