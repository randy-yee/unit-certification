read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
FIRST = 1;
VERIFY_GENERATORS = 1;
adjusttime1 = 0;


\\ INPUT:
\\ - G a number field
\\ - giant_divisor is a 3-element list consisting of
\\      an ideal, a vector of size G.r1+G.r1 corresponding to valuation of an element
\\      a compact representation
\\ - expected position is the coordinates in R^r of your giant lattice element
\\ - distance_ok is the limit of the acceptable distance from the divisor to the target position
adjust_giant_step_cpct(~G, ~giant_divisor, ~tracker, ~trackerLog, ~expected_position, eps)={
    a_time1 = getabstime();
    my(
        r = G.r1 + G.r2 -1,
        divisor_distance = expected_position - trackerLog~,
        adjustment_divisor, new_divisor,
        new_distance,
        newFactor,
        s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    );


    a_time2 = getabstime();
    adjusttime1 += (a_time2 - a_time1);
    if(sqrt(norml2(divisor_distance)) < s_radius,
        return(giant_divisor);
    , \\else
        \\print("\ntarget position = ", precision(expected_position,10), "\nOriginal Distance from Target ", precision(norml2(divisor_distance),10) );
        for(i=1, 2,

            adjustment_divisor = get_nearby_rdivisor(G, matid(poldegree(G.pol)), divisor_distance, i%2);
            if (sqrt(norml2(adjustment_divisor[3])) < eps,
                return(giant_divisor);
            ,
                new_divisor = [idealmul(G, giant_divisor[1], adjustment_divisor[1]),
                    pointwise_vector_mul(giant_divisor[2],adjustment_divisor[2] )~,
                    mul_compact(G, giant_divisor[3], [ [numerator(adjustment_divisor[4])],[denominator(adjustment_divisor[4])] ])  ];

                reduced_product = reddiv_compact(new_divisor[1], new_divisor[2],G, G[5][1] );

                reduced_product_cpct = mul_compact(G, new_divisor[3], [[numerator(reduced_product[4])],[denominator(reduced_product[4])]] );
                new_divisor[3] = reduced_product_cpct;

                newFactor = nfeltmul(G, adjustment_divisor[4], reduced_product[4]);
                logNewFactor = log(abs(nfeltembed(G,newFactor) ))[1..r];

                new_distance = norml2(expected_position - logNewFactor);
                \\print("Adjusted distance from target = ", precision(new_distance, 10) );
                if(new_distance < norml2(divisor_distance)+eps,
                    \\print("returning adjusted divisor");

                    listput(~tracker, [newFactor, 1]);
                    trackerLog += logNewFactor;
                    return(new_divisor);
                );
            );
        );
        \\ adjustment fails, so we should use jump to compute something close
        gstep_divisor = jump_compact(matid(length(G.zk)), expected_position, G, length(G.zk), eps);
        print("WARNING: recomputing tracking divisor with jump. Reset the tracker also");
        return(gstep_divisor);
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
\\ bmap is passed by reference, and any new minima are added to it
overlap_scanball(~G, ~bmap, ~y, ~u, ~log_distance_list, ball_distance, eps, ~repeated_minima)={

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


\\ #Subalgorithm of bsgs. computes baby steps incrementally using ideal
\\ #multiplication to obtain adjacent elements.
\\ #Compact Represenation Version
incremental_baby_steps(y, ~lattice_lambda, ~giant_legs, ~baby_hashmap, G, eps)=
{
    print("baby steps using increments");
    my(
        field_deg = poldegree(G.pol),
        r = G.r1 + G.r2 -1,
        zero_vec = vector(r, i , 0),
        web_coords = zero_vec,
        current_giant_vec = zero_vec,
        identity = matid(field_deg),
        place_marker = r,
        directions,
        denoms,
        scan_shrink_factor,
        babystock_t,
        expected_position
    );
    GP_ASSERT_EQ(r, length(giant_legs));

    \\ # setup up baby stock scanning region
    denoms = vector( length(lattice_lambda), i, ceil(norml2(giant_legs[,i])));
    scan_shrink_factor = 1;
    denoms*= scan_shrink_factor;
    babystock_t = giant_legs;
    for(i=1, length(babystock_t),
        babystock_t[,i]=babystock_t[,i]/denoms[i];
    );
    expected_position = vector(r, i, 0)~;

    directions = vector(length(giant_legs), i , 1);
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
        trackingLog = vector(r, i, 0);
    );
    \\ #for some reason I have to flip these. giantstep() may actually return an inverse
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, babystock_t, field_deg, eps);

    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );
    baby_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Additional timing variables and file write in jump");
    my(baby_t1, baby_tn, baby_tmid, ctr);
    baby_t1 = getabstime();
    baby_tmid = baby_t1;

    \\# increments are done modulo the product of avec, returns to zero when
    \\# all elements have been visited

    while(web_coords != zero_vec,
        if(directions[place_marker] == 1,
            baby_divisor = [idealmul(G, baby_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],direction_elements[place_marker][2] )~,
                mul_compact(G, baby_divisor[3],direction_elements[place_marker][3])  ];
            expected_position += babystock_t[,place_marker];
            compactTracking[place_marker][2] += 1;
            trackingLog += log_from_cpct(G, direction_elements[place_marker][3])[1..r];
        ,\\else
            baby_divisor = [idealmul(G, baby_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(baby_divisor[2],inverse_direction_elements[place_marker][2] )~,
                mul_compact(G, baby_divisor[3], inverse_direction_elements[place_marker][3])  ];
            expected_position -= babystock_t[,place_marker];
            compactTracking[place_marker][2] -= 1;
            trackingLog -= log_from_cpct(G, direction_elements[place_marker][3])[1..r];
        );

        baby_divisor = get_next_giant_divisor_cpct(G, ~baby_divisor, ~compactTracking);
        trackingLog += log(abs(nfeltembed(G, compactTracking[length(compactTracking)][1])))[1..r];

        baby_divisor = adjust_giant_step_cpct(~G, ~baby_divisor,~compactTracking, ~trackingLog, ~expected_position, eps);
        logdist = log_from_cpct(G, baby_divisor[3]);
        \\GP_ASSERT_VEC_NEAR(logdist[1..r], trackerLogarithm(G, ~compactTracking, r)[1..r], 0.0000001);
        \\print("Compare: \n", precision(logdist,10), "\n", precision(trackingLog,10));
        \\verify_generator(G, baby_divisor[1], baby_divisor[3]);

        \\# identify all ideal to be scanned plus the corresponding u
        if (mapisdefined(scanIdeals, baby_divisor[1], &distanceList),
            repeat_counter+=1;
            if(!is_repeat_babystock(distanceList, logdist),
                listput(~distanceList, logdist);
                mapput(~scanIdeals, baby_divisor[1], distanceList);
            );
        ,
            \\ # Important! scanIdeals pushes the nu value of the first
            \\ # occurrence. This is accounted for in overlap_scanball
            mapput(~scanIdeals, baby_divisor[1] ,List([baby_divisor[2], logdist]));
        );

        \\# increase the tracking variable and update directions
        place_marker = increment_with_place_marker(~denoms, ~web_coords);
        updateDirections(~directions, ~place_marker);

        ctr++;
        if((ctr % 1000) == 0,
            print(ctr);
            baby_tn = getabstime();
            baby_tmid = baby_tn;
        );
    );

    \\# go through each unique ideal and enumerate once. The results are multipled
    \\# by each associated u to account for all the distinct minima
    print("Babystock: enumerating");
    if (length(scanIdeals) < 1, return([lattice_lambda, []]); );
    scanIdealsMatrix = Mat(scanIdeals)[,1];
    for(i = 1, length(scanIdealsMatrix),
        distanceList = mapget(scanIdeals, scanIdealsMatrix[i]);
        print("elements per ideal: ", length(distanceList));
        nu = distanceList[1];
        overlap_scanball(~G, ~baby_hashmap, ~scanIdealsMatrix[i], ~nu, ~distanceList, scan_shrink_factor, eps, ~repeated_minima);
    );
    print("scan time: ", getabstime() - baby_t1);
    [lattice_lambda, newctr] = check_units_bstock(~baby_hashmap,~lattice_lambda,~G,eps);
    if(newctr != 0,
        print("Found a unit in babystock. new reg =", precision(abs(matdet(lattice_lambda)),10) );
        print("Babysteps stored ", length(Mat(baby_hashmap)~), "  Number of times ideals were repeats: ", repeat_counter);
        return([lattice_lambda, [1] ]);
    );
    return([lattice_lambda, []]);
}

\\#brief - Ensure that ideal = Ok/cpct
\\#param G is a number field
\\#param ideal is an ideal of G
\\#param cpct is a compact representation
\\#detail ASSERT will fail if the ideal Ok/cpct is not equal to input ideal
verify_generator(G, ideal, cpct)=
{
    if(VERIFY_GENERATORS,
        VERIFY_GENERATORS = 0;
        print("WARNING: Verifying ideal generators, expect performance decrease");
    );
    my(
        ndeg = poldegree(G.pol),
        Ok = matid(ndeg),
        tempIdeal = Ok,
        quotientIdeal = idealdiv(G, Ok, compact_reconstruct(G, cpct[1], cpct[2]))
    );

    for(i=1, length(cpct[1]),
        tempIdeal = idealmul(G, tempIdeal, tempIdeal);
        tempIdeal = idealdiv(G,tempIdeal, nfeltdiv(G, cpct[1][i], cpct[2][i]));
    );
    GP_ASSERT_EQ(tempIdeal, quotientIdeal);
    GP_ASSERT_EQ(ideal, quotientIdeal);
}
