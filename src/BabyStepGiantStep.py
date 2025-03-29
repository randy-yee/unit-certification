read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/CompactRepresentation.py");
read("src/Neighbours.py")
read("src/bounds.gp")
read("src/BabyStep.py")
read("src/BSGSHelper.py")
read("src/BSGS_ParameterSelection.py")
/*
\\ BSGS related routines
increment_coordinates
increment_with_place_marker
get_axis_aligned_box
update_directions

jump_giant_steps
get_nearby_rdivisor
multiply_and_reduce_divisors
adjust_giant_step
is_repeat_babystock
expand_babystock_region
get_giant_step_increment_vectors
get_next_giant_divisor
one_massive_qfminim

bsgs
*/
DEBUG_BSGS = 0;
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ BBBBBBB    SSSSSS    GGGGGG     SSSSSS
\\ BB    BB  SS        GG    GG   SS
\\ BBBBBBB    SSSSSS   GG          SSSSSS
\\ BB    BB        SS  GG   GGGG        SS
\\ BBBBBBB    SSSSSS    GGGGGG     SSSSSS
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\ Below are all function primarily related to the babystep giant step algorithm
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


/******************************************************************************/

\\# used to determine if an element's real embedding has already been encountered
\\# existing_entry is a list of embeddings
\\# candidate embedding is the embedding we wish to check is repeated
\\# eps is the precision for comparison
is_repeat_babystock(~existing_entry, ~candidate_embedding, eps)={
    my(repeatflag = 0);
    for(t=1, length(existing_entry),
        if(samevecs(existing_entry[t], candidate_embedding,eps), repeatflag = 1);
    );
    return(repeatflag);
}


babystock_scan_jump(y,~L,~giant_legs,~baby_hashmap,G, scanballRadius, eps)={
    print("Babystock: jump version");
    my(next_coords = [],
        zerovec = vector(length(L), i, 0),
        web_coords = zerovec,                                                   \\ used to track the coefficients for the web points
        web_step = [[],[]],                                                     \\ stores the fixed elements we use to modify the web point
        webpoint,
        web_distance,
        web_increments,
        box_volume,
        box_subdivisions,
        exp_webpoint,
        ideal_J, nu, logdist, beta,
        directions,
        ball_minima,
        newctr = 0,
        existing_entry,
        field_deg = poldegree(G.pol),
        identity_n = matid(field_deg),
        scan_radius = 0,
        scan_shrink_factor
    );
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Establish a number of parameters including the web of regularly distributed points,
    \\ and some values that allow us to move easily between each of the points sucessively.
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    repeat_counter = 0;
    \\ subdivide the babystock region based on the lengths of the defining vectors
    denoms = vector( length(L), i, ceil(norml2(giant_legs[,i])));
    scan_shrink_factor = scanballRadius;
    denoms*= scan_shrink_factor;
    babystock_t = giant_legs;
    for(i=1, length(babystock_t),
        babystock_t[,i]=babystock_t[,i]/denoms[i];
    );

    xvector = vector(matsize(giant_legs)[2], i,0)~;
    increment_coordinates(denoms, web_coords);

    scanIdeals = Map();
    my(
        repeated_minima = 0,
        distanceList = List(),
        scanIdealsMatrix
    );
    \\ #precomputation of ideals to scan plus log distances

    time2 = getabstime();
    print("Babystock: precompute ideals");
    while(web_coords != zerovec,
        xvector = babystock_t*web_coords~;
        increment_coordinates(denoms, web_coords);
        [ideal_J, nu, logdist] = giantstep(y, xvector, G, field_deg, eps);
        if (mapisdefined(scanIdeals, ideal_J, &distanceList),
            repeat_counter+=1;
            if(!is_repeat_babystock(distanceList, logdist),
                listput(~distanceList, logdist);
                mapput(scanIdeals, ideal_J, distanceList);
            );
        ,
            \\ # Important! scanIdeals pushes the nu value of the first
            \\ # occurrence. This is accounted for in overlap_scanball
            mapput(scanIdeals, ideal_J ,List([nu, logdist]));
        );
    );
    print("Babystock: enumerating");
    if(length(scanIdeals)>0,
        scanIdealsMatrix = Mat(scanIdeals)[,1];
    ,\\else
        \\# if there are no babystock elements, just return L and empty list
        return([L, []]);
    );

    for(i = 1, length(scanIdealsMatrix),
        distanceList = mapget(scanIdeals, scanIdealsMatrix[i]);
        nu = distanceList[1];
        overlap_scanball(G, ~baby_hashmap, scanIdealsMatrix[i], nu, distanceList, scan_shrink_factor, eps, repeated_minima);
    );
    \\print("new method time ", getabstime() - time2);
    /*
    time1 = getabstime();increment_coordinates(denoms, web_coords);
    while(web_coords != zerovec,                                               \\ loops over all of the web points in the babystock region

        xvector = babystock_t*web_coords~;                                     \\ this holds the web point to be scanned
        increment_coordinates(denoms, web_coords);
        [ideal_J, nu, logdist] = giantstep(y, xvector, G, field_deg, eps);

        ball_minima = scanball_map(~G, ~baby_hashmap, ideal_J, nu, logdist[1..length(L)], scan_shrink_factor, eps, ~repeated_minima);

        if (mapisdefined(scanIdeals, ideal_J),
            repeat_counter+=1;
            mapget(scanIdeals, ideal_J);
        ,
            mapput(scanIdeals, ideal_J ,1);
        );

        [L, newctr] = check_units_bstock(baby_hashmap,L,G,eps);

        if(newctr != 0,
            print("Found a unit in babystock. new reg =", precision(abs(matdet(L)),10) );
            print("Babysteps stored ", length(Mat(baby_hashmap)~), "  Number of times ideals were repeats: ", repeat_counter);

            return([L, [1] ]);
        );

        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\ Take the elements found on the ball near the current point in the web, and enter them into the hash map
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        for(ctr=1, length(ball_minima),
            if(mapisdefined(baby_hashmap, ball_minima[ctr][1], &existing_entry),
                repeatflag = is_repeat_babystock(existing_entry, ball_minima[ctr][2],eps);
                if(repeatflag==0,
                    listput( ~existing_entry, ball_minima[ctr][2]);
                    mapput(baby_hashmap, ball_minima[ctr][1], existing_entry );
                );
            ,\\else
                mapput(baby_hashmap, ball_minima[ctr][1], List([ball_minima[ctr][2]]));
            );
        );
    );
    */

    print("Babysteps stored ", length(Mat(baby_hashmap)~), "  Number of times ideals were repeats: ", repeat_counter);
    return([L, []]);
}


\\ subroutine of bsgs - performs giant steps using jump algorithm one at a time and checks with babystock
\\ INPUT:
\\ - G a number field and lat_lambda is unit lattice
\\ - gs_sublattice is a sublattice of the fundamental domain of lat_lambda
\\ - bstock is a Map object of the babystock. Keys are ideals, values are Psi values of generating elements
\\ - avec indicates the number of subdivisions of lat_lambda were taken in each coordinate to get gs_sublattice
\\ - eps is an error value
\\ RETURN: a potentially expanded unit_lattice.
jump_giant_steps(~G, ~lat_lambda, ~gs_sublattice, ~bstock, avec, eps)={
    my(
        field_deg = poldegree(G.pol), r = G.r1 + G.r2 -1,
        identity = matid(field_deg), zero_vec = vector(r, i , 0),
        giant_coeffs = zero_vec, current_giant_vec = zero_vec,
        giant_divisor,
        matches, new_vec,
        eps_sqr = (eps*10)^2
    );
    giant_coeffs[r] = 1;
    SCREEN(0, "Additional timing variables and file write in jump");
    my(giant_t1, giant_tn, giant_tmid, ctr);
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    while(giant_coeffs != zero_vec,

        current_giant_vec = (gs_sublattice*giant_coeffs~)~;
        giant_divisor = giantstep(identity, current_giant_vec, G, field_deg, eps);

        if(mapisdefined(bstock, giant_divisor[1]),
            matches = mapget(bstock, giant_divisor[1]);                         \\ list of existing babystock matched elements
            for(i=1, length(matches),
                new_vec = giant_divisor[3][1..r] - matches[i];                  \\ compute difference
                if(norml2(new_vec) > eps_sqr && is_vec_in_lattice(new_vec~,lat_lambda,eps)==0,
                    print("giant-step jump",matsize(lat_lambda));
                    lat_lambda = my_mlll( matconcat([lat_lambda, new_vec~]),eps);
                    if(1, print("New element found. New regulator = ", precision(matdet(lat_lambda),10)););
                );
            );
        );
        increment_coordinates(avec, giant_coeffs);

        ctr++;
        if(ctr %1000 ==0,
            print(ctr);
            giant_tn = getabstime();
            \\write("data/jump-timing.txt", ctr "jumps took: ", (giant_tn- giant_tmid), "milliseconds");
            giant_tmid = giant_tn;
        )
    );
    return(lat_lambda);
};

\\ used to compute the multiplying elements for the incremental giant steps
\\ strategy. Obtains r 'forward' elements, and r 'backward' elements
get_giant_step_increment_vectors_compact(G, giant_sublattice, field_deg, eps)={
    my(
        field_deg = poldegree(G.pol),
        giant_step_ideals = List(),
        giant_step_ideals_inverse = List(),
        new_vec, gstep_divisor,
        new_alpha, new_denom, inverse_elt,
        urank = G.r1 + G.r2-1
    );

    \\ loop through each dimension of the lattice
    \\#print("ISSUE: Sometimes when the vectors are too small, the minimum selected is 1, which causes issues");
    \\#after testing, it seems to still be faster than jump so I will leave it as is
    for(j = 1, length(giant_sublattice),
        new_vec = giant_sublattice[,j];

        gstep_divisor = jump_compact(matid(field_deg), new_vec, G, field_deg, eps);
        listput(~giant_step_ideals, gstep_divisor);

        GP_ASSERT_TRUE(type(DEBUG_BSGS)=="t_INT");
        if(DEBUG_BSGS,

            print("1:" , gstep_divisor[1], " ", precision(gstep_divisor[2],10));
            print(gstep_divisor[3]);
            print("BSGS-DEBUG ", precision(log_from_cpct(G, gstep_divisor[3]),10), "  ",precision(giant_sublattice[,j],10), "  ", precision(log(abs(G.disc)^(1/2)),10) );
            GP_ASSERT_TRUE(sqrt(norml2(log_from_cpct(G, gstep_divisor[3])~-giant_sublattice[,j]))<log(abs(G.disc)^(1/2)) );
        );

        gstep_divisor[3] = invert_compact(G, gstep_divisor[3]);
        listput(~giant_step_ideals_inverse, [idealinv(G,gstep_divisor[1]), invert_coordinates(gstep_divisor[2]),gstep_divisor[3] ] );


        \\debug_compare(idealdiv(G, matid(field_deg), compact_reconstruct(G,giant_step_ideals_inverse[j][3][1], giant_step_ideals_inverse[j][3][2])), giant_step_ideals_inverse[j][1]);
        \\ #check that the cpct reps are actually inverses.
        \\ #also check that the reps are near the intended log lattice point
        \\print("WARNING: Verifying increment ideals");
        \\myprod = mul_compact(G, giant_step_ideals[j][3], giant_step_ideals_inverse[j][3]);
        \\GP_ASSERT_EQ(compact_reconstruct(G, myprod[1], myprod[2]), 1);
        \\debug_compare(idealdiv(G, matid(field_deg), compact_reconstruct(G,giant_step_ideals[j][3][1], giant_step_ideals[j][3][2])), giant_step_ideals[j][1]);
    );
    return([giant_step_ideals, giant_step_ideals_inverse]);
}

updateDirections(~directions_vec, place_marker)=
{
    for(j = place_marker+1, length(directions_vec),
        directions_vec[j] = (directions_vec[j] %2) +1;
    );
}

check_lambda_is_valid(~latt)=
{
    if (matsize(latt)[1] != matsize(latt)[2],
        print("ERROR: mlll failed to properly determine new basis in giant steps.");
        breakpoint();
    );
    if(DEBUG_BSGS>0, print("new regulator = ", precision(matdet(latt),10)););
}

collision_check_cpct(G, ~lattice_lambda, r, eps,\
                 ~matches_cpct, ~giant_divisor, ~compactTracking)=
{
    gd_log = trackerLogarithm(G, ~compactTracking, r);
    if(normlp(gd_log)< 2^(-10),
        log_inf_v = 1,
        log_inf_v = log(normlp(gd_log)) );
    temp_precision = prec_compact(poldegree(G.pol), log(abs(G.disc))/log(2), log_inf_v);
    old_precision = change_precision(temp_precision);

    for(i=1, length(matches_cpct),

        baby_log = trackerLogarithm(G, ~matches_cpct[i], r);
        new_vec = gd_log - baby_log;
        if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec[1..r]~,lattice_lambda,eps)==0,
            if(DEBUG_BSGS>1,
                lattice_lambda_copy = lattice_lambda;
                verify_lattice_containment(~G, ~new_vec);
            );

            lattice_lambda = my_mlll( matconcat([lattice_lambda, new_vec[1..r]~]),eps);
            check_lambda_is_valid(~lattice_lambda);
        );
    );
    change_precision(old_precision);
    return(lattice_lambda);
}

\\ # this function takes as input the log of the current giant step element
\\ # the idea being that it is computed in the calling function
collision_check_bootstrap(G, ~lattice_lambda, r,\
                 ~matches_cpct, ~giant_divisor, ~giant_element_logarithm)=
{
    my(gd_log, eps, temp_precision, old_precision, baby_log, new_vec);
    gd_log = giant_element_logarithm;
    eps = 2^(- (ceil(log(4*poldegree(G.pol))/log(2))+1 + ceil(0.5*log(poldegree(G.pol))/log(2))) );

    if(normlp(gd_log)< 2^(-10),
        log_inf_v = 1,
        log_inf_v = log(normlp(gd_log)) );
    temp_precision = REQ_JUMP(G, gd_log);
    old_precision = change_precision(temp_precision);

    for(i=1, length(matches_cpct),

        baby_log = log_from_cpct(G,~matches_cpct[i][1])+trackerLogarithm(G, ~matches_cpct[i][2], r);
        \\print(precision(log_from_cpct(G,~matches_cpct[i][1]),10) , " + " precision(trackerLogarithm(G, ~matches_cpct[i][2], r),10) );
        new_vec = gd_log - baby_log;
        if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec[1..r]~,lattice_lambda,eps)==0,

            if(DEBUG_BSGS>1,
                lattice_lambda_copy = lattice_lambda;
                verify_lattice_containment(~G, ~new_vec);
            );

            lattice_lambda = my_mlll( matconcat([lattice_lambda, new_vec[1..r]~]),eps);
            check_lambda_is_valid(~lattice_lambda);
        );
    );
    change_precision(old_precision);
    return(lattice_lambda);
}
\\ # this function takes as input the log of the current giant step element
\\ # the idea being that it is computed in the calling function
collision_check_bootstrap_log(G, ~lattice_lambda, r,\
                 ~matches, ~giant_divisor, ~giant_element_logarithm)=
{
    my(gd_log, eps, temp_precision, old_precision, baby_div, baby_log, new_vec);
    gd_log = giant_element_logarithm;
    eps = 2^(- (ceil(log(4*poldegree(G.pol))/log(2))+1 + ceil(0.5*log(poldegree(G.pol))/log(2))) );

    if(normlp(gd_log)< 2^(-10),
        log_inf_v = 1,
        log_inf_v = log(normlp(gd_log)) );
    temp_precision = prec_compact(poldegree(G.pol), log(abs(G.disc))/log(2), log_inf_v);
    old_precision = change_precision(temp_precision);

    for(i=1, length(matches),
        \\# do this to recalculate the log vector to greater precision
        baby_div = compact_rep_full_input(G, matches[i],giant_divisor[1], eps);
        baby_log = log_from_cpct(G, baby_div);
        new_vec = gd_log - baby_log;
        if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec[1..r]~,lattice_lambda,eps)==0,

            if(DEBUG_BSGS>1,
                lattice_lambda_copy = lattice_lambda;
                verify_lattice_containment(~G, ~new_vec);
            );

            lattice_lambda = my_mlll( matconcat([lattice_lambda, new_vec[1..r]~]),eps);
            check_lambda_is_valid(~lattice_lambda);
        );
    );
    change_precision(old_precision);
    return(lattice_lambda);
}

\\ checks collisions by computing the log from a compact representation list
\\ but assumes the babystock as sufficiently accurate log approximations
collision_check1(G, ~lattice_lambda, r, eps,\
                 ~matches, ~giant_divisor, ~compactTracking)=
{
    for(i=1, length(matches),
        \\gd_log = log_from_cpct(G, giant_divisor[3])[1..r];
        gd_log = trackerLogarithm(G, ~compactTracking, r);
        new_vec = (gd_log - matches[i])[1..r];

        if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec~,lattice_lambda,eps)==0,

            if(DEBUG_BSGS>1,
                lattice_lambda_copy = lattice_lambda;
                verify_lattice_containment(~G, ~new_vec);
            );

            lattice_lambda = my_mlll( matconcat([lattice_lambda, new_vec~]),eps);
            if (matsize(lattice_lambda)[1] != matsize(lattice_lambda)[2],
                print("ERROR: mlll failed to properly determine new basis in giant steps.");
                breakpoint();
            );
            if(DEBUG_BSGS>0, print("new regulator = ", precision(matdet(lattice_lambda),10)););
        );
    );
    return(lattice_lambda);
}

\\ recomputes the babystock logarithm
\\ matches is the list of baby stock log vectors which share the matchIdeal
\\ compactTracking is a list of compact representations, which represent the giant divisor element
collision_check2(G, ~lattice_lambda, r, eps,\
                 matchIdeal, ~matches, ~giant_divisor, ~compactTracking)=
{
    my(
        deg = poldegree(G.pol),
        Ok = matid(deg),
        mainbitprecision = default(realbitprecision)
    );
    /*
        gd_precision = prec_compact(deg, log(abs(G.disc)), infinity_norm("giant divisor's log"); breakpoint();
        "Get a loose approximation of the log vector using the cpct representation, then determine gd_precision"
        \\ Decide if we can actually lower precision here. We might not
        default(realbitprecision, gd_precision);    \\ reduce precision
        \\\ ... do some computation
        default(realbitprecision, mainbitprecision); \\ restore original precision
        \\ s_radius = bitprecision(s_radius, mainbitprecision);
    */

    gd_log = trackerLogarithm(G, ~compactTracking, r);
    \\#bound based on the shortest vector of the ideal lattice (for a reduced ideal)
    compare_precision = 2^(-ceil((1/poldegree(G.pol))*log(abs(G.disc))+3/2)+1) ;
    for(i=1, length(matches),

        print("computing matched babystock element");
        baby_element = compact_rep_full_input(G, matches[i], matchIdeal, compare_precision,1, 2);
        print("done computing");
        psi_baby = log_from_cpct(G, baby_element);
        if(norml2(matches[i] - psi_baby) > 0.00000001,
            print("CC2 mismatch: ", precision(matches[i],20), "   ", precision(psi_baby,20) );
        );


        new_vec = gd_log - psi_baby;

        if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec[1..r]~,lattice_lambda,eps)==0,

            if(DEBUG_BSGS>1,
                lattice_lambda_copy = lattice_lambda;
                verify_lattice_containment(~G, ~new_vec);
            );
            lattice_lambda = my_mlll( matconcat([lattice_lambda, new_vec[1..r]~]),eps);
            if (matsize(lattice_lambda)[1] != matsize(lattice_lambda)[2],
                print("ERROR: mlll failed to properly determine new basis in giant steps.");
                breakpoint();
            );
            if(DEBUG_BSGS>0, print("new regulator = ", precision(matdet(lattice_lambda),10)););
        );
    );
    return(lattice_lambda);
}


\\ Subalgorithm of bsgs. computes giant steps incrementally using ideal
\\ multiplication to obtain adjacent elements.
\\ Compact Represenation Version
incremental_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, ~babystock, avec,\
                        eps, outFileInfo=[])=
{
    my(timeout, OUTFILE_GS,
        localbitprecision = 50;\\prec_compact(field_deg, log(abs(G.disc)), log(normlp(max_element)));
    );
    \\print(REQ_JUMP(G, column_sum(giant_sublattice)));
    if(length(outFileInfo) == 2,
        timeout = outFileInfo[2];
        OUTFILE_GS = outFileInfo[1];
    ,
        timeout = 0;
        OUTFILE_GS = 0;
    );

    my(
        LOG_REFRESH_RATE = 500,
        field_deg = poldegree(G.pol),   r = G.r1 + G.r2 -1,
        zero_vec = vector(r, i , 0),    identity = matid(field_deg),
        giant_coeffs = zero_vec,
        giant_divisor,                  expected_position,
        matches,
        total_steps,                    logging_interval,
        place_marker = r,               directions
    );
    GP_ASSERT_EQ(r, length(giant_sublattice));


    total_steps = 1; for(i =1, length(avec), total_steps*=avec[i]);
    logging_interval = max(floor(total_steps/20),1);
    print("Total number of giantsteps needed: ", total_steps);

    expected_position = vector(G.r1+G.r2, i, 0)~;
    directions = vector(length(giant_sublattice), i , 1);
    giant_coeffs[r] = 1;    \\# initialize variable to track giant elements at 1

    my(
        \\# vectors of length r which are used to compute an adjacent element in
        \\# a particular direction in the logarithm lattice
        \\# direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
        direction_elements,                 \\# set of corresponding to a giant step in each lattice dimension
        inverse_direction_elements,         \\# inverse elements of the above
        compactTracking = List(),           \\# list of cpct reps corresponding to current giant step
        trackingLog = vector(r+1, i, 0),    \\# log corresponding to current giant step. Precision drops fast
        s_radius                            \\# max acceptable distance from expected position in log lattice
    );

    mainbitprecision = default(realbitprecision);
    default(realbitprecision, localbitprecision);
    s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    max_element = lattice_lambda[,1];
    for(i=2, length(lattice_lambda), max_element += lattice_lambda[,i]);
    default(realbitprecision, mainbitprecision);
    s_radius = bitprecision(s_radius, mainbitprecision);

    \\# for some reason I have to flip these. giantstep() may actually return an inverse
    \\# initialize direction elements and inverses
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, giant_sublattice, field_deg, eps);

    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );

    emptyCompactList = compactTracking;
    giant_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Giant Steps LOG: With additional timing variables");
    my(giant_t1, giant_tn, giant_tmid, ctr,
        tDivisorCompute, tDivisorReduce,
        tCompare, tNext);
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    tDivisorCompute = 0;
    tDivisorReduce = 0;
    tCompare = 0;
    tNext = 0;
    \\ #increments are done modulo the product of avec, returns to zero when
    \\ #all elements have been visited
    while(giant_coeffs != zero_vec,
        tDC = getabstime();
        if(directions[place_marker] == 1,
            giant_divisor = [idealmul(G, giant_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],direction_elements[place_marker][2] )~];

            compactTracking[place_marker][2] += 1;
            default(realbitprecision, localbitprecision);
            expected_position += giant_sublattice[,place_marker];
            trackingLog += log_from_cpct(G, direction_elements[place_marker][3]); \\ MARK_OPTIMIZE
            default(realbitprecision, mainbitprecision);
        ,\\else
            giant_divisor = [idealmul(G, giant_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],inverse_direction_elements[place_marker][2] )~];

            compactTracking[place_marker][2] -= 1;
            default(realbitprecision, localbitprecision);
            expected_position -= giant_sublattice[,place_marker];
            trackingLog -= log_from_cpct(G, direction_elements[place_marker][3]); \\ MARK_OPTIMIZE
            default(realbitprecision, mainbitprecision);
        );

        bitprecision(expected_position, mainbitprecision);
        tDR = getabstime();
        tDivisorCompute += (tDR-tDC);

        my(wasModified);
        wasModified = get_next_giant_divisor_cpct(~G, ~giant_divisor, ~compactTracking);
        t_half = getabstime();
        tNext += (t_half-tDR);

        default(realbitprecision, localbitprecision);
        if(wasModified,
            trackingLog += embeddings_to_normalized_logvec(G, nfeltembed(G, compactTracking[length(compactTracking)][1]));
        );

        default(realbitprecision, mainbitprecision);
        trackingLog = bitprecision(trackingLog, mainbitprecision);

        [giant_divisor, tempLog] = adjust_giant_step_cpct(~G, ~giant_divisor, ~compactTracking, bitprecision(trackingLog,mainbitprecision), ~expected_position, bitprecision(s_radius, mainbitprecision), eps);

        tAdjust = getabstime();
        tDivisorReduce += (tAdjust - t_half);

        EXACT_CHECK = 1;
        \\ use babystock hashmap to check for collisions and update as needed
        if(mapisdefined(babystock, giant_divisor[1]),
            \\if(DEBUG_BSGS>0 ||1 ,print(ctr,"  baby-giant match. checking for new vector: ", giant_divisor[1]););
            matches = mapget(babystock, giant_divisor[1]);

            \\#use these as a template to debug giant steps if needed.
            \\if (giant_divisor[1] == [1, 0, 0, 0, 6/43; 0, 1, 0, 0, 128/129; 0, 0, 1, 0, 91/129; 0, 0, 0, 1, 89/129; 0, 0, 0, 0, 1/129],
            \\    print("gdbad ", precision(giant_divisor[2],50) );
            \\    bstock_vec = mapget(babystock, giant_divisor[1]);print("bstock ", precision(bstock_vec,50));breakpoint(););
            \\if (giant_divisor[1] == [1, 0, 0, 7804915/15049807; 0, 1, 0, 8863997/15049807; 0, 0, 1, 1043039/15049807; 0, 0, 0, 1/15049807],
            \\    print("gdbad ", precision(giant_divisor[2],50) );print("bstock ", precision(mapget(babystock, giant_divisor[1]),50));breakpoint(););

            if (!EXACT_CHECK,
                lattice_lambda = collision_check1(G, ~lattice_lambda, r, eps,
                                 ~matches, ~giant_divisor, ~compactTracking);
            ,
                 lattice_lambda = collision_check_bootstrap_log(G, ~lattice_lambda, r,
                                  ~matches, ~giant_divisor, trackingLog);
            );
        );

        \\ refresh the logarithm using compact reps to eliminate precision losses
        if ((ctr%LOG_REFRESH_RATE == 0) && ctr > 0,
            base_value = compact_rep_full_input(G, trackingLog, giant_divisor[1], eps,1,2);
            trackingLog = log_from_cpct(G, base_value);
        );
        tEnd = getabstime();
        tCompare += (tEnd - tAdjust);

        \\# increase the tracking variable and update directions
        place_marker = increment_with_place_marker(~avec, ~giant_coeffs);
        updateDirections(~directions, ~place_marker);

        if((timeout > 0)&&(tEnd - giant_t1 > timeout),
            write(OUTFILE_GS, "gstep computation ", precision((tEnd - giant_t1)/60000.0,10), " mins exceeds timeout. Halting gsteps");return(lattice_lambda);
        );
        ctr++;

        if(ctr %logging_interval ==0,
            print("Roughly ", 5*floor(ctr/logging_interval),"% complete. (", ctr, "/",total_steps, ")."  );
            print("Time breakdown: Compute: ", tDivisorCompute, " Reduce: ", tNext," Adjust: ", tDivisorReduce, " Compare: ", tCompare);
            giant_tn = getabstime();
            giant_tmid = giant_tn;
        );
    );
    print("Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);

    return(lattice_lambda);
}

\\ Subalgorithm of bsgs. computes giant steps incrementally using ideal
\\ multiplication to obtain adjacent elements.
\\ Compact Represenation Version
incremental_giant_steps_compact(~G, ~lattice_lambda, ~giant_sublattice, ~babystock, avec,\
                        eps, storage, cpct_bstock, outFileInfo=[])=
{

    if(storage != "COMPACT" && storage != "LOG", print("invalid storage format"); break;);
    if(storage == "LOG", print("Warning there is a bug in this setting, use the regular incremental giant steps instead"));
    my(timeout, OUTFILE_GS);
    if(length(outFileInfo) == 2,
        timeout = outFileInfo[2];
        OUTFILE_GS = outFileInfo[1];
    ,
        timeout = 0;
        OUTFILE_GS = 0;
    );

    my(
        field_deg = poldegree(G.pol),   r = G.r1 + G.r2 -1,
        zero_vec = vector(r, i , 0),    identity = matid(field_deg),
        giant_coeffs = zero_vec,
        giant_divisor,                  expected_position,
        matches,
        place_marker = r,               directions
    );
    GP_ASSERT_EQ(r, length(giant_sublattice));

    expected_position = vector(G.r1+G.r2, i, 0)~;       \\# r1+r2 length vec
    directions = vector(length(giant_sublattice), i , 1);

    \\# initialize variable to track giant elements at 1
    giant_coeffs[r] = 1;

    my(
        \\# vectors of length r which are used to compute an adjacent element in
        \\# a particular direction in the logarithm lattice
        \\# direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
        direction_elements,                 \\# set of corresponding to a giant step in each lattice dimension
        inverse_direction_elements,         \\# inverse elements of the above
        compactTracking = List(),           \\# list of cpct reps corresponding to current giant step
        trackingLog = vector(r+1, i, 0),    \\# log corresponding to current giant step. Precision drops fast
        s_radius,                            \\# max acceptable distance from expected position in log lattice
        one_element,
        base_value = [[1], [1]],
        max_element,
        localbitprecision = 30
    );

    mainbitprecision = default(realbitprecision);
    default(realbitprecision, localbitprecision);
    s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    max_element = column_sum(lattice_lambda);
    max2 = lattice_lambda[,1];
    for(i=2, length(lattice_lambda), max2 += lattice_lambda[,i]);
    GP_ASSERT_NEAR(norml2(max_element - max2), eps);
    default(realbitprecision, mainbitprecision);
    s_radius = bitprecision(s_radius, mainbitprecision);

    \\# initialize direction elements and inverses
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, giant_sublattice, field_deg, eps);

    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );

    emptyCompactList = compactTracking;

    giant_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Additional timing variables present");
    my(giant_t1, giant_tn, giant_tmid, ctr,
        tDivisorCompute, tDivisorReduce,
        tCompare, tNext
    );
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    tDivisorCompute = 0;
    tDivisorReduce = 0;
    tCompare = 0;
    tNext = 0;

    localbitprecision = REQ_JUMP(G, max_element);
    print("GS precision based on jump to largest element: ", localbitprecision);
    \\ #increments are done modulo the product of avec, returns to zero when
    \\ #all elements have been visited
    while(giant_coeffs != zero_vec,

        tDC = getabstime();
        if(directions[place_marker] == 1,
            giant_divisor = [idealmul(G, giant_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],direction_elements[place_marker][2] )~];
            if(storage == "COMPACT",
                compactTracking[place_marker][2] += 1;
                expected_position += giant_sublattice[,place_marker];
            );
            if(storage == "LOG",
                default(realbitprecision, localbitprecision);
                expected_position += giant_sublattice[,place_marker];
                trackingLog += log_from_cpct(G, direction_elements[place_marker][3]);
                default(realbitprecision, mainbitprecision);
            );

        ,\\else
            giant_divisor = [idealmul(G, giant_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],inverse_direction_elements[place_marker][2] )~];

            if(storage == "COMPACT",
                compactTracking[place_marker][2] -= 1;
                expected_position -= giant_sublattice[,place_marker];
            );
            if(storage == "LOG",
                default(realbitprecision, localbitprecision);
                expected_position -= giant_sublattice[,place_marker];
                trackingLog -= log_from_cpct(G, direction_elements[place_marker][3]);
                default(realbitprecision, mainbitprecision);
            );
        );

        bitprecision(expected_position, mainbitprecision);
        tDR = getabstime();
        tDivisorCompute += (tDR-tDC);
        get_next_giant_divisor_cpct(~G, ~giant_divisor, ~compactTracking);
        t_half = getabstime();
        tNext += (t_half-tDR);

        if(storage == "LOG",
            default(realbitprecision, localbitprecision);
            trackingLog += embeddings_to_normalized_logvec(G, nfeltembed(G, compactTracking[length(compactTracking)][1]));
            default(realbitprecision, mainbitprecision);
            trackingLog = bitprecision(trackingLog, mainbitprecision);
            if (ctr > 0 && (ctr%200 == 0),
                \\print(precision(trackingLog,10));
                default(realbitprecision, localbitprecision);
                renewElement = compact_rep_full_input(G, trackingLog, identity, eps);
                renewLog = log_from_cpct(G, renewElement);
                default(realbitprecision, mainbitprecision);
                \\print("renewed log: ",precision(renewLog,10)); breakpoint();
            );
        );
        if(storage == "COMPACT",
            \\trackingLog = vector(r+1, i, 0);
            trackingLog = log_from_cpct(G, base_value);
        );

        [giant_divisor, tempLog] = adjust_giant_step_cpct(~G, ~giant_divisor, ~compactTracking, bitprecision(trackingLog,mainbitprecision), ~expected_position, bitprecision(s_radius, mainbitprecision), eps,storage);

        tAdjust = getabstime();
        tDivisorReduce += (tAdjust - t_half);

        EXACT_CHECK = 1;
        \\ use babystock hashmap to check for collisions and update as needed
        if(mapisdefined(babystock, giant_divisor[1]),
            if(0 ,print("baby-giant match. checking for new vector:", giant_divisor[1]););

            matches = mapget(babystock, giant_divisor[1]);
            if (storage == "COMPACT",
                matches_cpct = mapget(cpct_bstock, giant_divisor[1]);
            );
            \\#if (giant_divisor[1] == [1, 0, 0, 7804915/15049807; 0, 1, 0, 8863997/15049807; 0, 0, 1, 1043039/15049807; 0, 0, 0, 1/15049807],
            \\#    print("gdbad ", precision(giant_divisor[2],50) );
            \\#    print("bstock ", precision(mapget(babystock, giant_divisor[1]),50));
            \\#    for(i=2, length(matches_cpct),
            \\#        print(precision(log_from_cpct(G,~matches_cpct[i][1])+trackerLogarithm(G, ~matches_cpct[i][2], r),50));
            \\#    );
            \\#    breakpoint();
            \\#);
            if (storage == "LOG",
                if (!EXACT_CHECK,
                    lattice_lambda = collision_check1(G, ~lattice_lambda, r, eps,
                                     ~matches, ~giant_divisor, ~compactTracking);
                ,
                    lattice_lambda = collision_check2(G, ~lattice_lambda, r, eps,
                                    giant_divisor[1], ~matches, ~giant_divisor, ~compactTracking);
                );
            );
            if (storage == "COMPACT",
                \\\ compute the log assuming modified strategy
                giant_step_element_log = log_from_cpct(G, base_value)+trackerLogarithm(G, ~compactTracking, r);
                lattice_lambda = collision_check_bootstrap(G, ~lattice_lambda, r,
                                 ~matches_cpct, ~giant_divisor, ~giant_step_element_log);
                \\lattice_lambda = collision_check_cpct(G, ~lattice_lambda, r, eps,
                \\                 ~matches_cpct, ~giant_divisor, ~compactTracking);
            );
        );
        if(ctr < 50, print(giant_coeffs, "   ", precision(giant_divisor[1],10)));
        tEnd = getabstime();
        tCompare += (tEnd - tAdjust);

        \\# increase the tracking variable and update directions
        place_marker = increment_with_place_marker(~avec, ~giant_coeffs);
        updateDirections(~directions, ~place_marker);

        if((timeout > 0)&&(tEnd - giant_t1 > timeout),
            write(OUTFILE_GS, "gstep computation ", precision((tEnd - giant_t1)/60000.0,10), " mins exceeds timeout. Halting gsteps");return(lattice_lambda);
        );
        ctr++;
        if(ctr %50 ==0,
            GP_ASSERT_NEAR(norml2(trackerLogarithm(G, ~compactTracking, r)+log_from_cpct(G, base_value)- tempLog), 0, 2^(-40));
            base_value = compact_rep_full_input(G, tempLog, matid(field_deg), eps);
            compactTracking = emptyCompactList;

            if(ctr %500 == 0,
                print(ctr, " giant steps taken.");
                print("Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);
                giant_tn = getabstime();
                \\write("data/jump-timing.txt", ctr, " jumps took: ", (giant_tn- giant_tmid), "milliseconds");
                giant_tmid = giant_tn;
            )
        );
    );
    print("Giant steps Completed: Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);

    return(lattice_lambda);
}

\\ gets a nearby reduced divisor to log_coordinates_Rn
\\ INPUTS:
\\ - G a number field, idealmat a coefficient matrix of an ideal
\\ - log_coordinates_Rn is the log vector (length = equal to unit rank)
\\ - inverse is an optional argument that will compute the inverse instead if set to 1
get_nearby_rdivisor(G, idealmat, log_coordinates_Rn, inverse=1)={
    exponentiated_coordinates = exponentiate_logvec(G.r1+G.r2-1, unsquare_log_embeddings(G, log_coordinates_Rn), inverse);
    \\exponentiated_coordinates = create_target(G, log_coordinates_Rn, inverse);

    return( reddiv_compact(idealmat,exponentiated_coordinates, G, G[5][1] ); );
}

multiply_and_reduce_divisors(G, ideal1, real_vec1, log_vec1, ideal2, real_vec2, log_vec2)={
    product_ideal = idealmul(G, ideal1, ideal2);
    product_real = pointwise_vector_mul(real_vec1,real_vec2);
    reduced_product = reddiv_compact(product_ideal, product_real~, G, G[5][1]);

    reduced_product_log = log_vec1 + log_vec2 + reduced_product[3];

    return([reduced_product[1], reduced_product[2], reduced_product_log]);
}
\\ INPUT:
\\ - G a number field
\\ - giant_divisor is a 3-element list consisting of
\\      an ideal, a vector of size G.r1+G.r1 corresponding to valuation of an element
\\      a vector of size G.r1+G.r1 corresponding to log distance of the divisor
\\ - expected position is the coordinates in R^r of your giant lattice element
\\ - distance_ok is the limit of the acceptable distance from the divisor to the target position
adjust_giant_step(G, giant_divisor, expected_position, distance_ok,eps)={
    my(
        r = G.r1 + G.r2 -1,
        divisor_distance = expected_position - giant_divisor[3][1..r]~,
        adjustment_divisor, new_divisor,
        new_distance
    );
    \\ if the distance is sufficiently close then no need to do anything
    if(norml2(divisor_distance) < distance_ok,
        return(giant_divisor);
    , \\else
        \\#print("\ntarget position = ", precision(expected_position,10), "\nOriginal Distance from Target ", precision(norml2(divisor_distance),10) );
        for(i=1, 2,
            adjustment_divisor = get_nearby_rdivisor(G, matid(poldegree(G.pol)), divisor_distance, i%2);
            if (norml2(adjustment_divisor[2]) < eps,
                return(giant_divisor);
            ,
                new_divisor = multiply_and_reduce_divisors(G, giant_divisor[1], giant_divisor[2], giant_divisor[3], adjustment_divisor[1], adjustment_divisor[2], adjustment_divisor[3]);
                new_distance = norml2(expected_position - (new_divisor[3][1..r]~));
                \\#print("Adjusted distance from target = ", precision(new_distance, 10) );
                if(new_distance < norml2(divisor_distance)+eps,
                    return(new_divisor);
                );
            );
        );
        return(giant_divisor);
    );
};


get_giant_step_increment_vectors(G, giant_sublattice, field_deg, eps)={
    my( giant_step_ideals = [], giant_step_ideals_inverse = [],
        new_vec,gstep_divisor);

    for(j = 1, length(giant_sublattice),
        new_vec = giant_sublattice[,j];
        gstep_divisor = giantstep(matid(field_deg), new_vec, G, field_deg, eps);
        \\debug_compare(gstep_divisor, giantstep_high_precision(matid(field_deg), new_vec, G, field_deg, eps));

        giant_step_ideals = concat(giant_step_ideals, [[gstep_divisor[1], gstep_divisor[2], gstep_divisor[3]]]);
        giant_step_ideals_inverse = concat(giant_step_ideals_inverse, [[idealinv(G,gstep_divisor[1]), invert_coordinates(gstep_divisor[2]),-gstep_divisor[3]  ]] );
    );

    return([giant_step_ideals,giant_step_ideals_inverse]);
}

get_next_giant_divisor(G, giant_divisor)={
    my(new_divisor);
    new_divisor = reddiv_compact(~giant_divisor[1], ~giant_divisor[2], ~G,~G[5][1]);
    if(type(new_divisor[3])== "t_COL", new_divisor[3] = new_divisor[3]~);
    return([new_divisor[1], new_divisor[2], giant_divisor[3]+new_divisor[3]]);
}


get_next_giant_divisor_cpct(~G, ~giant_divisor, ~tracker)={
    \\ MARK_OPTIMIZE if new_divisor was a reference input, would not need to create/destroy it
    my(new_divisor, oneVector = vector(poldegree(G.pol), i, 0)~);

    oneVector[1] = 1;
    new_divisor = reddiv_compact(~giant_divisor[1], ~giant_divisor[2], ~G,~G[5][1]);

    giant_divisor[1] = new_divisor[1];
    giant_divisor[2] = new_divisor[2];

    \\# if the coefficient vector corresponds to the 1st elementary vector
    \\# it's just 1, so don't store it.
    if(length(tracker) > poldegree(G.pol) && new_divisor[4] == oneVector,
        \\# do nothing
        return(0);
    ,
        listput(~tracker, [new_divisor[4], 1]);
        return(1);
    );
}

compute_default_scan_radius(G)=
{
    return ((sqrt(poldegree(G.pol))/4)*log(abs(G.disc)));
}
\\ INPUT:
\\ - G is a number field
\\ - cpct_reps is a list of cpct representations corresponding to a sublattice of the unit lattice
\\ - B is an integer that indicates a fraction of the fundamental region to search
\\   Divides the last coordinate vector by B.
\\ - babystock_scale_factor indicates the volume of the babystock region
\\ - scanballRadius indicates how large each of the ball regions is in scanball
\\ - eps is the error, REQ_BSGS is the calculated precision needed for accurate results
\\ - REQ_BSGS,FILE1 is a file for certain outputs
\\ - options is a vector of additional arguments.
\\ - options[1] may specify a different algorithm to scan with (not available yet)
\\ - or it specifies a timeout value, which will be checked periodically and
\\ -    the program terminated if exceeded
bsgs(G, cpct_reps, B, babystock_scale_factor, scanballRadius,eps, REQ_BSGS,FILE1="data/tmp-bsgs-log", options = [])={

    bsgs_start = getabstime();
    my(
        r = G.r1 + G.r2 -1,     field_deg = poldegree(G.pol),
        zero_vec = vector(r, i , 0),
        babystock = Map(),
        timeout = 0,
        alg = "SCAN",
        storage = "LOG"     \\storage variable has 1 valid option "LOG". "COMPACT" is now
    );
    if (storage == "COMPACT", cpct_bstock = Map());
    print("\nBSGS Algorithm Start. Storage type: ", storage);
    \\ check if an auxilliary string for output has been provided
    if (length(options)!=0,
        if(type(options[1]) == "t_STR",
            alg = options[1];,
            timeout = options[1];
        );
    );

    if(alg != "NEIGHBOURS" && alg != "SCAN", print("Unknown algorithm specified, returning -1", return(-1) ));
    if(DEBUG_BSGS>0 , print("Bound B = ", precision(B,10)); );

    my(lattice_lambda, S_radius, avec, giant_sublattice);
    lattice_lambda = log_lattice_from_compact_set(G, cpct_reps);
    S_radius = compute_default_scan_radius(G);
    [avec, giant_sublattice] = get_giant_step_params(G,lattice_lambda, r, B, babystock_scale_factor, REQ_BSGS);

    tb = getabstime();

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #THIS OPTION USES NEIGHBORS TO SEARCH THE BABYSTOCK
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(alg != "SCAN",
        \\# note this option has been disabled for some time so it no longer works
        print("SCAN is the only strategy implemented. Rerun using SCAN"); quit;
        neighbour_search(G,lattice_lambda, babystock_region, ~babystock, giant_sublattice, eps);
    );  \\# end neighbors for babystock

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #THIS SECTION SEARCHES USING THE SCAN ALGORITHM
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    minimal_vec_length = get_minimal_length(giant_sublattice);

    aprod = 1; for(i=1, length(avec), aprod*=avec[i]);
    if(minimal_vec_length < scanballRadius,
        scanballRadius = max(minimal_vec_length, 0.01);

        write(FILE1,strprintf("Adjusted scan radius: %-10.9F  Aproduct: %-10.9F", scanballRadius, aprod),"  aVec =", precision(avec, 10));
        ,
        write(FILE1,strprintf("Scan radius: %-10.9F  Aproduct: %-10.9F", scanballRadius, aprod),"  aVec =", precision(avec, 10));
    );

    if(alg == "SCAN",
        my(foundflag = 1, num_elts_found = []);

        if (storage == "COMPACT",
            print("Deprecated storage format");breakpoint();
            [lattice_lambda, num_elts_found] =
                incremental_baby_steps_compact(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, storage, ~cpct_bstock, [FILE1, timeout]);
        ,
            [lattice_lambda, num_elts_found] = b_scan(G, matid(field_deg),~lattice_lambda, ~giant_sublattice,
                ~babystock, scanballRadius, eps, [FILE1, timeout]);
            \\[lattice_lambda, num_elts_found] = b_scan_jump(G, matid(field_deg),~lattice_lambda, ~giant_sublattice,~babystock, scanballRadius, eps, [FILE1, timeout]);
        );

        foundflag = length(num_elts_found);
        print("babystock size: ", length(babystock));
        while(foundflag !=0,
            \\# if new elements found, regenerate lattice_lambda to remove precision loss
            print("Readjusting precision due to new element");
            \\newprec = REQ_GIANT(G, get_abs_determinant(lattice_lambda), column_sum(lattice_lambda));
            newprec = prec_bsgs(field_deg, log(abs(G.disc)), sum_inf_norm(lattice_lambda));
            \\#print(newprec1, "  ", newprec); breakpoint();
            REQ_BSGS = newprec; default(realbitprecision, REQ_BSGS);

            lattice_lambda = log_lattice_from_compact_set(G, cpct_from_loglattice(G, lattice_lambda, eps));
            [avec, giant_sublattice] = get_giant_step_params(G,lattice_lambda, r, B, babystock_scale_factor, REQ_BSGS);

            \\# then rerun the baby steps using the new lattice
            if (storage == "COMPACT",
                [lattice_lambda, num_elts_found] =
                    incremental_baby_steps_compact(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, storage, cpct_bstock, [FILE1, timeout]);
            ,
                    \\[lattice_lambda, num_elts_found] = b_scan_jump(G, matid(field_deg),~lattice_lambda, ~giant_sublattice,~babystock, scanballRadius, eps, [FILE1, timeout]);
                    [lattice_lambda, num_elts_found] = b_scan(G, matid(field_deg),~lattice_lambda, ~giant_sublattice,
                        ~babystock, scanballRadius, eps, [FILE1, timeout]);
            );

            GP_ASSERT_EQ(matsize(lattice_lambda)[1], matsize(lattice_lambda)[2]);
            foundflag = length(num_elts_found);
            if(DEBUG_BSGS>0 ,print("FoundFlag ",foundflag););
        );
    );
    if(0,
        baby_matrix = Mat(babystock);
        print(precision(baby_matrix,20));
        for(debug_ctr=1, matsize(Mat(babystock))[1],
            print(matsolve(giant_sublattice, mapget(babystock,baby_matrix[debug_ctr,1])[1]~))
        );
    );

    babytime = getabstime()-tb;
    write(FILE1, "Babystock time: ", precision(babytime, 10), "   ",  precision(babytime/60000.0, 15), ". Babystock elements: ", length(babystock));

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #Begin giant step computations
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    print("---------------\nGiant Loop\n" );
    my(aProduct = 1);
    for(i=1, length(avec), avec[i] += 1; aProduct*=avec[i];);                                       \\ adjust avec to include boundaries.
    avec = round(avec);
    print("-- a-vector prior to GS: ", precision(avec,10), " Lambda dimensions: ", matsize(lattice_lambda));
    \\# the search region area is the product of the norms of each sublattice vector
    \\#normproduct = 1; for(i=1, length(giant_sublattice), normproduct*=sqrt(norml2(giant_sublattice[,i])) ); print("region size: ", normproduct);
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    tg = getabstime();
    print("eps giant step start: ", eps);
    if(storage == "LOG",
        lattice_lambda =
            incremental_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, \
                                    ~babystock, avec, eps, [FILE1, timeout]);
    );
    if(storage == "COMPACT",
        lattice_lambda =
            incremental_giant_steps_compact(~G, ~lattice_lambda, ~giant_sublattice, \
                                    ~babystock, avec, eps, storage, cpct_bstock, [FILE1, timeout]);
    );
    \\#lattice_lambda = jump_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, ~babystock, avec, eps);

    bsgs_end = getabstime(); gianttime = bsgs_end-tg;

    write(FILE1,  "Giantstep time:  ", gianttime,  "   ",precision(gianttime/60000.0,15), " . Giant steps computed: ", precision(aProduct,10));

    return(cpct_from_loglattice(G,lattice_lambda,eps));
}
