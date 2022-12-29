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
get_subdivisions
increment_coordinates
increment_with_place_marker
get_axis_aligned_box
get_giant_step_params
check_babystock_for_units
scanball
update_directions
babystock_scan

jump_giant_steps
get_nearby_rdivisor
multiply_and_reduce_divisors
adjust_giant_step
is_repeat_babystock
expand_babystock_region
get_initial_neighbour
get_giant_step_increment_vectors
get_next_giant_divisor
one_massive_qfminim

bsgs
*/
DEBUG_BSGS = 1;
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
/*28. Find a new vector v \in \Lambda_K\Lambda' where (O_K, v) \in \ep_B: norml2, fun. in 13 (is_vec_in_lattice),  fun. in 11 (updatelamda)*/
\\ Check the set ball_minima for a new vector in Lambda not already in Lambda'
\\ INPUT:
\\ - ball_minima a list of minima, each of the form [ideal, log vector] = [1/mu*O_K, vector = log(mu)]
\\ - L is the lattice we are checking (Lambda')
\\ - G the number field
\\ - eps the error

/******************************************************************************/
check_babystock_for_units(ball_minima,L,G,eps)={

    my(new_counter:small, candidate, ideal_identity, eps_sqr = eps^2, n:small = poldegree(G.pol));
    ideal_identity = matid(n);
    new_counter=0;                                                              \\# track the number of new minima found.

    for(i=1,length(ball_minima),
        if(ball_minima[i][1] == ideal_identity,                                 \\# if the ideal (1/mu)*O_k = O_k then check further
            candidate=ball_minima[i][2];                                                \\ candidate = log(mu)
            if(norml2(candidate)>eps_sqr&&is_vec_in_lattice(candidate~,L,eps_sqr)==0,   \\ if nonzero and v is not already in L,
                new_counter+=1;
                print("Babystock unit found, " precision(L,10), "  ", precision(candidate,10));
                L = my_mlll(matconcat([L, candidate~]), eps);
            )
        )
    );
    return([L,new_counter]);        \\ modifed to return new_counter as well
}

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
\\ - a Set of minima in B_D, each represented as a pair (ideal, logvector)
/******************************************************************************/
scanball(~G, y,u,psimu,web,eps)={
    my(
        n = poldegree(G.pol),
        ball_minima=Set(),
        x, scan_bound,
        vecholder, gram_mat,
        scan_elements,
        vec_numerical,
        LLL_reduced_yu,
        new_y,
        scanball_map = Map(), repeat_ctr = 0
    );

    x = G[5][1]*y;                                                              \\ numerical representation of y (complex)
    x = mulvec(x,u);                                                            \\ compute y*u
    x = embed_real(G,x);
    LLL_reduced_yu = x*qflll(x);                                                \\ lll reduce y*u
    vecholder = LLL_reduced_yu[,1];                                             \\ short vector, 1st element of LLL basis
    scan_bound = sqrt(n)*exp(2*web)*sqrt(norml2(vecholder));                    \\ See schoof alg 10.7, e^(2*var_eps)*sqrt(n)*sqrt(norml2(col))
    gram_mat=LLL_reduced_yu~*LLL_reduced_yu;                                    \\ get the gram matrix

    scan_elements=qfminim(gram_mat,scan_bound^2,,2)[3];
    scan_elements=y*scan_elements;                                              \\ get scanned elements wrt integral basis

    my(norm_deltaK = ceil(((2/Pi)^(G.r2))*abs(G.disc)^(1/2)*idealnorm(G,y)));
    for(ii=1, length(scan_elements),
        eltnorm = abs(nfeltnorm(G,scan_elements[,ii] ));
        \\\ Easy necessary condition for minimum
        \\\ norm of a minimum should satisfy 1 < N(s_elt) < N(y)*delta_K
        if(eltnorm>=1 && eltnorm<=norm_deltaK,
            new_y=idealdiv(G,y,scan_elements[,ii]);                                                    \\ get the ideal (1/w_i)*y
            new_yLLL = embed_real(G, G[5][1]*new_y)*qflll(embed_real(G, G[5][1]*new_y));
            if(norml2(new_yLLL) > 1-eps,
                if(checkred_old(new_y,G,eps)==1,                                                  \\ check if its reduced, and if yes, add add to ball_minima
                    vec_numerical = (G[5][1]*scan_elements[,ii])~;                                            \\ get the numericals
                    \\print("SCANBALL: Psi of potential minima: ",precision(vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps), 10 ));
                    ball_minima = setunion(ball_minima, Set([[new_y,vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps)]] ));
                    if (mapisdefined(scanball_map, new_y),
                        repeat_ctr += 1;
                    ,
                        mapput(~scanball_map, new_y, vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps));
                    );
                );
            );
        );
    );
    \\print("scanball elts ", length(scan_elements),"  scanball minima: ", length(ball_minima));
    return(ball_minima);
}

\\ subroutine of baby-step giant-step
\\ based on the index of the left-most nonzero coordinate of change_vec, changed the values
\\ of the direction_vec mod 2. direction_vec is used to tell us whether we are moving in a positive/negative direction along the giant-step lattice
update_directions(change_vec , directions_vec)={
    for(i=1, length(change_vec),
        if(change_vec[i] != 0,
            for(j = i+1, length(change_vec),
                directions_vec[j] = (directions_vec[j] %2) +1;
            );
            i = length(change_vec)+1;
        );
    );
    return(directions_vec);
}

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


\\ this is Schoofs scan algorithm, modified so that if we find a new element of Lambda, we add it to Lambda' and recompute it
\\ Given an axis aligned box (b1, b2), where b1 = coords of the smallest corner, and b2 is the largest corner
\\ we divide the sides into m_i pieces so that the volume of each cell is 'small' .
\\ Schoofs algorithm suggests
\\ INPUT:
\\ - y is an ideal
\\ - L is the log lattice
\\ - babystock box are the two vectors which define the babystock region in L
\\ - G is the number field
\\ - eps is the error
babystock_scan(y,L,babystock_box,G,eps)={
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
        baby_hashmap = Map(),
        ideal_J, nu, logdist, beta,
        region_minima,
        directions,
        ball_minima,
        newctr = 0,
        field_deg = poldegree(G.pol),
        identity_n = matid(field_deg)
    );

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Establish a number of parameters including the web of regularly distributed points,
    \\ and some values that allow us to move easily between each of the points sucessively.
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    box_volume = vector( length(L), i, babystock_box[2][i]-babystock_box[1][i]);        \\ dimensions of the box
    box_subdivisions = ceil((field_deg)*box_volume);                                    \\ store m_i, the number of subdivisions on each side of the box

    web_increments = vector(length(L), i, box_volume[i]/box_subdivisions[i]);           \\ vector of distances to next web point in each direction
    \\ Contains two vectors of length r. Multiply by by web_step[1][i] to move forward in the ith direction, web_step[2][i] to move in reverse
    web_step = get_web_step_multipliers(G, length(L),web_increments);

    [webpoint, exp_webpoint, web_distance] = get_initial_webpoint(G, length(L), babystock_box[1], web_increments);

    if(DEBUG_BSGS>3,
        print("DEBUGGING Babystock 'web' parameters: ");
        print("Babystock subdivisions: ", box_subdivisions);
        print("max dist. betwn web points: ", precision(web_distance,10));
        print("1st web point in log lattice = ", precision(webpoint,10));
    );
    scanball_ctr = 0;
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Begin scanning the points in the web of regularly distributed points
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    [ideal_J, nu, logdist, beta] = reddiv_compact(y, exp_webpoint, G, G[5][1]);

    ball_minima = scanball(G, ideal_J, nu, logdist[1..length(L)], web_distance, eps);

    region_minima = Set(concat(region_minima,ball_minima));
    \\ check for new elements in Lambda

    scanIdeals = Map(); print("initializing scan ideal hash map");
    repeat_counter = 0;
    mapput(scanIdeals, ideal_J ,1);
    directions = vector(length(L), i , 1);                                              \\ initialize as all positive directions

    while(web_coords != zerovec || next_coords==[],                             \\ loops over all of the web points in the babystock region
        \\print(web_coords);
        increment_coordinates(box_subdivisions, web_coords);      \\ compute the next point's coefficients
        web_coords = web_coords - next_coords;                                  \\ obtain difference vector
        directions = update_directions(web_coords, directions);                 \\ based on the diff vector, change step directions

        for(k=1, length(L),
            if(web_coords[k] !=0,                                               \\ based on the largest place value that changed, update exp_webpoint and webpoint
                exp_webpoint = pointwise_vector_mul(exp_webpoint, web_step[directions[k]][k])~;
                if(directions[k] == 1, webpoint[k] += web_increments[k], webpoint[k] -= web_increments[k];);
                k = length(L)+1;
            );
        );
        \\print("coord loop ", next_coords, "  ", precision(exp_webpoint));
        web_coords = next_coords;
        \\output = giantstep(y, webpoint, G, field_deg, eps);
        \\[ideal_J, nu, logdist, beta] = reddiv_compact(ideal_J, exp_webpoint, G, G[5][1]);
        [ideal_J, nu, logdist, beta] = reddiv_compact(y, exp_webpoint, G, G[5][1]);

        if (mapisdefined(scanIdeals, ideal_J),
            print("Repeat ideal"); repeat_counter+=1;
            mapget(scanIdeals, ideal_J);
        ,
            mapput(scanIdeals, ideal_J ,1);
        );


        ball_minima = scanball(G, ideal_J, nu, logdist[1..length(L)], web_distance, eps);

        [L, newctr] = check_babystock_for_units(ball_minima, L, G, eps);
        if(newctr != 0,
            print("Found a unit in babystock. new reg =", precision(abs(matdet(L)),10) );
            babystock_box[1][1] += (web_coords[1]*web_increments[1]);
            return([L, baby_hashmap, babystock_box]);
        );
        scanball_ctr+=1;
        if(scanball_ctr%10 == 0, print(scanball_ctr, " scanballs performed"));
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\ Take the elements found on the ball near the current point in the web, and enter them into the hash map
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        region_minima = Set(concat(region_minima,ball_minima));
        for(ctr=1, length(ball_minima),
            if(mapisdefined(baby_hashmap, ball_minima[ctr][1]),
                existing_entry = mapget(baby_hashmap, ball_minima[ctr][1]);
                repeatflag = is_repeat_babystock(existing_entry, ball_minima[ctr][2],eps);

                if(repeatflag==0,
                    existing_entry = concat( mapget(baby_hashmap, ball_minima[ctr][1]), [ball_minima[ctr][2]]);
                    mapput(baby_hashmap, ball_minima[ctr][1], existing_entry );
                );
            ,\\else
                mapput(baby_hashmap, ball_minima[ctr][1], [ball_minima[ctr][2]]);
            );
        );
    );
    print("Babysteps stored"); print(repeat_counter);
    return([L,baby_hashmap, []]);
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
    \\print("babystock search time: ", getabstime() - time1);
    \\print("Number of repeated minima during scans: ", repeated_minima);
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
        field_deg = poldegree(G.pol),
        r = G.r1 + G.r2 -1,
        zero_vec = vector(r, i , 0),
        giant_coeffs = zero_vec,
        giant_divisor,
        current_giant_vec = zero_vec,
        matches, new_vec,
        identity = matid(field_deg)
    );
    giant_coeffs[r] = 1;
    SCREEN(0, "Additional timing variables and file write in jump");
    my(giant_t1, giant_tn, giant_tmid, ctr);
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    my(eps_sqr = (eps*10)^2);
    while(giant_coeffs != zero_vec,

        \\ Get the expected giant-step position, and then compute a nearby reduced divisor
        \\if(giant_coeffs[r] != 0,
        \\    current_giant_vec += gs_sublattice[,r]~;
        \\,
            current_giant_vec = (gs_sublattice*giant_coeffs~)~;
        \\);
        giant_divisor = giantstep(identity, current_giant_vec, G, field_deg, eps);

        if(mapisdefined(bstock, giant_divisor[1]),
            matches = mapget(bstock, giant_divisor[1]);                         \\ list of existing babystock matched elements
            for(i=1, length(matches),
                new_vec = giant_divisor[3][1..r] - matches[i];                  \\ compute difference
                if(norml2(new_vec) > eps_sqr && is_vec_in_lattice(new_vec~,lat_lambda,eps)==0,
                    print(matsize(lat_lambda));
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
    \\write("data/jump-timing.txt", "Total jump time: ", (giant_tn- giant_t1), "milliseconds\n");
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
    for(j = 1, length(giant_sublattice),
        new_vec = giant_sublattice[,j];
        gstep_divisor = jump_compact(matid(field_deg), new_vec, G, field_deg, eps);
        listput(~giant_step_ideals, gstep_divisor);
        \\debug_compare(idealdiv(G, matid(field_deg), compact_reconstruct(G,giant_step_ideals[j][3][1], giant_step_ideals[j][3][2])), giant_step_ideals[j][1]);
        GP_ASSERT_TRUE(sqrt(norml2(log_from_cpct(G, gstep_divisor[3])[1..urank]~-giant_sublattice[,j]))<log(abs(G.disc)^(1/2)) );
        gstep_divisor[3] = invert_compact(G, gstep_divisor[3]);
        listput(~giant_step_ideals_inverse, [idealinv(G,gstep_divisor[1]), invert_coordinates(gstep_divisor[2]),gstep_divisor[3] ] );

        myprod = mul_compact(G, giant_step_ideals[j][3], giant_step_ideals_inverse[j][3]);
        \\debug_compare(idealdiv(G, matid(field_deg), compact_reconstruct(G,giant_step_ideals_inverse[j][3][1], giant_step_ideals_inverse[j][3][2])), giant_step_ideals_inverse[j][1]);
        \\ #check that the cpct reps are actually inverses.
        \\ #also check that the reps are near the intended log lattice point
        print("WARNING: Verifying increment ideals");
        GP_ASSERT_EQ(compact_reconstruct(G, myprod[1], myprod[2]), 1);

    );
    return([giant_step_ideals, giant_step_ideals_inverse]);
}

updateDirections(~directions_vec, place_marker)=
{
    for(j = place_marker+1, length(directions_vec),
        directions_vec[j] = (directions_vec[j] %2) +1;
    );
}

collision_check_cpct(G, ~lattice_lambda, r, eps,\
                 ~matches_cpct, ~giant_divisor, ~compactTracking)=
{
    for(i=1, length(matches_cpct),
        gd_log = trackerLogarithm(G, ~compactTracking, r);
        baby_log = trackerLogarithm(G, ~matches_cpct[i], r);
        new_vec = gd_log - baby_log;
        print("match[i] ", precision(matches_cpct[i],10));
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

collision_check1(G, ~lattice_lambda, r, eps,\
                 ~matches, ~giant_divisor, ~compactTracking)=
{
    for(i=1, length(matches),
        \\gd_log = log_from_cpct(G, giant_divisor[3])[1..r];
        gd_log = trackerLogarithm(G, ~compactTracking, r);
        new_vec = gd_log - matches[i];
        \\print("match[i] ", precision(matches[i],10));
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

collision_check2(G, ~lattice_lambda, r, eps,\
                 matchIdeal, ~matches, ~giant_divisor, ~compactTracking)=
{
    my(
        deg = poldegree(G.pol),
        mainbitprecision = default(realbitprecision)
    );

    /*
        gd_precision = prec_compact(deg, log(abs(G.disc)), infinity_norm("giant divisor's log");
        "Get a loose approximation of the log vector using the cpct representation, then determine gd_precision"
        \\ Decide if we can actually lower precision here. We might not
        default(realbitprecision, gd_precision);    \\ reduce precision
        \\\ ... do some computation
        default(realbitprecision, mainbitprecision); \\ restore original precision
        \\ s_radius = bitprecision(s_radius, mainbitprecision); \\add zeroes so that pari doesn't do future calcs at low precision
    */

    for(i=1, length(matches),
        print(precision(matches[i],10));
        print(matchIdeal, "    ", giant_divisor[1]);
        baby_element = compact_rep_buchmann(G, matches[i], giant_divisor[1], eps);
        psi_baby = log_from_cpct(G, baby_element)[1..r];
        print(precision(matches[i],20), "   ", precision(psi_baby,20));
        gd_log = trackerLogarithm(G, ~compactTracking, r);
        new_vec = gd_log - psi_baby;

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


\\ Subalgorithm of bsgs. computes giant steps incrementally using ideal
\\ multiplication to obtain adjacent elements.
\\ Compact Represenation Version
incremental_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, ~babystock, avec,\
                        eps, outFileInfo=[])=
{
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
        giant_coeffs = zero_vec,        current_giant_vec = zero_vec,
        giant_divisor,                  expected_position,
        matches,                        new_vec,
        place_marker = r,               directions
    );
    GP_ASSERT_EQ(r, length(giant_sublattice));

    expected_position = vector(G.r1+G.r2-1, i, 0)~;
    directions = vector(length(giant_sublattice), i , 1);
    giant_coeffs[r] = 1;    \\# initialize variable to track giant elements at 1

    my(
        \\# vectors of length r which are used to compute an adjacent element in
        \\# a particular direction in the logarithm lattice
        \\# direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
        direction_elements,                 \\# set of corresponding to a giant step in each lattice dimension
        inverse_direction_elements,         \\# inverse elements of the above
        compactTracking = List(),           \\# list of cpct reps corresponding to current giant step
        trackingLog = vector(r, i, 0),      \\# log corresponding to current giant step. Precision drops fast
        s_radius                            \\# max acceptable distance from expected position in log lattice
    );

    mainbitprecision = default(realbitprecision);
    default(realbitprecision, localbitprecision);
    s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    default(realbitprecision, mainbitprecision);
    s_radius = bitprecision(s_radius, mainbitprecision);

    \\# for some reason I have to flip these. giantstep() may actually return an inverse
    \\# initialize direction elements and inverses
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, giant_sublattice, field_deg, eps);

    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );

    giant_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Additional timing variables and file write in jump");
    my(giant_t1, giant_tn, giant_tmid, ctr,
        tDivisorCompute, tDivisorReduce,
        tCompare, tNext);
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    tDivisorCompute = 0;
    tDivisorReduce = 0;
    tCompare = 0;
    tNext = 0;


    SCREEN("Consider checking for inverses when comparing giant step divisors");
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
            trackingLog += log_from_cpct(G, direction_elements[place_marker][3])[1..r];
            default(realbitprecision, mainbitprecision);
        ,\\else
            giant_divisor = [idealmul(G, giant_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],inverse_direction_elements[place_marker][2] )~];

            compactTracking[place_marker][2] -= 1;
            default(realbitprecision, localbitprecision);
            expected_position -= giant_sublattice[,place_marker];
            trackingLog -= log_from_cpct(G, direction_elements[place_marker][3])[1..r];
            default(realbitprecision, mainbitprecision);
        );
        bitprecision(expected_position, mainbitprecision);
        tDR = getabstime();
        tDivisorCompute += (tDR-tDC);
        get_next_giant_divisor_cpct(~G, ~giant_divisor, ~compactTracking);
        t_half = getabstime();
        tNext += (t_half-tDR);

        default(realbitprecision, localbitprecision);
        trackingLog += log(abs(nfeltembed(G, compactTracking[length(compactTracking)][1])))[1..r];
        default(realbitprecision, mainbitprecision);
        trackingLog = bitprecision(trackingLog, mainbitprecision);

        [giant_divisor, tempLog] = adjust_giant_step_cpct(~G, ~giant_divisor, ~compactTracking, bitprecision(trackingLog,REQ_BSGS), ~expected_position, bitprecision(s_radius, REQ_BSGS), eps);

        tAdjust = getabstime();
        tDivisorReduce += (tAdjust - t_half);

        EXACT_CHECK = 0;

        \\ use babystock hashmap to check for collisions and update as needed
        if(mapisdefined(babystock, giant_divisor[1]),
            \\if(DEBUG_BSGS>0 ,print("baby-giant match. checking for new vector:"););
            matches = mapget(babystock, giant_divisor[1]);

            if (!EXACT_CHECK,
                lattice_lambda = collision_check1(G, ~lattice_lambda, r, eps,
                                 ~matches, ~giant_divisor, ~compactTracking);
            ,
                lattice_lambda = collision_check2(G, ~lattice_lambda, r, eps,
                                giant_divisor[1], ~matches, ~giant_divisor, ~compactTracking);
                print("collision check success");
            );
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
        if(ctr %1000 ==0,
            print(ctr, " giant steps completed");
            print("Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);
            \\print("Adjust log: ", adjusttime1);
            giant_tn = getabstime();
            write("data/jump-timing.txt", ctr, " jumps took: ", (giant_tn- giant_tmid), "milliseconds");
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
        giant_coeffs = zero_vec,        current_giant_vec = zero_vec,
        giant_divisor,                  expected_position,
        matches,                        new_vec,
        place_marker = r,               directions
    );
    GP_ASSERT_EQ(r, length(giant_sublattice));

    expected_position = vector(G.r1+G.r2-1, i, 0)~;
    directions = vector(length(giant_sublattice), i , 1);
    giant_coeffs[r] = 1;    \\# initialize variable to track giant elements at 1

    my(
        \\# vectors of length r which are used to compute an adjacent element in
        \\# a particular direction in the logarithm lattice
        \\# direction_elements[i] and inverse_direction_elements[i] are multiplicative inverses
        direction_elements,                 \\# set of corresponding to a giant step in each lattice dimension
        inverse_direction_elements,         \\# inverse elements of the above
        compactTracking = List(),           \\# list of cpct reps corresponding to current giant step
        trackingLog = vector(r, i, 0),      \\# log corresponding to current giant step. Precision drops fast
        s_radius                            \\# max acceptable distance from expected position in log lattice
    );

    mainbitprecision = default(realbitprecision);
    default(realbitprecision, localbitprecision);
    s_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    default(realbitprecision, mainbitprecision);
    s_radius = bitprecision(s_radius, mainbitprecision);

    \\# for some reason I have to flip these. giantstep() may actually return an inverse
    \\# initialize direction elements and inverses
    [direction_elements, inverse_direction_elements] =
        get_giant_step_increment_vectors_compact(G, giant_sublattice, field_deg, eps);

    for( i=1, length(direction_elements),
        listput(~compactTracking, [direction_elements[i][3], 0]);
    );

    giant_divisor = [matid(field_deg), vector(r+1, i,1 ), [[1],[1]] ];

    SCREEN(0, "Additional timing variables and file write in jump");
    my(giant_t1, giant_tn, giant_tmid, ctr,
        tDivisorCompute, tDivisorReduce,
        tCompare, tNext);
    giant_t1 = getabstime();
    giant_tmid = giant_t1;

    tDivisorCompute = 0;
    tDivisorReduce = 0;
    tCompare = 0;
    tNext = 0;


    SCREEN("Consider checking for inverses when comparing giant step divisors");
    \\ #increments are done modulo the product of avec, returns to zero when
    \\ #all elements have been visited
    while(giant_coeffs != zero_vec,

        tDC = getabstime();
        if(directions[place_marker] == 1,
            giant_divisor = [idealmul(G, giant_divisor[1], direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],direction_elements[place_marker][2] )~];
            if(storage == "COMPACT",
                compactTracking[place_marker][2] += 1;
            );
            if(storage == "LOG",
                default(realbitprecision, localbitprecision);
                expected_position += giant_sublattice[,place_marker];
                trackingLog += log_from_cpct(G, direction_elements[place_marker][3])[1..r];
                default(realbitprecision, mainbitprecision);
            );

        ,\\else
            giant_divisor = [idealmul(G, giant_divisor[1], inverse_direction_elements[place_marker][1]),
                pointwise_vector_mul(giant_divisor[2],inverse_direction_elements[place_marker][2] )~];

            if(storage == "COMPACT",
                compactTracking[place_marker][2] -= 1;
            );
            if(storage == "LOG",
                default(realbitprecision, localbitprecision);
                expected_position -= giant_sublattice[,place_marker];
                trackingLog -= log_from_cpct(G, direction_elements[place_marker][3])[1..r];
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
            trackingLog += log(abs(nfeltembed(G, compactTracking[length(compactTracking)][1])))[1..r];
            default(realbitprecision, mainbitprecision);
            trackingLog = bitprecision(trackingLog, mainbitprecision);
        );

        [giant_divisor, tempLog] = adjust_giant_step_cpct(~G, ~giant_divisor, ~compactTracking, bitprecision(trackingLog,REQ_BSGS), ~expected_position, bitprecision(s_radius, REQ_BSGS), eps);

        tAdjust = getabstime();
        tDivisorReduce += (tAdjust - t_half);

        EXACT_CHECK = 0;

        \\ use babystock hashmap to check for collisions and update as needed
        if(mapisdefined(babystock, giant_divisor[1]),
            \\if(DEBUG_BSGS>0 ,print("baby-giant match. checking for new vector:"););
            matches = mapget(babystock, giant_divisor[1]);
            if (storage == "COMPACT",
                matches_cpct = mapget(cpct_bstock, giant_divisor[1]);
            );

            if (storage == "LOG",
                if (!EXACT_CHECK,
                    lattice_lambda = collision_check1(G, ~lattice_lambda, r, eps,
                                     ~matches, ~giant_divisor, ~compactTracking);
                ,
                    lattice_lambda = collision_check2(G, ~lattice_lambda, r, eps,
                                    giant_divisor[1], ~matches, ~giant_divisor, ~compactTracking);
                    print("collision check success");
                );
            );
            if (storage == "COMPACT",
                print("Compact collision checker");
                lattice_lambda = collision_check_cpct(G, ~lattice_lambda, r, eps,
                                 ~matches_cpct, ~giant_divisor, ~compactTracking)
            );
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
        if(ctr %1000 ==0,
            print(ctr, " giant steps completed");
            print("Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);
            \\print("Adjust log: ", adjusttime1);
            giant_tn = getabstime();
            write("data/jump-timing.txt", ctr, " jumps took: ", (giant_tn- giant_tmid), "milliseconds");
            giant_tmid = giant_tn;
        );
    );
    print("Compute: ", tDivisorCompute, " Reduce1 ", tNext," Adjust ", tDivisorReduce, " Compare: ", tCompare);

    return(lattice_lambda);
}


trackerLogarithm(~G, ~tracker, r)={
    logResult = vector(r,j,0);
    for(i=1, length(tracker),
        if(i <= r,
            logResult += tracker[i][2]*log_from_cpct(G, tracker[i][1])[1..r];
        ,\\else
            logResult+= tracker[i][2]*log(abs(nfeltembed(G, tracker[i][1])))[1..r];
            \\logResult+= tracker[i][2]*log(abs(G[5][1]*(tracker[i][1])))[1..r]~;
        );
    );
    return(logResult);
};

\\ gets a nearby reduced divisor to log_coordinates_Rn
\\ INPUTS:
\\ - G a number field, idealmat a coefficient matrix of an ideal
\\ - log_coordinates_Rn is the log vector (length = equal to unit rank)
\\ - inverse is an optional argument that will compute the inverse instead if set to 1
get_nearby_rdivisor(G, idealmat, log_coordinates_Rn, inverse=1)={
    exponentiated_coordinates = create_target(G, log_coordinates_Rn, inverse);
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
        \\print("\ntarget position = ", precision(expected_position,10), "\nOriginal Distance from Target ", precision(norml2(divisor_distance),10) );
        for(i=1, 2,
            adjustment_divisor = get_nearby_rdivisor(G, matid(poldegree(G.pol)), divisor_distance, i%2);
            if (norml2(adjustment_divisor[2]) < eps,
                return(giant_divisor);
            ,
                new_divisor = multiply_and_reduce_divisors(G, giant_divisor[1], giant_divisor[2], giant_divisor[3], adjustment_divisor[1], adjustment_divisor[2], adjustment_divisor[3]);
                new_distance = norml2(expected_position - (new_divisor[3][1..r]~));
                \\print("Adjusted distance from target = ", precision(new_distance, 10) );
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
    my(new_divisor);
    \\print("b ", giant_divisor[2]);
    new_divisor = reddiv_compact(~giant_divisor[1], ~giant_divisor[2], ~G,~G[5][1]);
    \\print("a ", new_divisor[2]);
    listput(~tracker, [new_divisor[4], 1]);
    giant_divisor[1] = new_divisor[1];
    giant_divisor[2] = new_divisor[2];

}

\\ INPUT:
\\ - G is a number field
\\ - cpct_reps is a list of cpct representations corresponding to a sublattice of the unit lattice
\\ - B is an integer that indicates a fraction of the fundamental region to search
\\   Divides the last coordinate vector by B.
\\ - babystock_scale_factor indicates the volume of the babystock region
\\ - eps is the error, REQ_BSGS is the calculated precision needed for accurate results
\\ - options is a vector of additional arguments.
\\ - options[1] may specify a different algorithm to scan with (not available yet)
\\ - or it specifies a timeout value, which will be checked periodically and
\\ -    the program terminated if exceeded
bsgs(G, cpct_reps, B, babystock_scale_factor, scanballRadius,eps, REQ_BSGS,FILE1="data/tmp-bsgs-log", options = [])={

    print("\nBSGS Algorithm Start"); bsgs_start = getabstime();
    my(S_radius, r = G.r1 + G.r2 -1,
        field_deg = poldegree(G.pol),
        zero_vec = vector(r, i , 0),                                            \\ length r zero vector, used to check when giantsteps are done
        int_ring_real = embed_real(G, G[5][1]),
        avec,
        giant_coeffs,                                                           \\ holds numerator of coefficients for the giant steps
        giant_sublattice,                                                             \\ matrix of columns of lambda, each scaled by a_i
        current_giant_vec,
        babystock_region,
        babystock = Map(),
        giant_divisor,
        matches,
        new_vec,lattice_lambda,
        timeout = 0,
        alg = "SCAN",
        storage = "COMPACT"
    );
    if (storage == "COMPACT", cpct_bstock = Map());

    \\ check if an auxilliary string for output has been provided
    if (length(options)!=0,
        if(type(options[1]) == "t_STR",
            alg = options[1];,
            timeout = options[1];
        );
    );

    if(alg != "NEIGHBOURS" && alg != "SCAN", print("Unknown algorithm specified, returning -1", return(-1) ));
    if(DEBUG_BSGS>0 , print("Bound B = ", precision(B,10)); );

    lattice_lambda = log_lattice_from_compact_set(G, cpct_reps);
    S_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));
    [avec, giant_sublattice] = get_giant_step_params(G,lattice_lambda, r, B, babystock_scale_factor, REQ_BSGS);


    tb = getabstime();

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #THIS OPTION USES NEIGHBORS TO SEARCH THE BABYSTOCK
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(alg == "NEIGHBOURS",
        \\# note this option has been disabled for some time so it no longer works
        print("NEIGHBOURS option has been disabled. Rerun using SCAN"); quit;
        neighbour_search(G,lattice_lambda, babystock_region, ~babystock, giant_sublattice, eps);
    );  \\# end neighbors for babystock

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #THIS SECTION SEARCHES USING THE SCAN ALGORITHM
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    minimal_vec_length = norml2(giant_sublattice[,1]);
    for(ii = 1, length(giant_sublattice),
        if (norml2(giant_sublattice[,ii]) < minimal_vec_length, minimal_vec_length = norml2(giant_sublattice[,ii]));
    );
    minimal_vec_length = sqrt(minimal_vec_length);
    if(minimal_vec_length < scanballRadius,
        scanballRadius = minimal_vec_length;
        write(FILE1,"Adjusted scan radius: ", precision(scanballRadius,10), "  aVec =", avec);
        ,
        write(FILE1, "Scan radius: ", precision(scanballRadius,10), "  aVec =", avec);
    );
    if(alg == "SCAN",
        my(foundflag = 1, num_elts_found = []);

        if (storage == "COMPACT",
            [lattice_lambda, num_elts_found] =
                incremental_baby_steps_storage(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, storage, ~cpct_bstock, [FILE1, timeout]);
                breakpoint();
        ,
                [lattice_lambda, num_elts_found] =
                    incremental_baby_steps(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, [FILE1, timeout]);
        );

        foundflag = length(num_elts_found);
        print("babystock size: ", length(babystock), " cpct bstock size: ", length(cpct_bstock));
        while(foundflag !=0,
            \\ if new elements found, regenerate lattice_lambda to remove precision loss
            lattice_lambda = log_lattice_from_compact_set(G, cpct_from_loglattice(G, lattice_lambda, eps));
            [avec, giant_sublattice] = get_giant_step_params(G,lattice_lambda, r, B, babystock_scale_factor, REQ_BSGS);

            if (storage == "COMPACT",
                [lattice_lambda, num_elts_found] =
                    incremental_baby_steps_storage(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, storage, cpct_bstock, [FILE1, timeout]);
            ,
                    [lattice_lambda, num_elts_found] =
                        incremental_baby_steps(matid(field_deg),~lattice_lambda, ~giant_sublattice, ~babystock, G, scanballRadius, eps, [FILE1, timeout]);
            );
            GP_ASSERT_EQ(matsize(lattice_lambda)[1], matsize(lattice_lambda)[2]);
            foundflag = length(num_elts_found);
            if(DEBUG_BSGS>0 ,print("FoundFlag ",foundflag););
        );
    );
    babytime = getabstime()-tb;
    write(FILE1, "Babystock time: ", precision(babytime, 10), "   ",  precision(babytime/60000.0, 15), ". Babystock elements: ", length(babystock));

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ #Begin giant step computations
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    print("\nGiant Loop\n");
    my(aProduct = 1);
    for(i=1, length(avec), avec[i] += 1; aProduct*=avec[i];);                                       \\ adjust avec to include boundaries.
    avec = round(avec);
    print("avec just before giant steps: ", precision(avec,10)); \\breakpoint();
    \\# the search region area is the product of the norms of each sublattice vector
    \\#normproduct = 1; for(i=1, length(giant_sublattice), normproduct*=sqrt(norml2(giant_sublattice[,i])) );
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    tg = getabstime();

    lattice_lambda =
        incremental_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, \
                                ~babystock, avec, eps, [FILE1, timeout]);
    \\lattice_lambda = jump_giant_steps(~G, ~lattice_lambda, ~giant_sublattice, ~babystock, avec, eps);

    bsgs_end = getabstime();
    gianttime = bsgs_end-tg;

    write(FILE1,  "Giantstep time:  ", gianttime,  "   ",precision(gianttime/60000.0,15), " . Giant steps computed: ", aProduct);

    return(cpct_from_loglattice(G,lattice_lambda,eps));
}
