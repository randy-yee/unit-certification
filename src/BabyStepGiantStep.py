read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/CompactRepresentation.py");
read("src/Neighbours.py")
read("src/bounds.gp")

/*
\\ BSGS related routines
get_subdivisions
increment_coordinates
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
minimum_check(ideal_matrix, elt)={


}



\\ Determines the size of the babystock area
\\ INPUT:
\\    dimensions, an integer equal to r, the rank of the lattice (in FJ, this is r1+r2-1)
\\    detLambda, the determininant of the input lattice Lambda
\\    B the prime search bound.
\\    delta is a free variable used in the exponent of the runtime expression.
\\ OUTPUT:
\\ - a vector of integers of length = length(Lambda)
\\ - The product of the integers is close to (detLambda/B)^(1-delta)
get_subdivisions(G,lattice_lambda, dimensions, detLambda, B, delta,REQ_BSGS)={
    my(smallsquare, sides, avec, a_product, ind=1, fullregion_to_babystock_ratio= 1);
    \\smallsquare = (abs(detLambda)/B)^(1-delta);                                 \\ target volume of our box

    smallsquare = sqrt(  (abs(detLambda)/B)*giant_n(poldegree(G.pol), log(abs(G.disc)), REQ_BSGS,real(log(detLambda)))/baby_n( poldegree(G.pol),log(abs(G.disc)),REQ_BSGS,real(log(detLambda)) )  );

    if(DEBUG_BSGS, print("Expected babystock region size: ", precision(smallsquare,10) ););

    babystock_scale_factor =2;                               \\ increase this to make the babystock smaller

    fullregion_to_babystock_ratio = (babystock_scale_factor)*abs(detLambda)/(smallsquare*B);


    if(DEBUG_BSGS,
        print("G_n, B_n  ",ceil(giant_n(poldegree(G.pol), log(abs(G.disc)), REQ_BSGS,real(log(detLambda)))), "   ", precision(baby_n( poldegree(G.pol),log(abs(G.disc)),REQ_BSGS,real(log(detLambda)) ) ,10));
        print("a_I product target: ", ceil(fullregion_to_babystock_ratio), "   " );
    );

    norm_vr = ceil( sqrt(norml2(lattice_lambda[,dimensions] )) );
    if(DEBUG_BSGS,print("norm of rth vector: ", norm_vr););

    if(dimensions == 1,
        avec = vector(dimensions, i , ceil(fullregion_to_babystock_ratio));
    , \\else
        sides = max(1,floor(sqrtn(fullregion_to_babystock_ratio, dimensions-1)));             \\ approximate the sides with integers

        if(DEBUG_BSGS,print("(r-1)th root of target area  ", sides););
        avec = [];
        my(diffvector = [],
            vecnorm =1,
            ai_product=1,
            areadiff = 0,
            partialproduct = 1);

        for(i=1, dimensions-1,
            vecnorm =  floor( sqrt(norml2(lattice_lambda[,i] )) );
            \\print(vecnorm , "  ", sides);
            avec = concat(avec, min(vecnorm, sides));
            \\print("vecnorm, sides ", vecnorm, "  ", sides );
            diffvector = concat(diffvector, vecnorm - sides );
            ai_product *= avec[i];
        );
        for(i=1, length(diffvector),

            if(diffvector[i] >0,
                \\print("considering increasing a_", i);
                partialproduct = (ai_product/avec[i]);
                areadiff = round((fullregion_to_babystock_ratio - ai_product)/(ai_product/avec[i]));
                if(areadiff >= 1, avec[i] += min(areadiff, diffvector[i]); ai_product += partialproduct*min(areadiff, diffvector[i]) );
            );
        );

        if(DEBUG_BSGS, print("product without  a_r ",ai_product); );
        avec = concat(avec ,  min(max(1, round(norm_vr/B)  ), ceil(fullregion_to_babystock_ratio/ai_product ))  );
    );



    \\ This is a special fix for a particular rank 2 example from Ha, try to generalize later!
    \\    if(G.r1 +G.r2-1 == 2,
    \\        avec = vector(dimensions, i , sides);
    \\        avec[1]=round(sqrt(fullregion_to_babystock_ratio*sqrt(norml2(lattice_lambda[,1]))/(sqrt(norml2(lattice_lambda[,2]))/B ) ));
    \\        avec[2]=ceil(1/(fullregion_to_babystock_ratio/avec[1]));
    \\        print("new avec ", avec[1], "  ", avec[2]); breakpoint();
    \\        ,
    \\        avec[1] = ceil(fullregion_to_babystock_ratio);
    \\        avec[2] =1;
    \\        while(a_product < fullregion_to_babystock_ratio,
    \\            avec[ind]+=1;
    \\            a_product/= sides;
    \\            a_product *= (sides+1);
    \\            ind += 1;
    \\        );
    \\    );

    finalproduct =1; for(i=1, length(avec), finalproduct*=avec[i]); print(avec, " actual a_i product: ", finalproduct, " target: ", round(fullregion_to_babystock_ratio));
    \\    print("vector of a_i ", avec);
        return(avec);
}


\\ INPUT:
\\ - A vector of integers A which defines the subdivision for the babystock
\\ - A vector of integers V that defines a particular 'giant step'
\\   for all i, V[i] < A[i]
\\ OUTPUT:
\\ - The vector V 'incremented' by 1, similar to a number in which each digit
\\  is a different base
increment_coordinates(a_vec, current_vec)={
    my(
        place = length(current_vec),
        carryflag = 0;
    );
    current_vec[place]+=1;
    if(current_vec[place] >= a_vec[place], carryflag =1);

    while(carryflag == 1,
        current_vec[place] = 0; carryflag = 0;
        if(place == 1, return(current_vec));
        place -= 1;
        current_vec[place]+=1;
        if(current_vec[place] >= a_vec[place], carryflag =1);
    );

    return(current_vec);
}

\\ obtains an axis aligned box which contains the box generated by the columns of basismat
\\ represented by two points in R^r, where r is the unit rank
get_axis_aligned_box(basismat)={
    \\print("Getting aligned box:");
    my(r = length(basismat), twovec, maxcoords, mincoords, vertex_coeffs, vertex );
    twovec = vector(length(basismat), i, 2);
    zerovec = vector(r, i, 0);
    vertex_coeffs = zerovec;
    maxcoords = vertex_coeffs;
    mincoords = vertex_coeffs;
    vertex_coeffs = increment_coordinates(twovec, vertex_coeffs);
    \\print(vertex_coeffs);
    while(vertex_coeffs != zerovec,

        vertex = basismat*(vertex_coeffs~);
        for(i=1, r,
            if(vertex[i] < mincoords[i], mincoords[i] = vertex[i]);
            if(vertex[i] > maxcoords[i], maxcoords[i] = vertex[i]);
        );
        vertex_coeffs = increment_coordinates(twovec, vertex_coeffs);
        \\print(vertex_coeffs);
    );
    \\print("length of box coords ",length(mincoords));
    return([mincoords, maxcoords]);
}

\\ subroutine of the bsgs method.
\\ Determines a basis for the giant steps (coefficients limited by avec)
\\ INPUT:
\\ - lattice_lambda is the current log lattice
\\ - r is the unit rank
\\ - B defines the fraction of the fundamental region F to search, 1/B
\\ - delta determines how big the babystep region is in relation to the search region
\\ OUTPUT:
\\ - avec is a r-vector of coefficient upperbounds for the giant steps
\\ - giant_legs is a basis that determines the giant steps
\\ - babystock_region is an axis aligned box that contains the babystock region. Given by two point vectors.
get_giant_step_params(G,lattice_lambda, r,B, delta,REQ_BSGS)={
    my(det_lambda, giant_legs, babystock_region);
    det_lambda = matdet(lattice_lambda);

    originalarea = 1;
    if(DEBUG_BSGS,
        for(i=1, length(lattice_lambda),
            print("length of v_",i,"   ",  precision(sqrt(norml2(lattice_lambda[,i])),10));
            originalarea*=sqrt(norml2(lattice_lambda[,i]))
        );
        print("original area by multiplying norms: ", precision(originalarea,10), "  det_lambda: ", precision(det_lambda,10) );
    );
    avec = get_subdivisions(G,lattice_lambda,r, det_lambda, B, delta,REQ_BSGS);

    giant_legs = matrix(r, r-1, s,t, -lattice_lambda[s,t]/avec[t]);             \\ giant_legs is a matrix with columns (-v_i /a_i) for i =1 ..r-1


    \\print("vectors only scaled by a", precision(norml2(giant_legs[,r]),10)); breakpoint();
    giant_legs = matconcat([giant_legs, -lattice_lambda[,r]/(avec[r]*B)] );     \\ column r is (-v_r /(a_r*B))
    if(DEBUG_BSGS,
    print("- avec ", avec, "  ", B);
    area = 1;
    for(i=1, length(giant_legs),
        print("Norms of scaled vectors (a_i, B) ", precision(sqrt(norml2(giant_legs[,i])),10));
        area*=sqrt(norml2(giant_legs[,i])));
    print("actual babystock area region (multiplying vector norms): ", precision(area,10),"\n \n");
    );
    babystock_region = get_axis_aligned_box(giant_legs);

    return([avec,giant_legs, babystock_region]);
}

/******************************************************************************/
/*28. Find a new vector v \in \Lambda_K\Lambda' where (O_K, v) \in \ep_B: norml2, fun. in 13 (is_vec_in_lattice),  fun. in 11 (updatelamda)*/
\\ Check the set ball_minima for a new vector in Lambda not already in Lambda'
\\ INPUT:
\\ - ball_minima a list of minima, each of the form [ideal, log vector] = [1/mu*O_K, vector = log(mu)]
\\ - L is the lattice we are checking (Lambda')
\\ - G the number field
\\ - eps the error

PRECISIONLOSS = 10^5
/******************************************************************************/
check_babystock_for_units(ball_minima,L,G,eps)={

    my(new_counter:small, candidate, ideal_identity, eps_sqr = eps^2, n:small = poldegree(G.pol));
    ideal_identity = matid(n);
    new_counter=0;                                                                      \\ track the number of new minima found.

    for(i=1,length(ball_minima),
        if(ball_minima[i][1] == ideal_identity,                                         \\ if the ideal (1/mu)*O_k = O_k then check further
            candidate=ball_minima[i][2];                                                \\ candidate = log(mu)
            if(norml2(candidate)>eps_sqr&&is_vec_in_lattice(candidate~,L,eps)==0,           \\ if nonzero and v is not already in L, then
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
scanball(G, y,u,psimu,web,eps)={
    my(n, ball_minima=[], x, scan_bound, vecholder, gram_mat, scan_elements, vec_numerical, LLL_reduced_yu, new_y);
    n = poldegree(G.pol);
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
        vecholder=scan_elements[,ii];                                                     \\ get the ith element w_i
        eltnorm = abs(nfeltnorm(G,vecholder ));
        if(eltnorm>=1 && eltnorm<=norm_deltaK,                                            \\ if scan_element[ii] is a minimum, its norm satisfies 1 < N(s_elt) < N(y)*delta_K
            new_y=idealdiv(G,y,vecholder);                                                    \\ get the ideal (1/w_i)*y
            new_yLLL = embed_real(G, G[5][1]*new_y)*qflll(embed_real(G, G[5][1]*new_y));
            if(norml2(new_yLLL) > 1-eps,
                if(checkred_old(new_y,G,eps)==1,                                                  \\ check if its reduced, and if yes, add add to ball_minima
                    vec_numerical = (G[5][1]*vecholder)~;                                            \\ get the numericals
                    \\print("SCANBALL: Psi of potential minima: ",precision(vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps), 10 ));
                    ball_minima = setunion(ball_minima, Set([[new_y,vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1]))+psimu,eps)]] ));
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

is_repeat_babystock(existing_entry, candidate_embedding,eps)={
    my(repeatflag = 0);
    for(t=1, length(existing_entry),
        if(samevecs(existing_entry[t], candidate_embedding,eps), repeatflag = 1);
    );
    return(repeatflag);
}

compute_max_point_distance(G, unit_rank, webpoint)={
    my(web_distance);
    web_distance = vector(unit_rank, i,if(i <= G.r1, webpoint[i], 2*webpoint[i]););     \\ length r vector, scale complex entries
    web_distance = sqrt(norml2(concat(web_distance,sumvec(web_distance))));
    return(web_distance)
}

get_web_step_multipliers(G, unit_rank,web_increments)={
    my(tempvec, zerovec = vector(unit_rank, i, 0), web_step = [[],[]]);
    for(i = 1, unit_rank,
        tempvec = zerovec; tempvec[i] = web_increments[i];                      \\ vector corresponding a move in the ith direction
        web_step[1] = concat(web_step[1], [create_target(G, tempvec)]);                 \\ get corresponding exponential
        web_step[2] = concat(web_step[2], [invert_coordinates(web_step[1][i])]);        \\ get corresponding inverse
        \\print("web_step[1]=",precision(web_step[1][i],10), "   \nweb_step[2]=", precision(web_step[2][i],10), "\n" );
    );
    return(web_step);
}

get_initial_webpoint(G, unit_rank, babystock_corner, web_increments)={
    my(webpoint, web_distance, exp_webpoint);
    webpoint = 1/2*(web_increments);
    web_distance = compute_max_point_distance(G, unit_rank, webpoint);              \\ gets the value epsilon from Alg10.7, the max distance between web points

    webpoint = webpoint + babystock_corner;                                             \\ define the first webpoint relative to the 'small corner'
    exp_webpoint = create_target(G, webpoint);
    return([webpoint, exp_webpoint, web_distance]);
};

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

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Begin scanning the points in the web of regularly distributed points
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    [ideal_J, nu, logdist, beta] = reddiv_compact(y, exp_webpoint, G, G[5][1]);

    ball_minima = scanball(G, ideal_J, nu, logdist[1..length(L)], web_distance, eps);

    region_minima = Set(concat(region_minima,ball_minima));
    \\ check for new elements in Lambda

    directions = vector(length(L), i , 1);                                              \\ initialize as all positive directions

    while(web_coords != zerovec || next_coords==[],                             \\ loops over all of the web points in the babystock region
        \\print(web_coords);
        next_coords = increment_coordinates(box_subdivisions, web_coords);      \\ compute the next point's coefficients
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

        ball_minima = scanball(G, ideal_J, nu, logdist[1..length(L)], web_distance, eps);

        [L, newctr] = check_babystock_for_units(ball_minima,L,G,eps);
        if(newctr != 0,
            print("Found a unit in babystock. new reg =", precision(abs(matdet(L)),10) );
            babystock_box[1][1] += (web_coords[1]*web_increments[1]);
            return([L, baby_hashmap, babystock_box]);
        );

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
    print("Babysteps stored");
    return([L,baby_hashmap, []]);
}


babystock_scan_jump(y,L,giant_legs,G,eps)={
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
        region_minima=Set(),
        directions,
        ball_minima,
        newctr = 0,
        existing_entry,
        field_deg = poldegree(G.pol),
        identity_n = matid(field_deg)
    );

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Establish a number of parameters including the web of regularly distributed points,
    \\ and some values that allow us to move easily between each of the points sucessively.
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    \\ subdivide the babystock region based on the lengths of the defining vectors
    denoms = vector( length(L), i, ceil(norml2(giant_legs[,i])));
    babystock_t = giant_legs;
    for(i=1, length(babystock_t), babystock_t[,i]/denoms[i] );

    xvector = vector(matsize(giant_legs)[2], i,0)~;
    web_coords = increment_coordinates(denoms, web_coords);

    while(web_coords != zerovec,                                               \\ loops over all of the web points in the babystock region

        xvector = babystock_t*web_coords~;                                     \\ this holds the web point to be scanned
        web_coords = increment_coordinates(denoms, web_coords);                \\ compute the next point's coefficients

        [ideal_J, nu, logdist] = giantstep(y, xvector, G, field_deg, eps);

        ball_minima = scanball(G, ideal_J, nu, logdist[1..length(L)], web_distance, eps);
        [L, newctr] = check_babystock_for_units(ball_minima,L,G,eps);
        if(newctr != 0,
            print("Found a unit in babystock. new reg =", precision(abs(matdet(L)),10) );
            return([L, baby_hashmap, [1] ]);
        );

        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\ Take the elements found on the ball near the current point in the web, and enter them into the hash map
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\region_minima = setunion(region_minima, ball_minima); \\ Set(concat(region_minima,ball_minima));
        for(ctr=1, length(ball_minima),
            if(mapisdefined(baby_hashmap, ball_minima[ctr][1], &existing_entry),
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
    print("Babysteps stored");
    return([L,baby_hashmap, []]);
}


\\ subroutine of bsgs - performs giant steps using jump algorithm one at a time and checks with babystock
\\ INPUT:
\\ - G a number field and lat_lambda is unit lattice
\\ - gs_sublattice is a sublattice of the fundamental domain of lat_lambda
\\ - bstock is a Map object of the babystock. Keys are ideals, values are Psi values of generating elements
\\ - avec indicates the number of subdivisions of lat_lambda were taken in each coordinate to get gs_sublattice
\\ - eps is an error value
\\ RETURN: a potentially expanded unit_lattice.
jump_giant_steps(G, lat_lambda, gs_sublattice, bstock, avec, eps)={
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
                if(norml2(new_vec) > (eps*10)^2 && is_vec_in_lattice(new_vec~,lat_lambda,eps)==0,
                    lat_lambda = my_mlll( matconcat([lat_lambda, new_vec~]),eps);
                    if(DEBUG_BSGS>0, print("New element found. New regulator = ", precision(matdet(lat_lambda),10)););
                );
            );
        );
        giant_coeffs = increment_coordinates(avec, giant_coeffs);

        ctr++;
        if(ctr %500 ==0,
            giant_tn = getabstime();
            write("data/jump-timing.txt", ctr "jumps took: ", (giant_tn- giant_tmid), "milliseconds");
            giant_tmid = giant_tn;
        )
    );
    write("data/jump-timing.txt", "Total jump time: ", (giant_tn- giant_t1), "milliseconds\n");
    return(lat_lambda);
};


\\ gets a nearby reduced divisor to log_coordinates_Rn
\\ INPUTS:
\\ - G a number field, idealmat a coefficient matrix of an ideal
\\ - log_coordinates_Rn is the log vector (length = equal to unit rank)
\\ - inverse is an optional argument that will compute the inverse instead if set to 1
get_nearby_rdivisor(G, idealmat, log_coordinates_Rn, inverse=0)={
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
        \\print("A closer vector was not found");
        return(giant_divisor);
    );
};


\\ subalgorithm for bsgs, adjust by the radius S
expand_babystock_region(babystock_region,S_radius )={
    my(S_radius_vector);
    if(DEBUG_BSGS>1 ,print("Babystock box: ", precision(babystock_region,10)););
    S_radius_vector = vector(length(babystock_region[1]), i , S_radius);
    babystock_region[1] = babystock_region[1] - S_radius_vector;
    babystock_region[2] = babystock_region[2] + S_radius_vector;

    return(babystock_region);
}

\\ subalgorithm in the bsgs method, get the initial neighbour to start
\\ the babystock search when using the neighbours method
get_initial_neighbour(G, giant_sublattice, field_deg,eps)={
    my(middle_babystock);
    middle_babystock = vector(length(giant_sublattice), i, 1/2);
    middle_babystock = giant_sublattice*(middle_babystock~);
    return(giantstep(matid(field_deg), middle_babystock, G, field_deg, eps);)
}

get_giant_step_increment_vectors(G, giant_sublattice, field_deg, eps)={
    my( giant_step_ideals = [], giant_step_ideals_inverse = [],
        new_vec,gstep_divisor);

    for(j = 1, length(giant_sublattice),
        new_vec = giant_sublattice[,j];
        gstep_divisor = giantstep(matid(field_deg), new_vec, G, field_deg, eps);
        \\debug_compare(gstep_divisor, giantstep_high_precision(matid(field_deg), new_vec, G, field_deg, eps));

        \\current_giant_vec = create_target(G, new_vec);                                    \\ when do I use the argument 1 here? Used in giantstep function: create_target(G, new_vec,1);
        \\gstep_divisor = reddiv_compact(matid(field_deg), current_giant_vec, G, G[5][1]);

        giant_step_ideals = concat(giant_step_ideals, [[gstep_divisor[1], gstep_divisor[2], gstep_divisor[3]]]);
        giant_step_ideals_inverse = concat(giant_step_ideals_inverse, [[idealinv(G,gstep_divisor[1]), invert_coordinates(gstep_divisor[2]),-gstep_divisor[3]  ]] );
    );

    return([giant_step_ideals,giant_step_ideals_inverse]);
}

get_next_giant_divisor(G, giant_divisor)={
    my(new_divisor);
    new_divisor = reddiv_compact(giant_divisor[1], giant_divisor[2], G,G[5][1]);
    if(type(new_divisor[3])== "t_COL", new_divisor[3] = new_divisor[3]~);
    return([new_divisor[1], new_divisor[2], giant_divisor[3]+new_divisor[3]]);
}

\\ Alternate method to compute the babystock, just using one huge ass qfminim, unit rank 1 only
one_massive_qfminim(G, giant_sublattice, S_radius)={
    my(x, LLL_reduced_yu, gram_mat, bound, bigqfminim,vecholder,new_y,vec_numerical,epB);
    x = embed_real(G, G[5][1]);
    LLL_reduced_yu = x*qflll(x);
    gram_mat=LLL_reduced_yu~*LLL_reduced_yu;                                    \\ get the gram matrix

    bound = create_target(G, giant_sublattice[,1]/6);
    print("Containg sphere bound1: ", precision(norml2(bound),10));

    tempvec = zero_vec; tempvec[1] = S_radius;
    tempvec = create_target(G, tempvec~);
    tempvec = norml2(tempvec);
    print(precision(tempvec,10));
    bound = norml2(bound)+tempvec;
    \\bound = norml2(bound);
    print("Containg sphere bound +S: ", precision(bound,10));

    bigqfminim = qfminim(gram_mat, bound/6,, 2)[3];
    babystock1 = Map();
    for(ii=1, length(bigqfminim),
        vecholder=bigqfminim[,ii];                                                        \\ get the ith element w_i
        new_y=idealdiv(G,matid(field_deg),vecholder);                                     \\ get the ideal (1/w_i)*y
        if(checkred_old(new_y,G,eps)==1,                                                  \\ check if its reduced, and if yes, add info_vec to eB
            vec_numerical = (G[5][1]*vecholder)~;                                         \\ get the numericals
            epB=Set(concat(epB,[[new_y,vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1])),eps)]]));
            if(mapisdefined(babystock1, new_y),
                mapput(babystock1, new_y, concat( mapget(babystock1, new_y), [vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1])),eps)]) );,
                mapput(babystock1, new_y, [vector_approximate(log(abs(vec_numerical[1..G.r1+G.r2-1])),eps)]);
            );
        );
    );

    return(babystock1);
}

\\ INPUT:
\\ - G is a number field
\\ - lattice_lambda is a full rank sublattice of the log-unit lattice
\\ - B indicates how much of the fundamental region of lattice_lambda we want to search
\\   it is a positive integer, and we search the region F/B, where F is the fundamental region
\\ - delta indicates the proportion to allocate for the babystock
\\ - eps is the error
\\ - alg can be used to change the babystock search algorithm, So far only two options
\\ - lattice scan (default), and neighbours ("NEIGHBOURS")
bsgs(G, cpct_reps, B, delta, eps, REQ_BSGS,FILE1="data/tmp-bsgs-log", alg="SCAN")={

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
        babystock,
        giant_divisor,
        matches,
        new_vec,lattice_lambda
    );
    lattice_lambda = log_lattice_from_compact_set(G, cpct_reps);
    if(alg != "NEIGHBOURS" && alg != "SCAN", print("Unknown algorithm specified, returning -1", return(-1) ));

    S_radius = (sqrt(poldegree(G.pol))/4)*log(abs(G.disc));                     \\ get radius S (may not need)
    if(DEBUG_BSGS>0 , print("Bound B = ", B, "\nradius S = ", precision(S_radius,10) ); );

    [avec, giant_sublattice, babystock_region] = get_giant_step_params(G,lattice_lambda, r, B, delta,REQ_BSGS);
    babystock_region = expand_babystock_region(babystock_region,S_radius);

    tb = getabstime();
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ THIS SECTION USES NEIGHBORS TO SEARCH THE BABYSTOCK
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(alg == "NEIGHBOURS",
        my(mid_ideal, nu, logdist);
        print("Using neighbors for babystock\n",precision(babystock_region,10));
        [mid_ideal, nu, logdist] = get_initial_neighbour(G, giant_sublattice, field_deg,eps);

        Mset2= COLLECT_BABYSTOCK(G, mid_ideal, logdist[1..length(lattice_lambda)], babystock_region, eps);
        [lattice_lambda, newctr] = check_babystock_for_units(Mset2,lattice_lambda,G,eps);

        babystock = Map();
        for(ctr=1, length(Mset2),
            if(mapisdefined(babystock, Mset2[ctr][1]),
                existing_entry = mapget(babystock, Mset2[ctr][1]);
                repeatflag = is_repeat_babystock(existing_entry, Mset2[ctr][2],eps);

                if(repeatflag==0, mapput(babystock, Mset2[ctr][1], concat( mapget(babystock, Mset2[ctr][1]), [Mset2[ctr][2]]) );),
                mapput(babystock, Mset2[ctr][1], [Mset2[ctr][2]]);
            );
        );

    );\\ end neighbors for babystock
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\
    \\ THIS SECTION SEARCHES USING THE SCAN ALGORITHM
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if(alg == "SCAN",
        my(foundflag = 1, num_elts_found = []);
        [lattice_lambda,babystock, num_elts_found] = babystock_scan_jump(matid(field_deg),lattice_lambda, giant_sublattice, G, eps);
        \\print("unscanned_babystock_region variable: ", num_elts_found);
        foundflag = length(num_elts_found);


        while(foundflag !=0,
            \\ this should regenerate lattice_lambda to remove any precision loss
            lattice_lambda = log_lattice_from_compact_set(G, cpct_from_loglattice(G, lattice_lambda, eps));
            [avec, giant_sublattice, babystock_region] = get_giant_step_params(G,lattice_lambda, r, B, delta,REQ_BSGS);
            babystock_region = expand_babystock_region(babystock_region,S_radius );

            [lattice_lambda, babystock, num_elts_found] = babystock_scan_jump(matid(field_deg),lattice_lambda, giant_sublattice, G, eps);
            foundflag = length(num_elts_found);
            if(DEBUG_BSGS>0 ,print("FoundFlag ",foundflag););
        );
    );

    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    babytime = getabstime()-tb;
    write(FILE1, "Babystock time: ", precision(babytime, 10), "   ",  precision(babytime/60000.0, 15));

    print("\nGiant Loop\n");
    for(i=1, length(avec), avec[i] += 1);                                       \\ adjust avec to include boundaries.

    \\normproduct = 1;
    \\for(i=1, length(giant_sublattice), normproduct*=norml2(giant_sublattice[,i]) );
    \\print("normproduct  ", precision(normproduct,10));
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\ Jump method for computation of the giant steps

        tg = getabstime();
        lattice_lambda = jump_giant_steps(G, lattice_lambda, giant_sublattice, babystock, avec, eps);
        gianttime = getabstime() -tg;
    bsgs_end = getabstime();
    write(FILE1,  "\nGiantstep time:  ", gianttime,  "   ",precision(gianttime/60000.0,15), "\nTotal: ",precision((bsgs_end-bsgs_start)/60000.0 ,15 ));

    return(cpct_from_loglattice(G,lattice_lambda,eps));
}
