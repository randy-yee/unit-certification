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


\\ old version that still contains code for using lll based selection of beta.
\\
reddiv_compact_old(~y,~u,~G,~M1, p_avoid=1)={
    my(y1, ideal_uY, numerical_mat_Y, red1, shortest_vec, nu, lmu,
        ideal_denom,vec_ctr,beta_found =0,
        comp = 2^(-ceil((poldegree(G.pol)^2 +2)*log(infinity_norm(u))+2*poldegree(G.pol)^2 +5))
    );

    numerical_mat_Y = M1*y;                                                     \\ complex embedding matrix of y
    ideal_uY = mulvec(~numerical_mat_Y, ~u);                                    \\ ideal u*y
    LLL_change_of_basis = get_LLL_basis_change(G, ideal_uY);                    \\ qflll output has columns which are coords wrt the input matrix NOT the integral basis
    LLL_numerical_uY = ideal_uY*LLL_change_of_basis;                            \\ Obtain LLL basis of u*y in numerical form (possibly complex)
    LLLcoeffmat = y*LLL_change_of_basis;                                        \\ LLL basis, coords wrt the integral basis

    beta= LLLcoeffmat[,1];                                                      \\ beta holds coordinates of mu wrt the integral basis
    \\ need to scan to make sure the first basis vector is a shortest one
    real_mat_uY = embed_real(~G, ~LLL_numerical_uY);
    enum_result = qfminim(real_mat_uY~*real_mat_uY, norml2(real_mat_uY[,1])-comp,,2);

    \\\ new part, where beta is selected using qfminim\'s shortest vector function
    true_shortest = qfminim(real_mat_uY~*real_mat_uY,,,2);
    beta = LLLcoeffmat*true_shortest[3][,1];
    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);                      \\ the reduced ideal y / mu, in hnf form
    nu = pointwise_vector_mul(abs(M1*beta),u)~;
    \\GP_ASSERT_VEC_NEAR(nu,abs(shortest_vec), comp  );
    lmu = log(nu)-log(u);
    return([new_y, nu, lmu, beta]);

    \\boolA = (enum_result[1]!=2 && !is_minimum(LLLcoeffmat,beta , G, comp));
    boolB = (length(enum_result[3])!=0 && !is_minimum(LLLcoeffmat,beta , G, comp));

    if(boolB,
        short_index =1;
        short_length = norml2(real_mat_uY*enum_result[3][,1]);
        for(j=1, length(enum_result[3]),
            iter_length = norml2(real_mat_uY*enum_result[3][,j]);
            if(iter_length < short_length,
                short_index = j; short_length = iter_length;
            );
        );
        beta = LLLcoeffmat*enum_result[3][,short_index];
        if(!is_minimum(LLLcoeffmat, beta, G, comp),
            print("elements found in normed body of supposed minimum!!!!");
            breakpoint()
        );
        shortest_vec = LLL_numerical_uY*enum_result[3][,short_index];
    , \\else
        shortest_vec = LLL_numerical_uY[,1];                                        \\ vector of complex embeddings for the shortest vector of u*y, denoted mu

    );

    new_y = idealdiv(G,y,beta); new_y = idealhnf(G,new_y);                      \\ the reduced ideal y / mu, in hnf form

    \\\#Ran into trouble using abs(shortest_vec) with massive precision loss
    \\\# instead use alternate formula u*psi(beta)
    \\nu=abs(shortest_vec)~;                                                      \\ nu is a t_VEC of dimension r, (complex coordinates are not squared)
    nu = pointwise_vector_mul(abs(M1*beta),u)~;
    \\GP_ASSERT_VEC_NEAR(nu,abs(shortest_vec), comp  );
    lmu = log(nu)-log(u);                                                       \\ expect equal to log(nfeltembed(G, beta))
    \\GP_ASSERT_VEC_NEAR(lmu,log(nfeltembed(G, beta) ),10);

    return([new_y, nu, lmu, beta]);                                                     \\
}



cpct_rep_final_enum(G, idealB, beta, log_beta, desired_betalog, eps, testFlag = 0)={
    if (type(testFlag) == "t_INT" && (testFlag == 1),
        print("final enum");
    );

    my(neighbours_output, boundary_vector, ctr=1, unit_rank = G.r1+G.r2-1,
        exp_boundary, latticeB, lll_basis_change_matrix, latticeB_lll, scan_elements);
    GP_ASSERT_EQ(length(log_beta), G.r1+G.r2);

    boundary_vector = (desired_betalog - log_beta);
    print("Fix this! ", boundary_vector);breakpoint();

    \\breakpoint(desired_betalog - log_beta);
    \\#setup ideal for qfminim
    alternate = 1;
    temp_precision = ceil(idealPrecision(G, idealB, normlp(boundary_vector)));
    oldbitprecision = change_precision(temp_precision);
    if(!alternate,
        exp_bvec = exponentiate_logvec(G.r1+G.r2-1, boundary_vector);
        exp_boundary = norml2(exp_bvec); print(precision(exp_boundary,10));
        latticeB = embed_real(G, G[5][1]*idealB);
        lll_basis_change_matrix = qflll(latticeB);
        latticeB_lll = latticeB*lll_basis_change_matrix;
        lll_ideal = idealB*lll_basis_change_matrix;
        scan_elements = qfminim((latticeB_lll~*latticeB_lll),exp_boundary^2+eps,,2)[3]; \\#check all elements within the bound
        if (samevecs(desired_betalog, log_beta, eps),  return([idealB, log_beta, beta] ) );
        my(check_beta, checkvec);

        for(i=1, length(scan_elements),

            check_beta = lll_ideal*scan_elements[,i];          \\# beta wrt integral basis
            checkvec = log_beta + log(abs( nfeltembed(G, check_beta)));

            if(samevecs(desired_betalog, checkvec, eps),
                change_precision(oldbitprecision);
                \\print(i, "  ", precision(checkvec,10));
                idealB = idealdiv(G, idealB, check_beta);
                return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
            );
        );
    );
    if(alternate,
        degree = poldegree(G.pol);
        print("boundary: ", precision(boundary_vector,10));
        print(precision(exponentiate_logvec(G.r1+G.r2-1, boundary_vector, 1),10));
        exp_bvec = exponentiate_logvec(G.r1+G.r2-1, boundary_vector, 1);
        if(matsize(G[5][1])[1] != G.r1+G.r2, print("final_enum vector length mismatch"); breakpoint(); );
        scaled_latticeB = mulvec(G[5][1]*idealB, exp_bvec);
        scaled_latticeB = embed_real(G, scaled_latticeB);

        lll_CoB = qflll(scaled_latticeB);
        lll_ideal = idealB*lll_CoB;

        scaled_lll_lattice = scaled_latticeB*lll_CoB;

        \\print("degree bound: ", degree+eps); for(s = 1, length(scaled_lll_lattice), print(precision(norml2(scaled_lll_lattice[,s]),10)));breakpoint();
        scan_elements = qfminim(scaled_lll_lattice~*scaled_lll_lattice,degree+eps*1.0,,2)[3];
        print(length(scan_elements));
        for(i=1, length(scan_elements),
            check_beta = lll_ideal*scan_elements[,i];          \\# beta wrt integral basis
            checkvec = log_beta + log(abs( nfeltembed(G, check_beta)));

            if(samevecs(desired_betalog, checkvec, eps),
                idealB = idealdiv(G, idealB, check_beta);
                change_precision(oldbitprecision);
                return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
            );
            checkvec = log_beta - log(abs( nfeltembed(G, check_beta)));
            if(samevecs(desired_betalog, checkvec, eps),
                idealB = idealmul(G, idealB, check_beta);
                change_precision(oldbitprecision);
                return([idealB, checkvec, nfeltmul(G,beta,check_beta)]);
            );
        );
    );
    print("No elements satisfy the condition. Error.");
    return(-1);
}

\\\ Use neighbor algorithm to search the babystock
neighbour_search(G,lattice_lambda, babystock_region, ~babystock, giant_sublattice, eps)=
{
    my(mid_ideal, nu, logdist);
    print("Using neighbors for babystock\n",precision(babystock_region,10));
    [mid_ideal, nu, logdist] = get_initial_neighbour(G, giant_sublattice, field_deg,eps);

    Mset2= COLLECT_BABYSTOCK(G, mid_ideal, logdist[1..length(lattice_lambda)], babystock_region, eps);
    [lattice_lambda, newctr] = check_babystock_for_units(Mset2,lattice_lambda,G,eps);

    for(ctr=1, length(Mset2),
        if(mapisdefined(babystock, Mset2[ctr][1]),
            existing_entry = mapget(babystock, Mset2[ctr][1]);
            repeatflag = is_repeat_babystock(existing_entry, Mset2[ctr][2],eps);

            if(repeatflag==0, mapput(babystock, Mset2[ctr][1], concat( mapget(babystock, Mset2[ctr][1]), [Mset2[ctr][2]]) );),
            mapput(babystock, Mset2[ctr][1], [Mset2[ctr][2]]);
        );
    );
}

\\ subalgorithm in the bsgs method, get the initial neighbour to start
\\ the babystock search when using the neighbours method
get_initial_neighbour(G, giant_sublattice, field_deg,eps)={
    my(middle_babystock);
    middle_babystock = vector(length(giant_sublattice), i, 1/2);
    middle_babystock = giant_sublattice*(middle_babystock~);
    return(giantstep(matid(field_deg), middle_babystock, G, field_deg, eps);)
}



\\ Determines the size of the babystock area
\\ INPUT:
\\    dimensions, an integer equal to r, the rank of the lattice (in FJ, this is r1+r2-1)
\\    detLambda, the determininant of the input lattice Lambda
\\    B the prime search bound.
\\    babystock_scale_factor is used as a coefficient to scale the babystock region.
\\    A larger number means a smaller babystock region
\\ OUTPUT:
\\ - a vector of integers of length = length(Lambda)
get_subdivisions(G, lattice_lambda, dimensions, detLambda, B, babystock_scale_factor, REQ_BSGS)={
    my(smallsquare,
        sides,
        avec,
        a_product,
        g_n, b_n,
        ind=1,
        unit_rank = G.r1 + G.r2-1,
        deg = poldegree(G.pol),
        fullregion_to_babystock_ratio= 1);

    g_n = giant_n( deg, log(abs(G.disc)), REQ_BSGS, real(log(detLambda)) );
    b_n =  baby_n( deg, log(abs(G.disc)), REQ_BSGS, real(log(detLambda)) );
    \\# This is the basic calculation for the babystock region
    smallsquare = sqrt(  (abs(detLambda)/B)*g_n/b_n  );
    \\print("Old babystock volume value:  ", precision(smallsquare,10));

    \\babystock_scale_factor =((2)^(unit_rank-1))*(log(abs(G.disc))/8)^(1/(2));              \\ increase this to make the babystock smaller
    \\babystock_scale_factor = (2^unit_rank) * log(abs(G.disc)) / 32;
    print("BSGS: Babystock scaling factor ", precision(babystock_scale_factor,10), " B ", B);

    \\# Compute ratio of (full search region / babystock region)
    \\# any modification of the scale factor will shrink the babystock region
    \\fullregion_to_babystock_ratio = (babystock_scale_factor)*abs(detLambda)/(smallsquare*B);

    \\ assigning a volume of the babystock
    smallsquare = babystock_scale_factor;
    fullregion_to_babystock_ratio = abs(detLambda)/(smallsquare*B);

    write(OUTFILE1, " babystock_vol = ", precision( smallsquare,10),".  Region/babystock ratio: ", strprintf("%10.5f",fullregion_to_babystock_ratio) );

    if(DEBUG_BSGS,
        print("G_n, B_n  ",ceil(g_n), "   ", ceil(b_n));
        print("a_I product target: ", ceil(fullregion_to_babystock_ratio), "   " );
    );

    norm_vr = ceil( sqrt(norml2(lattice_lambda[,dimensions] )) );
    if(DEBUG_BSGS,print("norm of rth vector: ", precision(norm_vr, 10)););

    if(dimensions == 1,
        avec = vector(dimensions, i , ceil(fullregion_to_babystock_ratio));
    , \\else
        if (norm_vr/B > 2, dimensions++);
        sides = max(1,floor(sqrtn(fullregion_to_babystock_ratio, dimensions-1)));             \\ approximate the sides with integers

        if(DEBUG_BSGS,print("(r-1)th root of target area  ", precision(sides,10)););
        avec = [];
        my(diffvector = [],
            vecnorm =1,
            ai_product=1,
            areadiff = 0,
            partialproduct = 1);

        \\# adjust in case the r-1 th root is larger than the vector's norm
        for(i=1, dimensions-1,
            vecnorm =  compute_vector_norm(i, lattice_lambda, unit_rank, B);
            avec = concat(avec, min(vecnorm, sides));
            diffvector = concat(diffvector, vecnorm - sides );
            ai_product *= avec[i];
        );
        for(i=1, length(diffvector),
            if(diffvector[i] > 0,
                partialproduct = (ai_product/avec[i]);
                areadiff = round((fullregion_to_babystock_ratio - ai_product)/(ai_product/avec[i]));
                if(areadiff >= 1,
                    avec[i] += min(areadiff, diffvector[i]);
                    ai_product += partialproduct*min(areadiff, diffvector[i]) );
            );
        );

        if(DEBUG_BSGS, print("product without  a_r ", precision(ai_product,10)); );

        if (dimensions !=unit_rank+1,
            print("exceptional condition");
            avec = concat(avec ,  min(max(1, round(norm_vr/B)  ), round(fullregion_to_babystock_ratio/ai_product ))  );
        );
    );

    finalproduct =1; for(i=1, length(avec), finalproduct*=avec[i]);
    for(i=1, length(avec),
        if(avec[i] == 0, avec[i]+=1);
    );
    print("BSGS: Babystock a_i=",avec, ". a_i product=", finalproduct, " target=", round(fullregion_to_babystock_ratio));

    return(avec);
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
            trackingLog += embeddings_to_normalized_logvec(G, nfeltembed(G, compactTracking[length(compactTracking)][1]));
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

    \\# go through each unique ideal and enumerate once. The results are multipled
    \\# by each associated u to account for all the distinct minima
    print("Babystock: enumerating");
    if (length(scanIdeals) < 1, return([lattice_lambda, []]); );

    scanIdealsMatrix = Mat(scanIdeals)[,1];
    print("Number of ideal to scan: ", length(scanIdealsMatrix));
    for(i = 1, length(scanIdealsMatrix),
        distanceList = mapget(scanIdeals, scanIdealsMatrix[i]);

        \\cpct_list = mapget(idealCompactGenerator, scanIdealsMatrix[i]);
        \\GP_ASSERT_EQ(length(distanceList)-1, length(cpct_list));
        if (length(distanceList) > 2, print("Collision found within babystocks"););
        nu = distanceList[1];
        \\overlap_scanball_DEBUG(~G, ~baby_hashmap, ~scanIdealsMatrix[i], ~nu, ~distanceList, scan_radius, eps, ~repeated_minima); breakpoint();
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
