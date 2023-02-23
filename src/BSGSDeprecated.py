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
