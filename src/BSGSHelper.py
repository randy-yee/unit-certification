read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/CompactRepresentation.py");
read("src/Neighbours.py")
read("src/bounds.gp")

/*
Helper functions for the main BSGS routine

compute_vector_norm
selectBStockSideLengths
increment_with_place_marker
update_directions
is_repeat_babystock
compute_max_point_distance
get_web_step_multipliers
get_initial_webpoint
*/

compute_vector_norm(i, lattice_lambda, unit_rank, B)=
{
    my(vecnorm);
    vecnorm =  max(1,floor( sqrt(norml2(lattice_lambda[,i] )) ));
    if (i == unit_rank, vecnorm = vecnorm/B);
    return(vecnorm);
}

\\# function intended to choose the a_i for defining the babystock region
selectBStockSideLengths(lattice_lambda, fullregion_to_babystock_ratio, hybridB)=
{
    my(
        avec,
        dimensions = length(lattice_lambda),
        norm_vr = sqrt(norml2(lattice_lambda[,dimensions] ));
    );

    if(dimensions == 1,
        avec = vector(dimensions, i , ceil(fullregion_to_babystock_ratio));
    , \\else
        if (norm_vr/hybridB > 2, dimensions++);
        sides = max(1,floor(sqrtn(fullregion_to_babystock_ratio, dimensions-1)));             \\ approximate the sides with integers

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

        if(DEBUG_BSGS, print("product without  a_r ",ai_product); );
        if (dimensions !=unit_rank+1,
            avec = concat(avec ,  min(max(1, round(norm_vr/B)  ), round(fullregion_to_babystock_ratio/ai_product ))  );
        );
    );
    return(avec);
}


\\ modifies the input current_vec in place by incrementing it by 1
\\ INPUT:
\\ - A vector of integers A which defines the subdivision for the babystock
\\ - A vector of integers V that defines a particular 'giant step'
\\   for all i, V[i] < A[i]
\\ OUTPUT:
\\ - place marker of position that was modified
increment_with_place_marker(~a_vec, ~current_vec)={
    my(
        place = length(current_vec),
        carryflag = 0;
    );
    current_vec[place]+=1;
    if(current_vec[place] >= a_vec[place], carryflag =1);

    while(carryflag == 1,
        current_vec[place] = 0; carryflag = 0;
        if(place == 1, return(place));
        place -= 1;
        current_vec[place]+=1;
        if(current_vec[place] >= a_vec[place], carryflag =1);
    );

    return(place);
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

\\# computes the length of the longest vector
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

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
which contains the box generated by the columns of basismat
\\ represented by two points in R^r, where r is the unit rank
get_axis_aligned_box(basismat)={
    my(r = length(basismat), twovec, maxcoords, mincoords, vertex_coeffs, vertex );
    twovec = vector(length(basismat), i, 2);
    zerovec = vector(r, i, 0);
    vertex_coeffs = zerovec;
    maxcoords = vertex_coeffs;
    mincoords = vertex_coeffs;
    increment_coordinates(twovec, vertex_coeffs);
    while(vertex_coeffs != zerovec,

        vertex = basismat*(vertex_coeffs~);
        for(i=1, r,
            if(vertex[i] < mincoords[i], mincoords[i] = vertex[i]);
            if(vertex[i] > maxcoords[i], maxcoords[i] = vertex[i]);
        );
        increment_coordinates(twovec, vertex_coeffs);
    );
    return([mincoords, maxcoords]);
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



\\ subalgorithm for bsgs, adjust by the radius S
expand_babystock_region(babystock_region,S_radius )={
    my(S_radius_vector);
    if(DEBUG_BSGS>1 ,print("Babystock box: ", precision(babystock_region,10)););
    S_radius_vector = vector(length(babystock_region[1]), i , S_radius);
    babystock_region[1] = babystock_region[1] - S_radius_vector;
    babystock_region[2] = babystock_region[2] + S_radius_vector;

    return(babystock_region);
}
