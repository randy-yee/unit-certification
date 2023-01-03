read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/bounds.gp")

/*
get_giant_step_params

*/

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
    print("Old babystock volume value:  ", precision(smallsquare,10));
    if(DEBUG_BSGS, print("Default babystock region size: ", precision(smallsquare,10) ););

    \\babystock_scale_factor =((2)^(unit_rank-1))*(log(abs(G.disc))/8)^(1/(2));              \\ increase this to make the babystock smaller
    \\babystock_scale_factor = (2^unit_rank) * log(abs(G.disc)) / 32;
    print("BSGS: Babystock scaling factor ", precision(babystock_scale_factor,10));

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

non_integral_subdivisions(G, lattice_lambda, dimensions, detLambda, B, babystock_scale_factor, REQ_BSGS)={
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
    \\ This is the basic calculation for the babystock region
    smallsquare = sqrt(  (abs(detLambda)/B)*g_n/b_n  );
    print("Old babystock volume value:  ", precision(smallsquare,10));
    if(DEBUG_BSGS, print("Default babystock region size: ", precision(smallsquare,10) ););

    \\babystock_scale_factor =((2)^(unit_rank-1))*(log(abs(G.disc))/8)^(1/(2));              \\ increase this to make the babystock smaller
    \\babystock_scale_factor = (2^unit_rank) * log(abs(G.disc)) / 32;
    print("BSGS: Babystock scaling factor ", precision(babystock_scale_factor,10));

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

    norm_vr = sqrt(norml2(lattice_lambda[,dimensions] ));

    if(dimensions == 1,
        avec = vector(dimensions, i , fullregion_to_babystock_ratio);
    , \\else
        if (norm_vr/B > 2, dimensions++);
        sides = max(1,sqrtn(fullregion_to_babystock_ratio, dimensions-1));             \\ approximate the sides with integers

        if(DEBUG_BSGS,print("(r-1)th root of target area  ", precision(sides,10)););
        avec = [];
        my(diffvector = [],
            vecnorm =1,
            ai_product=1,
            areadiff = 0,
            partialproduct = 1);

        \\# adjust in case the r-1 th root is larger than the vector's norm
        for(i=1, dimensions-1,
            vecnorm =  sqrt(norml2(lattice_lambda[,i] ));
            if (i == dimensions-1, vecnorm/=B);
            avec = concat(avec, min(vecnorm, sides));
            diffvector = concat(diffvector, vecnorm - sides );
            ai_product *= avec[i];
        );
        for(i=1, length(diffvector),
            if(diffvector[i] > 0,
                partialproduct = (ai_product/avec[i]);
                areadiff = (fullregion_to_babystock_ratio - ai_product)/partialproduct;
                if(areadiff >= 1,
                    avec[i] += min(areadiff, diffvector[i]);
                    ai_product += partialproduct*min(areadiff, diffvector[i]) );
            );
        );

        if(DEBUG_BSGS, print("product without  a_r ",precision(ai_product,10)); );

        if (dimensions !=unit_rank+1,
            avec = concat(avec ,  min(norm_vr/B , fullregion_to_babystock_ratio/ai_product )  );
        );
    );

    finalproduct =1; for(i=1, length(avec), finalproduct*=avec[i]);
    for(i=1, length(avec),
        if(avec[i] == 0, avec[i]+=1);
    );

    return(avec);
}


\\ subroutine of the bsgs method.
\\ Determines a basis for the giant steps (coefficients limited by avec)
\\ INPUT:
\\ - G is the number field
\\ - lattice_lambda is the current log lattice (r x r)
\\ - r is the unit rank
\\ - B defines the fraction of the fundamental region F to search, 1/B
\\ - babystock_scale_factor defines the volume of the babystock
\\ - REQ_BSGS is the precision used for comparison
\\ OUTPUT:
\\ - avec is a r-vector of coefficient upperbounds for the giant steps
\\ - giant_legs is a basis that determines the giant steps
\\ - babystock_region is an axis aligned box that contains the babystock region. Given by two point vectors.
get_giant_step_params(G, lattice_lambda, r, B, babystock_scale_factor, REQ_BSGS)={
    my(det_lambda, giant_legs, extra_log_coords, babystock_region);
    det_lambda = unscaled_determinant(G, lattice_lambda);

    originalarea = 1;
    if(1,
        for(i=1, length(lattice_lambda),
            print("- length of v_",i,"   ",  precision(sqrt(norml2(lattice_lambda[,i])),10));
            originalarea*=sqrt(norml2(lattice_lambda[,i]))
        );
        print("- product: ", precision(originalarea,10), "\n- det_lambda: ", precision(det_lambda,10) );
    );

    \\\# determine r independent vectors that will define the giant steps
    \\\# note that the fundamental region of these vectors is the babystock
    \\avec = get_subdivisions(G,lattice_lambda,r, det_lambda, B, babystock_scale_factor, REQ_BSGS);
    print("Using non-integer values for a-vector");
    avec = non_integral_subdivisions(G,lattice_lambda,r, det_lambda, B, babystock_scale_factor, REQ_BSGS);
    print(precision(non_integral_subdivisions(G,lattice_lambda,r, det_lambda, B, babystock_scale_factor, REQ_BSGS),10));

    extra_log_coords = vector(r, i, extra_log_coordinate(G.r1, G.r2, lattice_lambda[,i]));
    giant_legs = matconcat([lattice_lambda; extra_log_coords]);
    for (i=1, r-1, giant_legs[,i] = -giant_legs[,i]/avec[i]);
    giant_legs[,r] = -giant_legs[,r]/(avec[r]*B);

    \\giant_legs = matrix(r, r-1, s,t, -lattice_lambda[s,t]/avec[t]);             \\ giant_legs is a matrix with columns (-v_i /a_i) for i =1 ..r-1
    \\giant_legs = matconcat([giant_legs, -lattice_lambda[,r]/(avec[r]*B)] );     \\ column r is (-v_r /(a_r*B))

    if(1,
        \\print("- avec ", avec, "  ", B);
        area = 1;
        for(i=1, length(giant_legs),
            \\print("Norms of scaled vectors (a_i, B) ", precision(sqrt(norml2(giant_legs[,i][1..r])),10));
            area*=sqrt(norml2(giant_legs[,i][1..r])));
        print("actual babystock area region (multiplying vector norms): ", precision(area,10),"\n \n");
    );

    return([avec,giant_legs]);
}
