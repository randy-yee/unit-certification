read("src/VectorMethods.py");
read("src/LogarithmLattice.py");
read("src/bounds.gp")

/*
Subfunctions for the purpose of getting the giant step parameters
- get_giant_step_params

*/

\\ Determines the size of the babystock area
\\ INPUT:
\\    dimensions, an integer equal to r, the rank of the lattice (in FJ, this is r1+r2-1)
\\    detLambda, the determininant of the input lattice Lambda
\\    B the prime search bound.
\\    babystock_scale_factor sets the size of the bstock
\\    REQ_BSGS is the worst case computed precision required for BSGS
\\ OUTPUT:
\\ - a vector of integers of length = length(Lambda)
\\\# dimensions = unit rank
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

    \\# This is the basic calculation for the babystock region
    \\# It is currently overriden
    \\g_n = giant_n( deg, log(abs(G.disc)), REQ_BSGS, real(log(detLambda)) );
    \\b_n =  baby_n( deg, log(abs(G.disc)), REQ_BSGS, real(log(detLambda)) );
    \\smallsquare = sqrt(  (abs(detLambda)/B)*g_n/b_n  );
    \\print("Old babystock volume value:  ", precision(smallsquare,10));

    \\# Compute ratio of (full search region / babystock region)
    \\# any modification of the scale factor will shrink the babystock region
    \\#fullregion_to_babystock_ratio = (babystock_scale_factor)*abs(detLambda)/(smallsquare*B);

    \\# assigning a volume of the babystock. This is an input to bsgs
    smallsquare = babystock_scale_factor;

    \\# the product of the a_i should be roughly equal to the value below
    fullregion_to_babystock_ratio = abs(detLambda)/(smallsquare*B);
    print("babystock volume: ", precision(smallsquare,10), " B = ", floor(B) , ".  ratio: ", precision(fullregion_to_babystock_ratio,10));
    write(OUTFILE1, " babystock_vol = ", precision( smallsquare,10),".  Region/babystock ratio: ", strprintf("%10.5f",fullregion_to_babystock_ratio) );

    GP_ASSERT_TRUE(dimensions == unit_rank);
    GP_ASSERT_NEAR(detLambda , get_abs_determinant(lattice_lambda), 0.00001);

    norm_vr = sqrt(norml2(lattice_lambda[,dimensions] ));


    \\# product of the a_i should equal (size of search region)/(size of babystock region)
    \\# calculation below aims to achieve this
    if(unit_rank == 1,
        avec = vector(unit_rank, i , fullregion_to_babystock_ratio); \\#easy case
    , \\else
        \\# when this if fails, we exclude the last vector from the calculation
        \\# and set a_r = 1
        if(norm_vr/B > 1, dimensions++);

        sides = max(1,sqrtn(fullregion_to_babystock_ratio, dimensions-1));

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

            if (i == unit_rank, vecnorm/=B);
            print(i, " ", precision(vecnorm,10), "  ", precision(sides,10));
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

        \\# this will only be the case if (norm_vr/B <= 1), which meant we
        \\# excluded a_r from the calculation and set it to be 1
        if (dimensions !=unit_rank+1,
            print("Final a_i: ", precision(norm_vr/B,10), "   ",precision(fullregion_to_babystock_ratio/ai_product,10));
            avec = concat(avec ,  1  );
            print("fullregion_to_babystock_ratio/(a_i product) : ", precision(fullregion_to_babystock_ratio/ai_product,10));
            \\concat(avec ,  min(norm_vr/B , fullregion_to_babystock_ratio/ai_product )  );
        );
    );

    finalproduct =1; for(i=1, length(avec), finalproduct*=avec[i]);
    for(i=1, length(avec),
        if(avec[i] == 0, avec[i]=1);
        if(avec[i] < 1, print("warning, the a_i is less than one for i = ",i));
    );

    miniLambda = lattice_lambda;
    for(i=1, length(avec),
        miniLambda[,i] /= avec[i]);
    miniLambda[,unit_rank]/=B;
    print(precision(fullregion_to_babystock_ratio, 10), "  ", precision(finalproduct,10),"\nbabystock det: ", precision(abs(matdet(miniLambda)),10));
    return(avec);
}


\\ subroutine of the bsgs method.
\\ Determines a basis for the giant steps (coefficients limited by avec)
\\ giant_leg vectors are extended to vectors of length (r+1). The extra coord
\\ is computed so that the sum of the logs of each vector equals 0
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
    det_lambda = get_abs_determinant(lattice_lambda);

    originalarea = 1;
    if(1,
        for(i=1, length(lattice_lambda),
            print("- length of v_",i,"   ",  precision(sqrt(norml2(lattice_lambda[,i])),10));
            originalarea*=sqrt(norml2(lattice_lambda[,i]))
        );
        print("- product of v_i: ", precision(originalarea,10), ", - det_lambda: ", precision(det_lambda,10) );
    );

    \\\# determine r independent vectors that will define the giant steps
    \\\# note that the fundamental region of these vectors is the babystock

    avec = non_integral_subdivisions(G,lattice_lambda,r, det_lambda, B, babystock_scale_factor, REQ_BSGS);
    print("Using non-integer values for a-vector: ", precision(avec,10));

    extra_log_coords = vector(r, i, extra_log_coordinate(G.r1, G.r2, lattice_lambda[,i]));
    giant_legs = matconcat([lattice_lambda; extra_log_coords]);
    for (i=1, r-1, giant_legs[,i] = -giant_legs[,i]/avec[i]);
    giant_legs[,r] = -giant_legs[,r]/(avec[r]*B);

    if(1,
        print_area_info(giant_legs, r);
    );
    print("giant leg det: ", precision(abs(matdet(giant_legs[1..r,])),10));
    return([avec,giant_legs]);
}

print_area_info(~giant_legs, rank)=
{
    my(area = 1);
    for(i=1, length(giant_legs),
        print("Norms of scaled vectors (a_i, B) ", precision(sqrt(norml2(giant_legs[,i][1..rank])),10));
        area*=sqrt(norml2(giant_legs[,i][1..rank])));
    \\print("giant leg det: ", precision(abs(matdet(giant_legs[1..rank,])),10));
    print("actual babystock area region (multiplying vector norms): ", precision(area,10),"\n \n");
}
