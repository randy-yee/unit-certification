read("src/BabyStepGiantStep.py")
read("src/CompactRepresentation.py");
read("src/bounds.gp");
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);
DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;
DEBUG_BSGS = 0;
\\\\ Setup functions
setInstanceVariables(readData)={
  K = nfinit(readData[1]);
  lglat = process_complex_loglattice(K ,readData[3]);
  reg1 = unscaled_determinant(K, lglat);
  r = K.r1+K.r2 -1;
  return([K, lglat, reg1,r]);
}

compute_precision(~K, ~lglat, ~reg1)={
    my(X1,X2,sumv = lglat[,1]);
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K.pol), log(abs(K.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K.pol), log(abs(K.disc)),abs(reg1),infinity_norm(sumv) );
    print("BSGS precision ");print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));

    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    REQ_COMPARE = ceil((poldegree(K.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);
    return([REQ_BSGS, REQ_COMPARE,eps]);
}

scale_lattice_column(loglat, col, factor)={
    print("Scaling lattice column ", col, " by ", factor);
    copy_lattice = loglat;
    copy_lattice[,col] = factor*copy_lattice[,col];
    return(copy_lattice);
}

\\\
generateFileStrings(signature_string, suffixString, auxilliary)=
{
    outfilestring = strexpand("data/data-bsgs-",signature_string,suffixString);

    infilestring = concat(concat("input/test-poly-",signature_string),".gp");
    \\infilestring = concat(concat("input/polynomial-",signature_string),".gp");
    OUTFILE1 = outfilestring;

    if(length(aux) >1 && (type(aux[2]) == "t_STR"),
        OUTFILE1 = aux[2];

    );
    if(length(aux) >2 && (type(aux[3]) == "t_STR"),
        infilestring = aux[3];
    );
    print("Output directed to file ", OUTFILE1);
    print("Input File is: ",infilestring );
    return([OUTFILE1, infilestring]);
}

generateFileStringsGeneral(signature_string, experimentString, suffixString, aux)=
{
    outfilestring = strexpand("data/data-",experimentString, "-",signature_string,suffixString);
    print("Output directed to file ", outfilestring);
    \\infilestring = concat(concat("input/test-poly-",signature_string),".gp");
    infilestring = concat(concat("input/polynomial-",signature_string),".gp");
    OUTFILE1 = outfilestring;

    if(length(aux) >1 && (type(aux[2]) == "t_STR"),
        OUTFILE1 = aux[2];
        print("GFS: ",OUTFILE1 );
    );
    if(length(aux) >2 && (type(aux[3]) == "t_STR"),
        infilestring = aux[3];
        print("GFS: ",infilestring );
    );
    return([OUTFILE1, infilestring]);
}

outputInstanceInfo(fNum, K, lglat_new, reg1, signature_string, prec)={
    \\inputreg = unscaled_determinant(K,lglat_new);
    print("Input determinant ", precision(unscaled_determinant(K,lglat_new),10));
    write(OUTFILE1, "\n--------------------------\n", fNum, " Field pol: ", K.pol,
    ".  Sig: (", K.r1, ",", K.r2, ") -- Precision: ", ceil(REQ_BSGS));
    write(OUTFILE1, strprintf("%-20s %-20s %s\n%-20.9F %-20.9F %d\n", "Log(Disc) ", "Regulator: ", "Disc:", log(abs(K.disc))/log(2), reg1, K.disc), " ",precision(avec,10));
    write(concat("data/table-bsgs-", signature_string), strprintf("%-20.9F %-20.9F %d", log(abs(K.disc))/log(2), reg1, K.disc));
}

get_baby_stock_fit_size(urank, deg_n, detLambda)=
{
    my(coeffA, coeffB, coeffC, sqrt_reg, log_sqrt_reg);
    sqrt_reg = sqrt(abs(detLambda));
    log_sqrt_reg = log(sqrt_reg);

    \\ a*x*log(x) + b*sqrt(log(x)) + c

    coeffA = 0.000368562*(deg_n^2) -0.002715*deg_n - 0.0000533746*(urank^2)+0.00687096;
    coeffB = (262.421 - 96.1079*deg_n + 10.4862*(deg_n^2) + 171.561*urank - 16.4892*deg_n*urank - 30.3957*(urank^2));
    coeffC = (-110.308 *(2^urank)+ 221.575*deg_n -14*(deg_n^2)-779.873*sqrt(urank) + 159.539 *(urank^2));
    \\-15.2457*(deg_n^2)
    \\# absolute value on sqrt(log_sqrt_reg) prevents complex numbers. In this case
    \\# the discriminant is so small the calculation barely matters
    estimate = floor(coeffA*sqrt_reg*log_sqrt_reg+coeffB*abs(sqrt(log_sqrt_reg))+coeffC);
    return(estimate);

}

hybrid_balance_calculator(urank, deg_n, detLambda)=
{
    \\# these coefficients are in terms of sqrtdetLambda
    \\# (see estimate formula in get_baby_stock_fit_size)
    bsgsCoeffA = 0.000368562*(deg_n^2) -0.002715*deg_n - 0.0000533746*(urank^2)+0.00687096;
    bsgsCoeffB = (262.421 - 96.1079*deg_n + 10.4862*(deg_n^2) + 171.561*urank - 16.4892*deg_n*urank - 30.3957*(urank^2));
    bsgsCoeffC = (-110.308 *(2^urank)+ 221.575*deg_n -15.2457*(deg_n^2)-779.873*sqrt(urank) + 159.539 *(urank^2));

    a_1 = 3.18006*10^(-11);
    b_1 = -7.98212*10^(-7);
    c_1 = - 2.8872*10^(-7);
    d_1 = 2.05083*10^(-6);
    a_2 = -5.17116*10^(-15); b_2 = 0.0000156797; c_2 = - 0.000028529;
    d_2 = -0.0000485853; e_2 = -0.0000149705;
    a_3 = 0.0178487;
    b_3 = 0.0000275878;
    c_3 = -0.461283 ;
    d_3 = -0.289706 ;
    e_3 = 1.45846;
    pmaxCoeffA = a_1*deg_n^5 + b_1*sqrt(deg_n) + c_1*(1/(urank^5)) + d_1;
    pmaxCoeffB = a_2*(deg_n^5*urank^10)  + b_2*urank + c_2*log(urank) + d_2*(2^(-deg_n)) + e_2;
    pmaxCoeffC = a_3*(deg_n^3/urank^3) + b_3*(deg_n^3*urank^3) + c_3*(deg_n/urank) + d_3*urank + e_3;

    first_guess = -pmaxCoeffC + bsgsCoeffC + bsgsCoeffA*sqrt(detLambda)+bsgsCoeffB*sqrt(log(detLambda));
    first_guess = first_guess/(bsgsCoeffA+pmaxCoeffA);
    return(first_guess);
}
\\# loop_range is a triple indicating the start, end and increment of the loop
\\# note tha the input files have a specific form, and the vector read in is always called data

\\# param: signature_string is used to specify the read in file and the output file
\\#     ex. "1-1", "1-2", "0-3" etc
\\
\\# param: loop_range is a length 3 vector which specifies which fields to compute on
\\#     they are the arguments for forstep: ex. [3,7,2], will run fields 3,5,7 from data
\\# \param: b_ranges is a matrix of size n x 3, where each row is a triple of forstep arguments
\\# that specifies a range of babystock volumes
\\# each corresponds to the field number's ceiling mod 3 , only because this corresponds to a rough
\\# indicator of the discriminant size.
\\# auxilliary is an additional vector that will hold options.
\\# - auxilliary[1] if present will correspond to the scanball size.
\\# - auxilliary[2] can specify the output file (string)
\\# - auxilliary[3] can specify the input file  (string)
run_bsgs_experiment(signature_string, loop_range, b_ranges, auxilliary)=
{
    GP_ASSERT_EQ(length(loop_range),3);

    suffix = strexpand("(", loop_range[1], ",", loop_range[2], ")");
    [OUTFILE1, infilestring] = generateFileStrings(signature_string, suffix, auxilliary);

    \\\ Note that variable _data_ must be defined in the input file
    read(infilestring);

    GP_ASSERT_TRUE(type(data) == "t_VEC");
    GP_ASSERT_TRUE(loop_range[2] <= length(data));
    if(type(b_ranges) == "t_MAT", GP_ASSERT_TRUE(loop_range[2] <= matsize(b_ranges)[1]););
    if(loop_range[2] > length(data), loop_range[2] = length(data));

    forstep(i=loop_range[1],loop_range[2],loop_range[3],
        timeout = 12*60*60*1000; \\12 hours
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        [K, lglat, reg1, r] = setInstanceVariables(data[i]);
        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        \\# compute expected precision requirements
        [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat, ~reg1);
        default(realprecision, ceil(REQ_BSGS));

        lglat_new = lglat;

        outputInstanceInfo(i, K, lglat_new, reg1, signature_string, REQ_BSGS);

        cpct_units = cpct_from_loglattice(K, lglat_new, eps);

        \\# Fraction of fundamental region to search. Standalone BSGS should use 2 (1/2)
        \\# 1 means you scan the whole region, but this is not needed
        scaleB = 2;

        original_precision = default(realbitprecision);

        if(auxilliary[1] == 0,
            print("No scan radius specified. Auto-selecting");
            default(realbitprecision, 30);
            scanBallRadius = log( 7/sqrt( poldegree(K.pol) ))/2;
            default(realbitprecision, original_precision);
        ,
            scanBallRadius = auxilliary[1];
        );
        if(b_ranges == [],
            \\# behaviour when 3rd arg is the empty vector
            print("Auto-selecting babystock region size");
            sqrt_disc = sqrt(abs(K.disc));
            estimate = floor(0.4811*(sqrt_disc^0.3222));
            init = estimate - 2*floor(estimate/3);
            end = estimate + 1*floor(estimate/3);;
            step = floor((end-init)/10);
            write(OUTFILE1,"babystock-range: ", init, " ", end, " ", step);
        ,
        (length(b_ranges)==3) && (type(b_ranges)!="t_MAT"),    \\elif
            print("Auto-selecting babystock region size based on coeffs");
            sqrt_reg = sqrt(abs(reg1)); \\ this is 'X' in the curve fit fcn
            log_sqrt_reg = log(sqrt_reg); \\ this is Log[x] in the curve fit fcn
            coeff_a = b_ranges[1];
            coeff_b = b_ranges[2];
            coeff_c = b_ranges[3];

            \\ a*x*log(x) + b*sqrt(log(x)) + c
            estimate = max(floor(coeff_a*sqrt_reg*log_sqrt_reg+coeff_b*sqrt(log_sqrt_reg)+coeff_c),10);
            init = estimate - 1*floor(estimate/5);
            end = estimate + 1*floor(estimate/5);
            step = max(floor((end-init)/8),1);
            write(OUTFILE1,"babystock-range: ", estimate, "  ", init, " ", end, " ", step);
        ,
        (length(b_ranges)==3) && (type(b_ranges)=="t_MAT"),
            init = b_ranges[i,1]; end = b_ranges[i,2]; step = b_ranges[i,3];
        ,
            print("Invalid selection of babystock size ranges");breakpoint();
        );

        timeVector =List();         \\ use to track timing changes
        mintime = 0;
        minTimeIndex = 0;
        forstep (j = init, end, step,
            if(length(timeVector)>0,
                trials = length(timeVector);
                if(timeVector[trials][2] > timeout,
                    write(OUTFILE1, "runs exceed timeout. Break to next field");
                    break;
                );
            );

            \\scaling_variable = ((2^r)* log(abs(K.disc))^(1+j/den))/constscale ;
            scaling_variable = j;
            \\write(OUTFILE1, "\nscaling var = log(abs(disc))^(1+",j, "/",den,")*(2*r)/",constscale , "=",precision(scaling_variable,10));

            t9 = getabstime();
            bsgs_output= bsgs(K,cpct_units, scaleB, scaling_variable, bitprecision(scanBallRadius, REQ_BSGS), eps,REQ_BSGS, OUTFILE1, [timeout]);
            t10 = getabstime();

            bsgs_out_lattice = log_lattice_from_compact_set(K,bsgs_output);
            print("result regulator: ", precision(unscaled_determinant(K, bsgs_out_lattice),10));
            print("actual regulator: ", precision(reg1,10));
            write(OUTFILE1, "Overall   time: ",precision(t10-t9,10), "  In mins: " ,precision((t10-t9)/60000.0,10),"  reg ratio: ", precision(unscaled_determinant(K, bsgs_out_lattice)/reg1, 10),"\n");
            overallTime = t10-t9;
            listput(~timeVector, [j,overallTime]);
            if (j == init,
                timeout = min(timeout, max(2*(overallTime), 5*60000)  );
                minTimeIndex = j;
                mintime = overallTime;
            ,\\else
                if (overallTime > mintime,
                    minTimeIndex = j;
                    mintime = overallTime;
                    timeout = min(timeout, max(2*(mintime), 5*60000)  );
                );
            );
        );
    );
}

\\# loop_range is a triple indicating the start, end and increment of the loop
\\# note tha the input files have a specific form, and the vector read in is always called data
run_bsgs_experiment_scaling(signature_string, loop_range, b_ranges, auxilliary)=
{
    print("Function incomplete"); breakpoint();

    GP_ASSERT_EQ(length(loop_range),3);

    suffix = strexpand("(", loop_range[1], ",", loop_range[2], ")");
    [OUTFILE1, infilestring] = generateFileStrings(signature_string, suffix, auxilliary);

    read(infilestring);

    forstep(i=loop_range[1],loop_range[2],loop_range[3],
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [K, lglat, reg1] = setInstanceVariables(data[i]);

        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat, ~reg1);
        default(realprecision, ceil(REQ_BSGS));

        \\
        \\  This takes the log lattice and modifies it so that we get a sublattice.
        \\  The modification depends on the 'latticetype', see function compute_sublattice
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\latticetype = 0;
        \\[lglat_new,modpair1]=compute_sublattice(lglat, OUTFILE1, latticetype);
        \\print(precision(lglat_new,10), "\nMODPAIR", modpair1);

        \\ # can scale lattice columns to create index divisors
        \\ lglat_new = scale_lattice_column(lglat, 1, 15);
    );

}


\\ aux is for variable argument
\\ aux[1] specifies the epsilon value
\\ aux[2] specifies output string
\\ and empty list in single_range will sweep the range up to roughly sqrt(reg)
run_bsgs_experiment_single(signature_string, fieldnum, single_range, auxilliary)=
{
    timeout = 12*60*60*1000;
    suffix = strexpand("(", fieldnum, "_", single_range, ")");
    [OUTFILE1, infilestring] = generateFileStrings(signature_string, suffix, auxilliary);
    read(infilestring);

    forstep(i=fieldnum,fieldnum,1,
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [K, lglat, reg1,r] = setInstanceVariables(data[i]);

        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat, ~reg1);
        default(realbitprecision, ceil(REQ_BSGS));

        print("REQBSGS", default(realbitprecision), "  ", ceil(REQ_BSGS));
        lglat_new = lglat;

        outputInstanceInfo(i, K, lglat_new, reg1, signature_string, REQ_BSGS);

        cpct_units = cpct_from_loglattice(K, lglat_new, eps);
        scaleB = 2;          \\ 1 means you scan the whole region

        \\scanBallRadius =1;
        original_precision = default(realbitprecision);

        \\# note that while scanBallRadius is computed to low precision,
        \\# when input to bsgs, it is wrapped in bitprecision(), which copies
        \\# it as a higher precision (extends with 0). This is needed otherwise
        \\# computations with scanBall radius get truncated
        if(auxilliary[1] == 0,
            print("No scan radius specified. Auto-selecting");
            default(realbitprecision, 20);
            scanBallRadius = log( 7/sqrt( poldegree(K.pol) ))/2;
            default(realbitprecision, original_precision);
        ,
            scanBallRadius = auxilliary[1];
        );
        my(init, end, stepsize);
        if (length(single_range) == 0,
            end = floor(sqrt(reg1));
            init = floor(end/50);
            stepsize = max(init, 1);
        ,
            init = single_range[1];
            end = single_range[2];
            stepsize = single_range[3];
        );
        print(init, " ", end, " ", stepsize);
        forstep (j = init, end, stepsize,
            \\scaling_variable = ((2^r)* log(abs(K.disc))^(1+j/den))/constscale ;
            scaling_variable = j;
            \\write(OUTFILE1, "\nscaling var = log(abs(disc))^(1+",j, "/",den,")*(2*r)/",constscale , "=",precision(scaling_variable,10));

            t9 = getabstime();
            bsgs_output= bsgs(K,cpct_units, scaleB, scaling_variable, bitprecision(scanBallRadius, REQ_BSGS), eps,REQ_BSGS,OUTFILE1, [timeout]);
            t10 = getabstime();

            bsgs_out_lattice = log_lattice_from_compact_set(K,bsgs_output);
            print("result regulator: ", precision(unscaled_determinant(K, bsgs_out_lattice),10));
            print("actual regulator: ", precision(reg1,10));
            write(OUTFILE1, "Overall   time: ",precision(t10-t9,10), "  In mins: " ,precision((t10-t9)/60000.0,10),"  reg ratio: ", precision(unscaled_determinant(K, bsgs_out_lattice)/reg1, 10));
            write(strexpand("bsgs-b-",sigstring, suffix), j, " , ",precision(t10-t9,10));
            if (j == init, timeout = min(timeout, max(2*(t10-t9), 30*60000)); );
        );
    );

}

pmax_log_experiment(signature_string, loop_ranges, auxilliary) =
{
    GP_ASSERT_EQ(length(loop_ranges),3);
    suffix = strexpand("(", loop_ranges[1], ",", loop_ranges[2], ")");
    [OUTFILE1, infilestring] = generateFileStringsGeneral(signature_string, "pmax-log",suffix, auxilliary );

    table_outfile = concat(OUTFILE1, "-table");
    if ((length(auxilliary) >1) && type(auxilliary[1]) == "t_STR", infilestring = auxilliary[1]);
    if ((length(auxilliary) >1) && type(auxilliary[2]) == "t_STR", OUTFILE1 = auxilliary[2]);

    \\\ Note that the input file must define the variable data
    read(infilestring);
    if(loop_ranges[2] > length(data), print("Invalid fields selected: ");loop_ranges[2] = length(data));
    GP_ASSERT_TRUE(loop_ranges[2] <= length(data));
    print(loop_ranges[2], "  ", loop_ranges[3]);
    forstep(i=loop_ranges[1], loop_ranges[2], loop_ranges[3],
        print("field ", i);
        timeout = 24*60*60*1000; \\12 hours

        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;
        \\lglat = process_complex_loglattice(K ,data[i][3]);
        \\reg1 = unscaled_determinant(K, lglat);

        [K, lglat, reg1, r2] = setInstanceVariables(data[i]);
        [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat, ~reg1);

        cpct_units = cpct_from_loglattice(K, lglat, eps);
        fDegree = poldegree(K.pol);
        sumv = lglat[,1];
        for(j=2, length(lglat), sumv+=lglat[,j]);
        X = prec_rigorous(fDegree, log(abs(K.disc)), log(infinity_norm(sumv)) ,abs(reg1));
        default(realprecision, ceil(X));
        outputInstanceInfo(i, K, lglat, reg1, signature_string, X);

        \\
        \\  This is a good spot to modify the log lattice to test sublattice performance.
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        latticetype = 0;
        \\[lglat_new,modpair1]=compute_sublattice(lglat, OUTFILE1, latticetype);
        \\print(precision(lglat_new,10), "\nMODPAIR", modpair1);

        \\lglat_new = lglat; lglat_new[,1] = 3*lglat_new[,1]; print(precision(lglat_new,10), "\nMODPAIR", modpair1);

        lglat_new = lglat; \\modpair1[2] =1;
        unitvector_cpct = cpct_from_loglattice(K, lglat_new, eps);                  \\ computation of compact reps

        memLimit = 1000000000; \\#MEM
        for (k = 1, 1,
        tbefore = getabstime();
            embeddedIntegralBasis = embed_real(K, K[5][1]);
            orthoBasis = gram_schmidt(embeddedIntegralBasis);
            prod1 = gamma(fDegree/2 +1)*memLimit;  \\# Gamma(n+1/2)*MEM
            for(l=1, length(orthoBasis),
                prod1*=norml2(orthoBasis[,l]);
            );
            prod1 = prod1/(8*fDegree*Pi^(fDegree/2));

            newKLimit = floor(prod1^(1/fDegree));
            \\# -1 indicates to use the usual j-value, the last argument says to limit the size of the lowerbound unit search area
            t_x = getabstime();
            indexbound = get_index_bound2(K, lglat_new, eps,-1, newKLimit, OUTFILE1);
            t_y = getabstime(); boundtime = (t_y-t_x)/60000.0;
            print("Indexbound: ", indexbound);
            \\indexbound = k*2000000;
            \\time_estimate = pmax_time_estimate(5,4, indexbound);

            write(OUTFILE1, "Index bound: ", indexbound, ".   bound calc time: ", precision(boundtime,15)  );

            \\print("Algorithm has been commented out to test the index bound");

            logout = log_pohst_pari(K,lglat_new,unitvector_cpct, indexbound, eps);
            tafter = getabstime();
            outreg = unscaled_determinant(K,logout);
            \\write(OUTFILE1,"Output Regulator: ", precision(outreg,10 ), "  quot: ", precision(inputreg/outreg,10), "YN? ",norml2(outreg*quot - inputreg) < eps, ". Ratios: ", (modpair1[2]-inputreg/outreg)< eps);
            write(OUTFILE1, "Output Regulator: ", precision(outreg,10 ), "\n  lpohst time ",precision((tafter-tbefore),10), " Below is in Minutes: ", precision((tafter-tbefore)/60000.0 ,10));
            if (i == loop_ranges[1],
                write(table_outfile, strprintf("%-20s %-20s %-20s %-20s %s\n", "Log(Disc) ", "Regulator: ", "Disc:", "Time(min)", "IndexBound"));
            );
            write(table_outfile, strprintf("%-20.9F %-20.9F %-20d %-20.9F %d", log(abs(K.disc)), reg1, K.disc, precision((tafter-tbefore)/60000.0 ,10), indexbound ));

        );
    );
}

pmax_time_estimate(deg, rank, prime_bound_B)=
{
    my(constant_term, linear_term, b_logb_term, time_estimate, n, r);
    n = deg; r = rank;
    constant_term = 0.0000275878*n^3*r^3 + 1.45846 +(0.0178487*(n^3)/(r^3))-(0.461283*(n/r))-(0.289706*rank);

    linear_term = -0.0000149705-0.0000485853*(2^-n)+0.0000156797*r-5.17116*(10^(-15))*((n^5)*(r^10))-0.000028529*log(r);

    b_logb_term = 2.05083*(10^-6)-7.98212*(10^-7)* sqrt(n)+3.18006*(10^-11)*(n^5)-(2.8872*10^-7)/(r^5);

    time_estimate = b_logb_term*prime_bound_B*log(prime_bound_B)+linear_term*B+constant_term;
    return(time_estimate);
}

guess_function(disc, deg, rank)=
{
    coeff = 1;
    ldisc = log(abs(disc));
    print(precision(ldisc,10));
    print(precision(coeff*exp(ldisc/(2^rank))^(1/4),10));
    print("\n");
    return();
}
skew_lattice(lattice, balanceB, scaling_value)=
{
    my(
        r = length(lattice),
        last_vector = lattice[,length(lattice)],
        large_length = sqrt( norml2(last_vector)),
        target_length = balanceB/scaling_value,
        second_largest_length,
        next_last_vector,
        add_multiple,
        skew_vector
    );
    if(target_length > large_length && (r>1),

        next_last_vector = lglat_new[,length(lglat_new)-1];
        second_largest_length = sqrt( norml2(next_last_vector));
        print("lengths: ", precision(large_length,10), "  ", precision(second_largest_length,10));

        add_multiple = ceil((target_length - large_length)/second_largest_length);
        skew_vector = add_multiple*next_last_vector+last_vector;
        print("target length: current_length " , precision(target_length,10), " ", precision(sqrt(norml2(skew_vector)),10) );

        lattice[,r] = skew_vector;
        GP_ASSERT_NEAR(reg1, unscaled_determinant(K, lglat_new), 0.000001);
        largest_dimension = sqrt(norml2(lglat_new[,r]/balanceB));

    );
    return(lattice);
};

hybrid_experiment(start, end, inputFile, outputFile)=
{
    GP_ASSERT_TRUE((start > 0) && (end > 0));
    GP_ASSERT_TRUE(end - start >= 0);
    read(inputFile);
for(i=start, end,

    \\
    \\ INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    [K, lglat, reg1, r] = setInstanceVariables(data[i]); \\ data must be defined in the inputfile
    n = poldegree(K.pol);
    logdisc = log(abs(K.disc));

    \\# initial precision calculations
    sumv = column_sum(lglat);
    [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat, ~reg1);
    REQ_RIG = prec_rigorous(n, logdisc, log(infinity_norm(sumv)),log(abs(reg1))  );
    default(realprecision, ceil(REQ_RIG));

    \\# If this line uncommented, then the input lattice is just the GRH-assumed unit lattice
    \\# for future reference, this would be the place to modify the lattice into a sublattice to test
    \\# index divisor finding
    lglat_new = lglat;

    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value ", ceil(REQ_BSGS));

    \\write(OUTFILE1,""\nModified lglat ", precision(lglat_new,10));
    \\inputreg = unscaled_determinant(K,lglat_new);
    \\write(OUTFILE1," Input Regulator: ", precision(inputreg,10), "  Original Regulator: ", precision(reg1,10)  );

    p1 = pmax_p1(n,logdisc, log(abs(reg1)) );
    p2 = pmax_p2(n, REQ_RIG, logdisc, log(abs(reg1)));
    g_n = giant_n(n, logdisc, REQ_BSGS,log(abs(reg1)));
    b_n = baby_n( n,logdisc,REQ_BSGS,log(abs(reg1)));
    \\pchoice = p1+p2;
    pchoice = p1;

    \\# identify the index of the basis element with largest norm
    maxnorm_index = 1;
    for(i=1, length(lglat_new),
        if(norml2(lglat_new[,i])> norml2(lglat_new[,maxnorm_index]),
            maxnorm_index = i;
        );
    );
    print("max vector norm index: ", maxnorm_index, " of ", length(lglat_new), ". Length: ",precision( sqrt(norml2(lglat_new[,maxnorm_index])),10) );
    if (maxnorm_index != r, print("warning, the last element of the integral basis is not the longest."););

    \\# old formula for B
    \\balanceB = ((abs(reg1))^(1/3)) * (g_n*b_n)^(1/3); balanceB /= (pchoice^(2/3));

    oldB = abs(log(reg1))*2^poldegree(K.pol)*reg1^(1/3);
    balanceB = hybrid_balance_calculator(r, n, reg1);
    balanceB = min(reg1, balanceB);

    print("Degree: ", n, " UnitRank: ",r, " Original B ", floor(oldB), " Fitted B: ", precision(hybrid_balance_calculator(r, n, reg1),10));

    \\balanceB = min(balanceB, sqrt( norml2(lglat_new[,length(lglat_new)])  ) ); balanceB*=2;

    \\balanceB = 1900; print("warning, hard coded value for B");
    largest_dimension = sqrt(norml2(lglat_new[,maxnorm_index]/balanceB));

    last_vector = lglat_new[,length(lglat_new)];
    large_length = sqrt( norml2(last_vector));
    target_length = balanceB/3.0;


    lglat_new = skew_lattice(lglat_new, balanceB, 3.0);
    largest_dimension = sqrt(norml2(lglat_new[,r]/balanceB));

    [REQ_BSGS, REQ_COMPARE, eps] = compute_precision(~K, ~lglat_new, ~reg1);
    REQ_RIG = prec_rigorous(n, logdisc, log(infinity_norm(sumv)),log(abs(reg1))  );

    write(OUTFILE1, "hybrid B ", precision(balanceB,10), "  Old B:  ", precision(oldB, 10));


    print("Running Pohst Algorithm");                                           \\ lglat_new is the input lattice, pohst_out_lattice is the result after ruling out index divisors up to pohstB

    unitvector_cpct = cpct_from_loglattice(K, lglat_new, eps);                  \\ computation of compact reps

    tbefore = getabstime();
    pohst_out_lattice = log_pohst_pari(K, lglat_new, unitvector_cpct, balanceB, eps);
    \\pohst_out_lattice = lglat_new; print("warning, pohst step commented out ");breakpoint();
    stage1_units = cpct_from_loglattice(K, pohst_out_lattice,eps);
    tafter = getabstime();

    lptime = tafter-tbefore;
    \\ Just checking the regulator of the output from the p-maximization
    \\write(OUTFILE1,"Pohst Output Regulator: ", precision(outreg,10 ), ". Ratios: ", (modpair1[2]-inputreg/outreg)< eps);
    write(OUTFILE1, "pmax time ",precision(lptime,10), " In minutes: ", precision(lptime/60000.0,15) );

    print("Running BSGS Algorithm");
    default(realprecision, ceil(REQ_BSGS));
    detLambda = unscaled_determinant(K, pohst_out_lattice);
    \\print("REQ_BSGS ",floor(REQ_BSGS) );
    t9 = getabstime();

    \\\ reduced precision to compute scanball radius
    original_precision = default(realbitprecision);
    timeout = 12*60*60*1000;
    default(realbitprecision, 30);
    scanBallRadius = log( 7/sqrt( poldegree(K.pol) ))/2;
    scanBallRadius = min(largest_dimension, log( 7/sqrt( poldegree(K.pol) ))/2);
    default(realbitprecision, original_precision);

    g_n = giant_n( n, logdisc, REQ_BSGS, real(log(detLambda)) );
    b_n =  baby_n( n, logdisc, REQ_BSGS, real(log(detLambda)) );
    \\ This is the basic calculation for the babystock region
    scaling_variable = sqrt(  (abs(reg1)/balanceB)*g_n/b_n  )/2;

    detLambda = unscaled_determinant(K, pohst_out_lattice);
    fitted_scale_variable = get_baby_stock_fit_size(K.r1+K.r2 -1, poldegree(K.pol), detLambda/balanceB);

    print("default Bstock Size: ", precision(scaling_variable,10), ". Fitted Size: ", precision(fitted_scale_variable,10));
    write(OUTFILE1, "babystock size via fit: ", precision(fitted_scale_variable, 10), "  Default value was: ", precision(scaling_variable,10));
    if(fitted_scale_variable > 0,
        scaling_variable = fitted_scale_variable;
        , \\else
        print("fit value produced a bad babystock size. Using default value");
    );

    bsgs_output= bsgs(K,stage1_units, balanceB, scaling_variable, bitprecision(scanBallRadius, REQ_BSGS), eps,REQ_BSGS, OUTFILE1, [timeout]);
    t10 = getabstime();
    bsgstime = t10-t9;
    bsgs_out_lattice = log_lattice_from_compact_set(K, bsgs_output);
    outreg = unscaled_determinant(K,bsgs_out_lattice);

    write(OUTFILE1, "bsgs time ",precision(bsgstime,10), " In minutes: ", precision(bsgstime/60000.0,15) );
    write(OUTFILE1,"Overall time: ", precision(bsgstime+lptime , 10) , " In minutes: ", precision((bsgstime+lptime)/60000.0,15) );

);
}


specialField11 = [x^3 - 14105*x^2 + 190241206*x - 1434835874037, -22729311273742752391460247, \
[472127480231.19058564549342808743779377712964461705822041039616921390834216514742089426164608700285224539653227784071320317983861427662412196184167337903859099635926348062696842208555837034430238487927123938442366161078114770082570477123549304314589048301328471704975818132364412646228885384176662444724730776560094482623627653752559248545297505006989109270575275833095985203116970896510434270727118438348758227920519567635721951088609069754806222097445919409332199539366007538296928313392888844385975 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037715054497824558763660238983*I; -472127480231.19058564549342808743779377712964461705822041039616921390834216514742089426164608700285224539653227784071320317983861427662412196184167337903859099635926348062696842208555837034430238487927123938442366161078114770082570477123549304314589048301328471704975818132364412646228885384176662444724730776560094482623627653752559248545297505006989109270575275833095985203116970896510434270727118438348758227920519567635721951088609069754806222097445919409332199539366007538296928313392888844385975 + 11.526005437528880337782384385123881089305298628901839535194975994927139034023654804794481258060837178537264612315098878730595190738165646846292171903515321479572585447781821209483527129178419194846600993280189034503894833762207375482528939811527790171178887104147603199155128340767552782934004910134236850448872574047836334223126587696428309062265188586148249612363345771382750511289965920273871747740817197899585290137242538418885207949771433323603174743205482689468297695044620521350987572083827597*I]];
{

specialField30 = [x^3 - 75498*x^2 + 1680279841*x - 11308481410487, 20588425792664626980132449, \
[371148.25103310256966453964730327803928545022729495216302508628256872297237054788826750776010068758174518332163918389780075653472927630658880199153790572579499382394088516041352584450520705255894997851522945434126458369063039810862003079819746758150660317138165348546466707783135896377638484615850216474219250126419124400044753632092648289195992376434626566424506030283155296912092231674605031939516977340036852053253137610661428292603623427259981203254434514572899444978878676659495592867357793115136 + 1.1253584242118130320 E-694*I, 507707.91102222347698083838234323742028154550500841537861200856150736619812197190172131611521661609275373190868698653129734217516817366138142931800877010442640044244037722181157916314244096578793391875439852004223386900445298957546938578896524629268664510857718560839135553301248338386914250136378297262099699463212035088437511341921041856834788539488484036995335123934567184829285429206816115086439505475996273662233774015008215903710985686023318599530692603238885868422343167303817762023164535725665 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491*I; -387290.38195408775202175803377179786936418551463946781908833212430459669789095559821563927726326509671636620038078654816653761739696625149724845407997445145536008106355250782263243557007160749524176354854573415461244799104107844757583893341912102826249356851138815736611930946825860972523742237915634501102487499440221736466151476585951660524877106242321280183408754949035486348873367652729587196914982980190403030082867633125888196525310247965974213243472079330219233333159459028488150181732978779962 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491*I, 453272.87391710002175566353703854009408358985520738640982244571359076944409648085384198003798103429716194550996518009540235136147498080199772861331771299778765697369882614194005561890540355850694911710354529407088813261872541074831463786660769780248536403244254231660887070338111054975278062573298255672582769292792322884045836155223485732124109237250868236831806087850149436729824264124652810522282706041460398269529977798669063784768577201869642776229992966363213474156349420136789021653428932779297 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491*I; 16142.130920985182357218386468519830078735287344515656063245841735873725520407709948131517162577514971182878741602650365781082667689944908446462542068725660366257122667347409106591064864554936291785033316279813347864300410680338955808135221653446755890397129734671901452231636899645948852576220654180268832373730210973364213978444933033713288847298076947137589027246658801894367811359781245552573980056401535509768297300224644599039216868207059930099890375647573197883542807823689925573143751856648261 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491*I, -960980.78493932349873650191938177751436513536021580178843445427509813564221845275556329615319765038991567741865216662669969353664315446337915793132648310221405741613920336375163478204784452429488303585794381411312200162317840032378402365557294409517200914101972792500022623639359393362192312709676552934682468756004357972483347497144527588958897776739352273827141211784716621559109693331468925608722211517456671931763751813677279688479562887892961375760685569602099342578692587440606783676593468504962 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491*I]];

specialField40 = [x^4 - 716*x^3 + 164632*x^2 - 12722131*x + 212364731, 73207246861319109539698517, \
[6061.7777125926362445835851632875416221836779355940788152677969404825737213904323773851610570690040080007169619760344243486728128709246534926253623069571075500674188394196958571061213750786030299049581755803140908828040473046315101940542492684450971791131812450532037084031873389855170145765380630981419974628989955671746665102209201428716287674004632746508899628568955111778993222614552620265267980186635949186007144443699780237678155385204791 + 1.08583322703270845029911236857200407395 E-405*I, -1481.98398295182931832213286526231571870708628986524490535857825291199235255542286845683774492101925125042541246431380293315543660349853573396560029937871055283129840344971796127571004122255175738943615806930287535318347628157839540417396201782077681449618762128334738842621139142463476736040320064846807465272361785913747614774715887805571577678378163157620151687468685824050039866669197242802424342999064480097342183125885973619381782626975521 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661114513780091951036380490021936107158751*I, 5615.0581550445460424672157125478128682189830114540621313045488887518983952875308952729337586528724642328966186856615616755711766023230282815078496159292531656597186263233584569236899252853160061235456599943358314995856696022590505458768228095048268510443209022929774505325855184563085738552397783474554679030999202015036355455051475736996180477927222971366693659636047445559411885668048869525273284295038890592645008764527459801430706251641224 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661130107117832683801057859335688990611821*I; 1041.646497325773175749487370064423931631205042159724854527699606239576262241970046540824579994578216387700974731773560062938969058742040279754194485972440118107757217767721439768276758001956343540352780034239037470739650547393827151302499956752500902588568763588895106324260410051208950088132811079025919232879778910349997859971113423520162395542515157833152012095718062923277236999212497525607568610742725162956198394765410 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330580287415677794164745190650006844106702*I, 6868.386210004924028941409762633679813570046045605573460304429354627870361630700498007024918136525990117502615857468582287874794893810504502651893992172965248301729267466586261673278039460069452646818883711827937182380248778390198484192761982286760601741063053611770506007682178115092881205726830801084374434195946474379207215713246291995450614962514100715438013358136631040141693510604541521417541015793811803318503781989868 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330551126762483609982149115789108955944749*I, -6090.332365537015515990541570202574565027593193407062265205498289208553845960196748525001727008851608838921794852712579134691036786652365904279570589729568739803535641921697362502731377954726484766146361641179199942509349492580618191107778905078934552727353995594076505987657894862805003386063576368166596757386761735788713967416631585597346104132904449306954162675436685301007881890379565348360635979766157052896422834163891 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330562017578694663431799952629417086514249*I; -6282.6184937227751316622013702951655255702699916376863866167269885659215082594514790266467601956198636695359374416128776304209659438806677825706507009158287855250262310293405627027924393901539242799279546788047859095326059654247059194033443364559109850754415813302147461482298440143324293693632154671346422347531650664040436884082912088541179800912823145100176795739660797732635496609270899129223495470017163391524952899817134805013879481776435 + 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092230474532810252977406661881550170096703*I, 790.34933969676748330124724569850185025977669021129116727007642253723850432703857143258125905409658536744104701933472304757479775928589421716532174165998627680092772612615314235910770401562196995572993849505249887412091058841547343132499662798072692756955781317889261746038069367945517785947989402685529689074220736196661346293493388611662193939942180755843524481857101343012288292905860776733644974633328241141734574173573914908513092195686481 + 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413596429617302656461329418768921910116446345071881625696223490056820540387704221111928924589790986076392885762195133186689225695129646757356633054240381829129713384692069722090865329642678721452049828254744917401321263117634976304184192565850818343072873578518072002266106109764093304276829390388302321886611454073151918390618424603221604724657577954658556263442978*I, 5547.5058933133894439202397651430798724437068150811626576122859110738101161192586748962346349923099949015007427845960221532839634076750309694198641533687160716822269789649467291182254596658287884330455990194488222055827582098975693705858672527702023169171992854319423247749577269658320447804517388015842210147114801076586167122989778811212643114814172074453557654879500864811306373678305607876175675383664948778629831430144582307013604224044994 + 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413596429617302656461329418768921910116446345071881625696223490056820540387704221111928924589790986076392885762195133186689225695129646757356633054240381829129713384692069722090865329642678721452049828254744917401321263117634976304184192565850818343072873578518072002266106109764093304276829390388302321886611454073151918390618430973224084220195734715045040375988561*I; -820.8057161956342886708711630568000282446129861161172831787695581562284753729509448993388768679623607188819992661951067811908159857860259898089060920137188826501498261580767341716056936904054491653830009357483424440110918866006314259534048887416870966263084273118840685792179050223935352953076587100332744610256094111206206817837423575376731828516961179740242953786474943279130095997406696392120170824046037424044175491536744 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530922102631442353878030*I, -6176.751566749862193920524143069865945122736445951619722215927524253116513402316200982768432269603324234518250412489502402294156049597862985851615434454240972271358590143021442756675702253139665213112664137577560703317683085227276511343796592446710714814433245507315735041851480369913291704803524179471596672214535977208344530901021300056356777578154276697671741302020786229764177772971176860729747332136449413762427692466748 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843021390217995818798*I, -5072.231682820919970396913907488318175635096633128162523711336510617154665446592821644166666636330850295475566617545004694164103223345693346648143179568400497538409963366607823539184006996418309790444897372605453762659078319576001725354911157196094615234166192130843269319885350559335615249627940780873092160424638573373538290387493869223536255141235055275070968776118145736063944044255882391784259988104226884231061185303313 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843398379726471003131*I]];

specialField50 = [x^5 - 151*x^4 + 8417*x^3 - 215351*x^2 + 2487446*x - 10000191, 71778511193858521448191097, \
[-198.19322487602258281527102909827398107350284621631742572143956007314647468297311202670951984439916387146216016569580061281522740486822648296163869327298071490431583060958066301089952399268737604360719931890504792266321113531257438934103198240211691309167061396102958403293679854324455163912936507085560742911379693410475520494518882797077214198260560727949568651276873774203579955048809445969523419265376187991408861127113163190956055849044177505653236845674982134900303701626143458015268411966313990 - 8.3027583391042336762001816252565291060 E-477*I, 300.00109521623450079820139797920496858610136180816497614278601807589800958066142957411695325045128890711065181739690984861303871183941226418773083413328473413819815174089730655061552822978902305388428328202501196344340105898501625074225446844169919875627297228572470993307895944564197883663918540335050224645610585751672960827778516733140253126778539417738153182985701239731853576205246631115552969279730084620997581216971232337019830032335464761857216984048882530376766043098629235387688434598862228 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857651790287365944932061784*I, -937.70125739136196330768435326403533420708127072166961928596636585322397263941858869439056017856840777342928264526543749350068103597113841996883461063111831181814694778192909993295451356527653673182540443536155280781376909405222592719895613618698048686163314908724591748782251848072521156360725003036749666363264018403620791516445973454431904056521664192273955220504010487969342108475963143527342057675342386199401397786015298134514629173593663081918785514699931778998791466071762875055056117586184388 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857215895474562972674739687*I, 588.15821527381287285079405399835657569929409433113226290908885489099893172651524150580021303720156910841001217702147670900045790116923073310251644233810977096647424954774590395831118986044668204761742669514762895958153921837643720489999866259694554668375794142477029331615479777556055968650246372894719429880001444541035504362024492666058839532110198708699799150524469014587943926324175402046000142014753915963957052457966824318931902621665837985865137486886490609557727867896473573117880838438747372 - 1.66055166782085522512152837011240954122 E-476*I; 209.48944033020894483271595530176736455101529233590847211445357705933417070956964719928975168499268774095187376680637351577624923313543210860586933354862311506001481743134685534445941137855014233164262722954609161373198935486875219487263024091942393114758749042989818620622906920786994413046375719760950200405861258922089721776825856402011249013327106066887637691493218248362088429248025232798320492808644631608212623897556815371372442633538880121615579094240959554170781328340161640915842656923860621 - 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518858660990907621302224422634*I, -786.13967505520855599372547568784136698872555012375900309529489232296119631054434185341884386732344743838775310277641757358007839850251270301255417951088310037660130262932177910849243092216911326708273948677412737382247136591713717423517550899207728391928212975094279732306497963896038105355983400818649867857203019568104853496229332283903858232930044821505802991399985364973178118758817185616265751870630961124100007647783945380454510978695489017038706131820087549866525054308611973024727204117904794 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518859227861905292256915840207*I, 246.37406681039903843902746787987995425834109719470805621599207772259936070187483021104184250491499417554399262977835617107903957122383763011401693544178224728194165886193289782942026061467777633811174630688430259795344271895141679367505295223709008696581675289034620717179054920946679759612307765652851141324862337642770006878550175754773440875193659523097149456078314747930995994530757531818501240511827541348659985303083054112242993098954200362501315677316225528910817730188359080573820399991461128 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037710802965342026556821669642*I, -392.48054363708696162711038728395473115929253025517340703028803075697957504906266403206697255336713083131383342000098020176456804188739595131553933450126532208136772503618556886937476227301569214681714524960448380976334282846201743234743902737749210193786082135841660287366078630492538582492398708503214703841395109663300445860562326298617445954628415139495608706498176691233376120915141057450102819889124208942774677707093965346343230644776512705456687047399968943536806877293469255360532763342440704 - 2.26748399068385238161018427610330178742 E-475*I; -487.68906267267981082433818858972465061900034128103459347957132843550573946051155665407191533125054605384128157076588378540657250384935203418181786922988749813479303035932393166710495000023711281791061183085063764944581941397311427885402357447223832913692617243136834868897529790462858603592581342944099031637176357456023162106940017983108105223927353252804424272335111344501309052727877861435874031960587789187360040095884060306517921897946035742717752468358069668823731930923583250737391971818554065 - 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037715005355412405397506755189*I, 93.740573446922549534754407179936844007918160265600276737572197520059823294088997479980545167866054580575340084429570679358138098341654965106193803287677822663376723481052202280224766732651246315291164240992415598271553691428086928305032305447673596047669719362006886958264313176787325893698213121244187088307129647005687223561466639759805759810312198746734490585015543081031830757934347565989543614052788389809503536744490533516818338371793842921452519323689962332152425077184830440454077678272445812 + 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857453535294049332838543765*I, 303.75580048159592027654401169789257683402348789686708129815354989329527318163687618364910941979421666964282446735670166749698678538252809989504342093005294954741204759568634990032081031852310266722536133180089876272275260741071261421827582666373478190471086060666463106123249304422765794408731116986713256301703478275412102768183228830127898057016623823814967595579430813787863376103262983410907039788616799163105719190308926686566069472098107106061459991287766805533784540086215970280470175758777553 + 1.84284045575123307886520025537874057418 E-476*I, -39.003639829475272815901278804796805371621036813352893399759916909973483751721190154363670560769322855232396151357507399594411035716290968554033780537738978331864163762120783987494942789818547719745912457960706008680659685815469548611440171045335263898891860387819686193698503530704416126732535733432458513489026310300160618530164951966089040137803249923727636709271478156008655238980645883309771985801738724281665233463741057006860222491890146494817954957595085061600166414839145927040760373880910696 - 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857428964087972649020660615*I; 125.71410900772053960452573522038383910240082204267965688162167458955533568709810013666428142083096952836623063180246911287485997574426075537912807915717152964896962675542831830959447805087509552893114025508363273092673617883914743041619831611373959656709497284851986926190647566417649389596888481138190923520697336250544118584955175833422440913700348522470937082700830971546915339187886839945418310349925763715389839479105657497906942042772629830862301204911428729085018726585413178865839704809980649 - 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518858097628464101564276283537*I, 292.06773606196060167034954069574373806114806297498293065614125437588399867259407095015349344176493265305684998523101909667581190673389970578289195231973100155211966641177802463939393155186228028807086244871771995145080631771010928971342688408706573084016881116657438315061529660927430875025335136364465929499342358207840425328304141879442672424184414374967205178353955095091401462568729665416810375309332321972936283280135715892545029762705828813069067662294545374213908484543106682092246434312373384 + 8.5556932773327089933669827670030323241 E-476*I, 16.629737468374052436489464723852858001121627336065170338023195436577169406650577789997416774787669037685314494382891429367869306058682346649862886404004890556353797646279955903528676230221062258958831143840588031067995310317797082822351522050663697679797261124903796357112066523792025593260769146425680348630469265396303149935491774007638216822521533258410927298477166806093130205386043518728365619049387081176249404521811729061781488272837759944322357760571299915425694296339819088388536835300097743 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037712915574505225576790581289*I, -639.83389648884306945166757412852752867297703650850003572208692916952270018219578039932255653780703771513009607183923898532145052522533305142264437344488699463076413864328135002274102658650209179544156706108472470695905997173228797098183313200809706277498784548718832911287207580582311748294999432261826524034503207435012525541428328164370487128540302466537725501568017623997500980736464122041538175343351890700677360551588615174556046581946719494806068591467213206041498841166957333845538120925972437 - 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518858668008015923752554802447*I; 350.67873821077290920236752716584742803908707311876389020493563685976270774681692134482740206982605265598533733785284176957069069983788565315845914979707356833012441678212942102395058456349925100094404366512596122745030501557778904290622699984119171451391432311397987725377655157582669964862253649130518650621997455693864842239677868544751629495160459391395418149417935898795885239340775234661658648067393581855166437846334750628194593070678703295893109014880663520468235577624151888970978022051026786 - 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037714662641759439269120571306*I, 100.33027033009090399042012983295581633355796507501081955879542235111936476319984384916785200724117129764491121571891794893308968158754576793573758977018954202290676099559424563825820440786656360983642951503897986065671029779392470547446185101563875827517062693663681728110641040725676757296908411994715004881537110908022744984000009695340356700935871154126995571558774722046740004191406132484948045876289715549215789476227943799207817346474811149967169553107663412060608018948393011499384567379424602 + 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037714466713726879523753688686*I, 370.94165263099295215562340896240994511359505829402931143379754280075216934925630450970219147907152789055715105374748822555678537330609034330991136785527822443243944367802989629968476640185459546752946565283576341606957845737229943648327583523549192031130827446533128289768740970323873043013609205754617233873651275945808366876163391468766743442059227519520745438998548245641169717303338276425097215469959337570010752840442144429527417775257579618923774070038809453011619766163205915361911858305935931 + 1.46946024419810637630914057839555201523 E-475*I, 483.15986468159243104388518621892248950459650924589407324304602194547682725646439307995298661474192229326631346617624987767997170165978923818970104614578152407752177789384179892129954178888964961438719807350228556582152326763333774704071366783397888192798258580865432486407656786589235974810405341213567649344799503587293528892982656993537997564838843889706298728468873116243798699225494365776618051797896056107661509147089861902653396854246408863879413647740200046180594492047867608792266083217756838 - 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037714270785694319770574012546*I]];
