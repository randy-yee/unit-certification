read("src/BabyStepGiantStep.py")
read("src/CompactRepresentation.py");
read("src/bounds.gp");
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);
DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;

\\\\ BSGS functions
setInstanceVariables(readData)={
  K = nfinit(readData[1]);
  lglat = process_complex_loglattice(K ,readData[3]);
  reg1 = unscaled_determinant(K, lglat);
  return([K, lglat, reg1]);
}

compute_precision(~K, ~lglat, ~reg1)={
    my(X1,X2,sumv = lglat[,1]);
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K.pol), log(abs(K.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K.pol), log(abs(K.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));

    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    REQ_COMPARE = ceil((poldegree(K.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K.pol)^2 +5);

    return([REQ_BSGS, REQ_COMPARE]);
}

scale_lattice_column(loglat, col, factor)={
    print("Scaling lattice column ", col, " by ", factor);
    copy_lattice = loglat;
    copy_lattice[,col] = factor*copy_lattice[,col];
    return(copy_lattice);
}

\\# loop_range is a triple indicating the start, end and increment of the loop
\\# note tha the input files have a specific form, and the vector read in is always called data

\\# \param: signature_string is used to specify the read in file and the output file
\\# ex. "1-1", "1-2", "0-3" etc
\\# \param: loop_range is a length 3 vector which specifies which fields to compute on
\\# they are the arguments for forstep: ex. [3,7,2], will run fields 3,5,7 from data
\\# \param: b_ranges is a matrix of size n x 3, where each row is a triple of forstep arguments
\\# that specifies a range of babystock volumes
\\# each corresponds to the field number's ceiling mod 3 , only because this corresponds to a rough
\\# indicator of the discriminant size.
\\# auxilliary is an additional vector that will hold options.
\\ # currently, the only 1st position of aux, if present will correspond to the
\\# scanball size.
run_bsgs_experiment(signature_string, loop_range, b_ranges, auxilliary)=
{
    GP_ASSERT_EQ(length(loop_range),3);
    timeout = 172800000; \\ 2 days

    outfilestring = strexpand("data/data-bsgs-",signature_string,"(", loop_range[1], ",", loop_range[2], ")" );
    print("Output directed to file ", outfilestring);
    \\infilestring = concat(concat("input/test-poly-",signature_string),".gp");
    infilestring = concat(concat("input/polynomial-",signature_string),".gp");
    OUTFILE1 = outfilestring;

    if(length(aux) >1 && (type(aux[2]) == "t_STR"),
        OUTFILE1 = aux[2];
    );
    read(infilestring);

    forstep(i=loop_range[1],loop_range[2],loop_range[3],
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [K, lglat, reg1] = setInstanceVariables(data[i]);

        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        [REQ_BSGS, REQ_COMPARE] = compute_precision(~K, ~lglat, ~reg1);
        default(realprecision, ceil(REQ_BSGS));
        print("REQBSGS", "  ", ceil(REQ_BSGS));
        eps = 2^(-REQ_COMPARE);

        lglat_new = lglat;

        \\inputreg = unscaled_determinant(K,lglat_new);
        print("Input determinant ", precision(unscaled_determinant(K,lglat_new),10));
        write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol,
        ".  Sig: (", K.r1, ",", K.r2, ") -- Precision: ", ceil(REQ_BSGS));
        write(OUTFILE1, strprintf("%-20s %-20s %s\n%-20.9F %-20.9F %d\n", "Log(Disc) ", "Regulator: ", "Disc:", log(abs(K.disc)), reg1, K.disc));
        write(concat("data/table-bsgs-", signature_string), strprintf("%-20s %-20s %s\n%-20.9F %-20.9F %d\n", "Log(Disc) ", "Regulator: ", "Disc:", log(abs(K.disc)), reg1, K.disc));

        cpct_units = cpct_from_loglattice(K, lglat_new, eps);

        scaleB = 2;          \\ 1 means you scan the whole region
        r = K.r1+K.r2 -1;

        original_precision = default(realbitprecision);

        if(auxilliary[1] == 0,
            print("No scan radius specified. Auto-selecting");
            default(realbitprecision, 30);
            scanBallRadius = log( 7/sqrt( poldegree(K.pol) ))/2;
            default(realbitprecision, original_precision);
        ,
            scanBallRadius = auxilliary[1];
        );

        if (ceil(i/3)== 1, init = b_ranges[1,1]; end = b_ranges[1,2]; step = b_ranges[1,3];);
        if (ceil(i/3)== 2, init = b_ranges[2,1]; end = b_ranges[2,2]; step = b_ranges[2,3];);
        if (ceil(i/3)== 3, init = b_ranges[3,1]; end = b_ranges[3,2]; step = b_ranges[3,3];);
        if (ceil(i/3)== 4, init = b_ranges[4,1]; end = b_ranges[4,2]; step = b_ranges[4,3];);
        if (ceil(i/3)== 5, init = b_ranges[5,1]; end = b_ranges[5,2]; step = b_ranges[5,3];);
        if (ceil(i/3)== 6, init = b_ranges[6,1]; end = b_ranges[6,2]; step = b_ranges[6,3];);

        timeVector =List();         \\ use to track timing changes
        forstep (j = init, end, step,
            if(length(timeVector)>0,
                trials = length(timeVector);
                if(timeVector[trials][2] > timeout,
                    write(OUTFILE1, "runs exceed 2 days. Ending job");
                    quit;
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
            write(OUTFILE1, "Overall   time: ",precision(t10-t9,10), "  In mins: " ,precision((t10-t9)/60000.0,10),"\n" );

            listput(~timeVector, [j,t10-t9]);
        );
    );
}

\\# loop_range is a triple indicating the start, end and increment of the loop
\\# note tha the input files have a specific form, and the vector read in is always called data
run_bsgs_experiment_scaling(signature_string, loop_range)=
{
    print("Function incomplete"); breakpoint();
    GP_ASSERT_EQ(length(loop_range),3);

    outfilestring = strexpand("data/data-bsgs-",signature_string,"(", loop_range[1], ",", loop_range[2], ")" );
    print("Output directed to file ", outfilestring);
    \\infilestring = concat(concat("input/test-poly-",signature_string),".gp");
    infilestring = concat(concat("input/polynomial-",signature_string),".gp");
    OUTFILE1 = outfilestring;
    read(infilestring);

    forstep(i=loop_range[1],loop_range[2],loop_range[3],
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [K, lglat, reg1] = setInstanceVariables(data[i]);

        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        [REQ_BSGS, REQ_COMPARE] = compute_precision(~K, ~lglat, ~reg1);
        default(realprecision, ceil(REQ_BSGS));
        eps = 2^(-REQ_COMPARE);
        print("REQBSGS", "  ", ceil(REQ_BSGS));

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
run_bsgs_experiment_single(signature_string, fieldnum, single_range, auxilliary)=
{

    outfilestring = strexpand("data/data-bsgs-",signature_string,"(", fieldnum, ",", ")" );
    print("Output directed to file ", outfilestring);
    \\infilestring = concat(concat("input/test-poly-",signature_string),".gp");
    infilestring = concat(concat("input/polynomial-",signature_string),".gp");
    OUTFILE1 = outfilestring;

    if(length(aux) >1 && (type(aux[2]) == "t_STR"),
        OUTFILE1 = aux[2];
    );
    read(infilestring);

    forstep(i=fieldnum,fieldnum,1,
        \\# INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
        \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        [K, lglat, reg1] = setInstanceVariables(data[i]);

        \\# in case bnfinit is needed
        \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;

        [REQ_BSGS, REQ_COMPARE] = compute_precision(~K, ~lglat, ~reg1);
        default(realbitprecision, ceil(REQ_BSGS));
        eps = 2^(-REQ_COMPARE);
        print("REQBSGS", default(realbitprecision), "  ", ceil(REQ_BSGS));
        lglat_new = lglat;

        \\inputreg = unscaled_determinant(K,lglat_new);
        print("input determinant ", precision(unscaled_determinant(K,lglat_new),10));
        write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol,
        ".  Sig: (", K.r1, ",", K.r2, ") -- Precision: ", ceil(REQ_BSGS));
        write(OUTFILE1, strprintf("%-20s %-20s %s\n%-20.9F %-20.9F %d\n", "Log(Disc) ", "Regulator: ", "Disc:", log(abs(K.disc)), reg1, K.disc));
        write(concat("data/table-bsgs-", signature_string), strprintf("%-20s %-20s %s\n%-20.9F %-20.9F %d\n", "Log(Disc) ", "Regulator: ", "Disc:", log(abs(K.disc)), reg1, K.disc));
        cpct_units = cpct_from_loglattice(K, lglat_new, eps);

        scaleB = 2;          \\ 1 means you scan the whole region
        r = K.r1+K.r2 -1;

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

        forstep (j = single_range[1], single_range[2], single_range[3],
            \\scaling_variable = ((2^r)* log(abs(K.disc))^(1+j/den))/constscale ;
            scaling_variable = j;
            \\write(OUTFILE1, "\nscaling var = log(abs(disc))^(1+",j, "/",den,")*(2*r)/",constscale , "=",precision(scaling_variable,10));

            t9 = getabstime();
            bsgs_output= bsgs(K,cpct_units, scaleB, scaling_variable, bitprecision(scanBallRadius, REQ_BSGS), eps,REQ_BSGS,OUTFILE1, [172800000]);
            t10 = getabstime();

            bsgs_out_lattice = log_lattice_from_compact_set(K,bsgs_output);
            print("result regulator: ", precision(unscaled_determinant(K, bsgs_out_lattice),10));
            print("actual regulator: ", precision(reg1,10));
            write(OUTFILE1, "Overall   time: ",precision(t10-t9,10), "  In mins: " ,precision((t10-t9)/60000.0,10),"\n" );
        );
    );

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
\\guess_function(9699284,4,3);
\\guess_function(5901732332357, 4,3);

\\guess_function(2113524069716037, 4, 3);
