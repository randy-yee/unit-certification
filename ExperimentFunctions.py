read("src/BabyStepGiantStep.py")
read("src/CompactRepresentation.py");
read("src/bounds.gp");


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
run_bsgs_experiment(signature_string, loop_range)=
{
    GP_ASSERT_EQ(length(loop_range),3);

    outfilestring = strexpand("data/data-bsgs-",signature_string,"(", loop_range[1], ",", loop_range[2], ")" );
    print("Output directed to file ", outfilestring)
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
