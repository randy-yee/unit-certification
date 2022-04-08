read("src/PmaxLog.py")
read("src/PmaxNormal.py")
read("src/BabyStepGiantStep.py")
read("src/bounds.gp")
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
\p500
default(parisizemax, 8000000000);
\\ Global variables
eps = 10^(-100);      \\ error tolerance
sqrt2 = sqrt(2);
DEBUG_CPCT = 0;
DEBUG_REDDIV = 0;


OUTFILE1 = "data/tmp-experiment-hybrid.txt";

\\read("input/test-poly-1-1.gp"); data = data1_1;
read("input/test-poly-3-0.gp"); data = data3_0;
\\read("input/test-poly-4-0.gp"); data = data4_0;
\\read("input/test-poly-2-1.gp"); data = data2_1;
\\read("input/test-poly-0-2.gp"); data = data0_2;

\\read("input/test-poly-1-2.gp"); data = data1_2;
\\read("input/test-poly-3-1.gp"); data = data3_1;
\\read("input/test-poly-5-0.gp"); data = data5_0;
\\read("input/test-poly-0-3.gp"); data = data0_3;
\\read("input/test-poly-2-2.gp"); data = data2_2;

\\read("input/test-poly-4-1.gp"); data = data4_1;
\\read("input/test-poly-0-4.gp"); data = data0_4;
\\read("input/test-poly-1-3.gp"); data = data1_3;


{

for(i=1, 10,

    \\
    \\ INSTANTIATES THE FIELD AND THE LOGLATTICE OF UNITS AND CPCT REPS
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    K = nfinit(data[i][1]);
    \\K1 = bnfinit(data[i][2],1); unit_index = random(length(K1.fu))+1;
    lglat = process_complex_loglattice(K ,data[i][3]);                          \\
    cpct_units = cpct_from_loglattice(K, lglat, eps);
    reg1 = unscaled_determinant(K, lglat);

    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K.pol), log(abs(K.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K.pol), log(abs(K.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));


    REQ_COMPARE = ceil((poldegree(K.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);

    REQ_RIG = prec_rigorous(poldegree(K.pol), log(abs(K.disc)), log(infinity_norm(sumv)),log(abs(reg1))  );
    print("REQ_RIG ",ceil(REQ_RIG));
    default(realprecision, ceil(REQ_RIG));
    \\
    \\  This takes the log lattice and modifies it so that we get a sublattice.
    \\  The modification depends on the 'latticetype', see function compute_sublattice
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \\latticetype = 0;
    \\[lglat_new,modpair1]=compute_sublattice(lglat, OUTFILE1, latticetype);
    \\print(precision(lglat_new,10), "\nMODPAIR", modpair1);



    \\ Determine if you want to create a sublattice or just use the normal one.
    \\compute_sublattice(lglat, OUTFILE1, 1);
    \\lglat_new = lglat; lglat_new[,1] = 3*lglat_new[,1]; print(precision(lglat_new,10), "\nMODPAIR", modpair1);

    \\ If this line uncommented, then the input lattice is just the GRH-assumed unit lattice
    lglat_new = lglat;
    write(OUTFILE1, "\n--------------------------\n", i, " Field pol: ", K.pol, "Disc: ", K.disc, ".      Signature: ", K.r1, " ", K.r2);
    write(OUTFILE1, "\nRegulator: ", precision(reg1,10),"--------------------------precision value ", ceil(REQ_BSGS));


    \\write(OUTFILE1,"\n--- Field number ", i, " " , K.pol, "\nModified lglat ", precision(lglat_new,10));
    \\inputreg = unscaled_determinant(K,lglat_new);
    \\write(OUTFILE1," Input Regulator: ", precision(inputreg,10), "  Original Regulator: ", precision(reg1,10)  );
    n = poldegree(K.pol);
    logdisc = log(abs(K.disc));



    p1 = pmax_p1(n,logdisc, log(abs(reg1)) );
    p2 = pmax_p2(n, REQ_RIG, logdisc, log(abs(reg1)));
    g_n = giant_n(n, logdisc, REQ_BSGS,log(abs(reg1)));
    b_n = baby_n( n,logdisc,REQ_BSGS,log(abs(reg1)));
    \\pchoice = p1+p2;
    pchoice = p1;

    print("gn ",  ceil(g_n), "   bn", ceil(b_n));
    print("pn ", ceil(pchoice));
    balanceB = ((abs(reg1))^(1/3)) * (g_n*b_n)^(1/3);
    balanceB /= (pchoice^(2/3));

    maxnorm_index = 1;
    for(i=1, length(lglat_new),
        if(norml2(lglat_new[,i])> norml2(lglat_new[,maxnorm_index]),
            print("vector norms: ", precision( norml2(lglat_new[,i]),10));
            maxnorm_index = i;
        );
    );
    print("max vector norm: ", maxnorm_index, "   ",precision( sqrt(norml2(lglat_new[,maxnorm_index])),10) );



    balanceB = abs(log(reg1))*2^poldegree(K.pol)*reg1^(1/3);
    balanceB = min(reg1, balanceB);
    balanceB = min(balanceB, sqrt( norml2(lglat_new[,length(lglat_new)])  ) );


    write(OUTFILE1, "Chosen bound B ", balanceB);

    \\ Should define this based on the optimal balance of pohst and bsgs
    \\pohstB = 10000;
    \\pohstB = inputreg^(1/3);

    write(OUTFILE1, "B = ", precision(balanceB,10));
    print("Running Pohst Algorithm");                                           \\ lglat_new is the input lattice, pohst_out_lattice is the result after ruling out index divisors up to pohstB
    unitvector_cpct = cpct_from_loglattice(K, lglat_new, eps);                  \\ computation of compact reps
    tbefore = getabstime();
    pohst_out_lattice = log_pohst_pari(K,lglat_new,unitvector_cpct, balanceB, eps);

    stage1_units = cpct_from_loglattice(K, pohst_out_lattice,eps);
    tafter = getabstime();

    lptime = tafter-tbefore;
    \\ Just checking the regulator of the output from the p-maximization
    \\write(OUTFILE1,"Pohst Output Regulator: ", precision(outreg,10 ), ". Ratios: ", (modpair1[2]-inputreg/outreg)< eps);
    write(OUTFILE1, "pmax time ",precision(lptime,10), " In minutes: ", precision(lptime/60000.0,15) );

    print("Running BSGS Algorithm");
    default(realprecision, ceil(REQ_BSGS));
    \\print("REQ_BSGS ",floor(REQ_BSGS) );
    t9 = getabstime();
    bsgs_out= bsgs(K,stage1_units, balanceB, 1/2, eps, REQ_BSGS, OUTFILE1);
    t10 = getabstime();
    bsgstime = t10-t9;
    bsgs_out_lattice = log_lattice_from_compact_set(K, bsgs_out);
    outreg = unscaled_determinant(K,bsgs_out_lattice);

    \\write(OUTFILE1,"BSGS Output Regulator: ", precision(outreg,10 ), ". Ratios: ", (modpair1[2]-inputreg/outreg)< eps);
    write(OUTFILE1, "bsgs time ",precision(bsgstime,10), " In minutes: ", precision(bsgstime/60000.0,15) );
    write(OUTFILE1,"Overall time: ", precision(bsgstime+lptime , 10) , " In minutes: ", precision(bsgstime+lptime/60000.0,15) );

);
}
