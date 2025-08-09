
/**
* startvalue -- beginning value for the primes
* upperbound -- limit to the considered primes
* ap_count   -- number of primes to seek in arithmetic progression
*/

setInstanceVariables(readData)={
  K = nfinit(readData[1]);
  lglat = process_complex_loglattice(K ,readData[3]);
  reg1 = abs(matdet(lglat));
  r = K.r1+K.r2 -1;
  return([K, lglat, reg1,r]);
}

process_complex_loglattice(~K, ~lglat)={
    my(r, complexlogunit, lambda1, LambdaK);
    r = K.r1 + K.r2 -1;
    lambda1 = real(lglat);
    LambdaK = lambda1[1..r,];
    return(LambdaK);
}


prime_experiment(startvalue, upperbound, ap_count)={
    p = nextprime(startvalue);
    prime_limit = upperbound;
    values_per_prime = ap_count;
    while (p < prime_limit,

        ctr = 0;
        k = 1;
        q= 1;
        forprimestep(r = 1, +oo , p,
            ctr++;
            if(ctr >= ap_count,
                q = r;
                break;
            );
        );
        big_k = (q-1)/p; print("Prime p = ", p, ". Max k = ", big_k);
        threshold = log(k/values_per_prime)/log(p);
        if(threshold > 0.4,
          print(p, "  largest k = ", k, " Log_p(k)=", precision(log(k/values_per_prime)/log(p),10) );
        );
        p = nextprime(p+1);
    );
}

inertial_one_prime_experiment(OUTPUT_FILE, nf, startvalue, upperbound, ap_count)={
    K = nf;

    p = nextprime(startvalue);
    prime_limit = upperbound;
    values_per_prime = ap_count;
    maxlogp = 0;
    maxk = 0;
    maxp = 1;
    prime_counter = 1;

    while (p < prime_limit,


        ctr = 0;
        k = 1;
        q= 1;
        forprimestep(r = 1, +oo , p,
            prid1 = idealprimedec(K, r);
            for(i=1, length(prid1),
                if(prid1[i].f == 1,   \\ check if inertial degree is 1
                    \\print("inertial 1 ", prid1[i]);
                    ctr++;
                    break;
                ,
                    \\print("inertial >1 ", prid1[i]);
                );
            );

            if(ctr >= ap_count,
                q = r;
                break;
            );
        );
        big_k = (q-1)/p;
        if((prime_counter % 10000) == 0,
          print(p, ". Max k = ", big_k)
        );
        log_p_k = log(big_k)/log(p);

        if(log_p_k > maxlogp,
            maxlogp = log_p_k;
        );
        if(big_k > maxk,
          maxk = big_k;
          maxp = p;
        );

        \\threshold = log(k/values_per_prime)/log(p);
        \\if(threshold > 0.4,
        \\  print(p, "  largest k = ", k, " Log_p(k)=", precision(log(k/values_per_prime)/log(p),10) );
        \\);

        prime_counter++;
        p = nextprime(p+1);
    );
    write(OUTPUT_FILE, "Max_k= " ,maxk, " N=", ap_count," p in: ", startvalue, "-",upperbound,  " :: ",precision(log(nf.disc)/log(2),10), " , ", nf.pol );
}

{
    sigstring = "4-0";
    prefix = "../input/experiment-polynomials-";
    INPUT_FILE = concat(prefix, sigstring);
    OUTPUT_FILE = concat("Q-experiment-", sigstring);
    read(INPUT_FILE);

    startvalue = 1;
    endvalue = startvalue+1000000;
    N = 2;
    for(j=1, 8,
      for(i=26, length(data)-5,
          [K, lglat, reg1, urank] = setInstanceVariables(data[i]);
            inertial_one_prime_experiment(OUTPUT_FILE, K, startvalue, endvalue, N);
      );
      N+=2;
    );

}
