process_complex_loglattice(~K, ~lglat)={
    my(r, complexlogunit, lambda1, LambdaK);
    r = K.r1 + K.r2 -1;
    lambda1 = real(lglat);                                                       \\ equivalent to getting the log of the abs vals
    LambdaK = lambda1[1..r,];
    return(LambdaK);
}

verify_field(filen, r1, r2)={
	read(filen);
	print("There are ", length(data), " (", r1, ",", r2, ")", " fields in this file.
		 Verifying signatures...");
	for(i=1, length(data),
		K = nfinit(data[i][1]); \\print(K.r1);
		if(K.r1 !=r1 || K.r2 != r2,
				print (i, " BAD FIELD!! Expected Signature ", r1, " ", r2, "\n");
				print("Got:  ", K.r1, " ", K.r2, "  ", K.pol, " ", real(log(K.disc)/log(10))););
		print(floor(log(abs(K.disc))/log(2)), "  ", floor(log(abs(K.disc))/log(10)));
		\\print(i,  "  ", K.r1, "  ", K.pol, " ", precision(real(log(K.disc)/log(10)),10));
	);
}

verify_field_bitsize(filen, r1, r2)={
	read(filen);
	print("There are ", length(data), " (", r1, ",", r2, ")", " fields in this file.
		 Verifying signatures...");
	for(i=1, length(data),
		K = nfinit(data[i][1]); \\print(K.r1);

		lglat = process_complex_loglattice(K ,data[i][3]);
  		reg1 = abs(matdet(lglat));

		if(K.r1 !=r1 || K.r2 != r2,
				print (i, " BAD FIELD!! Expected Signature ", r1, " ", r2, "\n");
				print("Got:  ", K.r1, " ", K.r2, "  ", K.pol, " ", real(log(K.disc)/log(2))););
		print(K.disc, "   ", round(log(abs(K.disc))/log(2)), "  reg: ", round(log(abs(reg1))/log(2)));

	);
}
{


verify_field_bitsize("new-experiment-polynomials-0-2", 0,2);
/*
verify_field_bitsize("experiment-polynomials-0-2", 0,2);
verify_field_bitsize("experiment-polynomials-0-3", 0,3);
verify_field_bitsize("experiment-polynomials-0-4", 0,4);

verify_field_bitsize("experiment-polynomials-1-1", 1,1);
verify_field_bitsize("experiment-polynomials-1-2", 1,2);
verify_field_bitsize("experiment-polynomials-1-3", 1,3);

verify_field_bitsize("experiment-polynomials-2-1", 2,1);
verify_field_bitsize("experiment-polynomials-2-2", 2,2);

verify_field_bitsize("experiment-polynomials-3-0", 3,0);
verify_field_bitsize("experiment-polynomials-3-1", 3,1);

verify_field_bitsize("experiment-polynomials-4-0", 4,0);
verify_field_bitsize("experiment-polynomials-5-0", 5,0);
*/
/*
\\verify_field_bitsize("bsgs-r3-polynomials-0-4", 0,4);
\\verify_field_bitsize("experiment-polynomials-4-0", 4,0);
\\verify_field_bitsize("bsgs-r2-polynomials-5-0", 5,0);

verify_field_bitsize("bsgs-r2-polynomials-0-2", 0,2);
verify_field_bitsize("bsgs-r2-polynomials-0-3", 0,3);
verify_field_bitsize("bsgs-r2-polynomials-0-4", 0,4);
\\verify_field_bitsize("experiment-polynomials-0-4", 0,4);

verify_field_bitsize("bsgs-r2-polynomials-1-1", 1,1);
verify_field_bitsize("bsgs-r2-polynomials-1-2", 1,2);
verify_field_bitsize("bsgs-r2-polynomials-1-3", 1,3);

verify_field_bitsize("bsgs-r2-polynomials-2-1", 2,1);
verify_field_bitsize("bsgs-r2-polynomials-2-2", 2,2);

verify_field_bitsize("bsgs-r2-polynomials-3-0", 3,0);
verify_field_bitsize("bsgs-r2-polynomials-3-1", 3,1);

verify_field_bitsize("bsgs-r2-polynomials-4-0", 4,0);
*/

}
