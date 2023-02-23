

verify_field(filen, r1, r2)={
	read(filen);
	print("There are ", length(data), " (", r1, ",", r2, ")", " fields in this file.");
	for(i=1, length(data),
		K = nfinit(data[i][1]); \\print(K.r1);
		if(K.r1 !=r1 || K.r2 != r2, print ("\nBAD FIELD!!!\n"); print(i,  "  ", K.r1, "  ", K.pol, " ", real(log(K.disc)/log(10))););
		print(floor(log(abs(K.disc))/log(10)));
		\\print(i,  "  ", K.r1, "  ", K.pol, " ", precision(real(log(K.disc)/log(10)),10));
	);
}
{
/*
verify_field("test-poly-0-2.gp", 0,2);
verify_field("test-poly-0-3.gp", 0,3);
verify_field("test-poly-0-4.gp", 0,4);


verify_field("test-poly-1-1.gp", 1,1);
verify_field("test-poly-1-2.gp", 1,2);
verify_field("test-poly-1-3.gp", 1,3);

verify_field("test-poly-2-1.gp", 2,1);
verify_field("test-poly-2-2.gp", 2,2);

verify_field("test-poly-3-0.gp", 3,0);
verify_field("test-poly-3-1.gp", 3,1);
verify_field("test-poly-4-0.gp", 4,0);
verify_field("test-poly-4-1.gp", 4,1);
*/
verify_field("polynomial-0-3.gp", 0,3);
verify_field("polynomial-0-2.gp", 0,2);
verify_field("polynomial-0-4.gp", 0,4);
verify_field("polynomial-1-1.gp", 1,1);
verify_field("polynomial-1-2.gp", 1,2);
verify_field("polynomial-2-1.gp", 2,1);
verify_field("polynomial-2-2.gp", 2,2);
verify_field("polynomial-4-0.gp", 4,0);
verify_field("polynomial-3-1.gp", 3,1);
verify_field("polynomial-1-3.gp", 1,3);
\\verify_field("test-poly-5-0.gp", 5,0);
\\verify_field("test-poly-6-0.gp", 6,0);

verify_field("test-poly-0-5.gp", 0,5);
\\verify_field("test-poly-0-6.gp", 0,6);
\\verify_field("test-poly-1-4.gp", 1,4);

}
