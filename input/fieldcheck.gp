\\read("test-poly-1-1.gp"); data = data1_1;
\\read("test-poly-3-0.gp"); data = data3_0;
\\read("test-poly-4-0.gp"); data = data4_0;
read("test-poly-2-1.gp"); data = data2_1;
\\read("test-poly-0-2.gp"); data = data0_2;

\\read("test-poly-1-2.gp"); data = data1_2;
\\read("test-poly-3-1.gp"); data = data3_1;
\\read("test-poly-5-0.gp"); data = data5_0;
\\read("test-poly-0-3.gp"); data = data0_3;
\\read("test-poly-2-2.gp"); data = data2_2;

\\read("test-poly-4-1.gp"); data = data4_1;
\\read("test-poly-0-4.gp"); data = data0_4;
\\read("test-poly-1-3.gp"); data = data1_3;
{

r1 =4;
for(i=1, length(data),
	K = nfinit(data[i][1]); \\print(K.r1);
	\\if(K.r1 !=r1, print ("\nBAD FIELD!!!\n"); print(i,  "  ", K.r1, "  ", K.pol, " ", real(log(K.disc)/log(10))););
	print(i,  "  ", K.r1, "  ", K.pol, " ", precision(real(log(K.disc)/log(10)),10));
);
}
