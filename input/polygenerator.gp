\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Functions for generation polynomials with a particular discriminant size
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\p1000

default(parisize,"12G")

print("Loaded polygenerator");
/************************************************************/
\\ OUTPUT: a random approximation of a real number,
\\ interval is a range in which you want the number, dec is a decimal number whose precision will dictate the output
\\ precision

random_rroot(interval, dec)={
	my(intpart, floatpart);
	intpart = random(interval);
	floatpart = random(dec);

	return(intpart+floatpart);
}
/************************************************************/
/************************************************************/
\\ OUTPUT: a random approximate complex number
\\ intervalre and intervalim indicate the range you want the integer part of the real and complex coefficients respectively
\\ dec is as before, just some decimal to indicate the precision.
random_croot(intervalre, intervalim, dec)={

  my(intpart, floatpart, croot);

  intpart = random(intervalre);
	floatpart = random(dec);
	croot = intpart+floatpart;

	intpart = random(intervalre);
	floatpart = random(dec);
	croot = croot + I*(intpart+floatpart);
	return(croot);
}

\\ OUTPUT: a random approximate complex number
\\ intervalre and intervalim indicate the range you want the integer part of the real and complex coefficients respectively
\\ dec is as before, just some decimal to indicate the precision.
random_croot_norm(normie, dec)={

  my(intpart, floatpart, croot, imag_bound);

  intpart = random( floor(sqrt(normie)) );
	floatpart = random(dec);
	croot = intpart+floatpart;

	intpart = ceil(sqrt((normie - croot^2)));
	imag_bound = [floor(intpart/10), intpart*10];
	intpart = random( imag_bound );
	floatpart = random(dec);
	croot = croot + I*(intpart+floatpart);
	return(croot);
}
/************************************************************/
/************************************************************/
\\ OUTPUT: generates an index 1 real cubic polynomial.
\\ magnitude indicates how large of a discriminant we want, as a power of 10
generate_real_cubic(magnitude)={

    my(r1,r2,r3, candidate, flag, distance1, distance2, sgn);
    flag = 0;
    magnitude = ceil(magnitude/4);
    if(magnitude%2 == 0, distance1 = magnitude/2; distance2 = distance1;);
    if(magnitude%2 == 1, distance1 = magnitude\2 ; distance2 = distance1+1);

    \\ distance one indicates how far to choose the 2nd root from r1, distance2 is for choosing r3
    distance1 = 10^distance1;
    distance2 = 10^distance2;
    while(flag == 0,
      r1 = random_rroot([0,1000],1.0);
      sgn = random(2);
      r1 = (-1^sgn)*r1;

      sgn = random(2);
      \\case 0: r1 will be the smallest absval root
      if(sgn == 0,
        sgn = (-1)^random(2);
        r2 = sgn*random_rroot([distance1, 2*distance1],1.0);
        r3 = sgn*random_rroot([distance1+distance2,2*distance1+distance2 ],1.0);
      );
      \\case 2: r1 will be the middle abs val root
      if(sgn== 1,
        r2 = random_rroot([distance1, 2*distance1],1.0);
        r3 = -random_rroot([distance2,2*distance2 ],1.0);
      );

      \\print(r1,"   ",r2, "   ",r3);
      candidate = round((x-r1)*(x-r2)*(x-r3) );
      \\print(candidate, "  " , poldisc(candidate));
      flag = 1;
      if(flag == 1,
          K = nfinit(candidate);
          if (K.index == 1,
            flag = 1,
            flag = 0
          );
      )
    );
    return(candidate);
}
/************************************************************/
/************************************************************/
\\ OUTPUT: A mixed quartic with discriminant 'close' to (10^magnitude)
\\ I can probably make these more general, but I wasn't too sure how to generalize the restriction on how to pick the roots
\\ generates a mixed quartic who |discriminant| is close to 10^mag
\\ and whose index is 1.
mixed_quartic(magnitude)={

    my(r1,r2,s1, candidate, flag, intervalre, real_distance, sgn, range);

    flag = 0;

    magnitude = ceil(magnitude/2);
    print("magnitude", magnitude);
    real_distance = 10^(ceil(magnitude/4));
    real_distance = real_distance^0.5;
    range = real_distance/2;
    while(flag == 0,
      sgn = (-1)^(random(2));
      r1 = sgn*random_rroot([0,range],1.0);
      r2 = sgn*random_rroot([range+real_distance,real_distance+2*range],1.0);
      if(r1<r2, intervalre = [r1,r2], intervalre = [r2,r1]);
      \\print(intervalre);

      sgn = (-1)^(random(2));
      s1 = random_croot(intervalre, [sgn*realdistance/2, sgn*(realdistance/2) + range], 1.0);
      s2 = conj(s1);

      candidate = round((x-r1)*(x-r2)*(x-s1)*(x-s2) );
      /*print(candidate, "  " , poldisc(candidate));
      print("factors:");
      print(r1-r2);
      print(2*(s1-s2));
      print(r1-s1, " ", norm(r1-s1));
      print(r2-s1, " ", norm(r2-s1));
      print((r1-r2)*2*(s1-s2)*norm(r1-s1)*norm(r2-s1) );
      */
      flag = 1;
      if(flag == 1,
          K = nfinit(candidate);
          if (K.index == 1,
            flag = 1,
            flag = 0
          );
      )
    );
    return(candidate);
}

/************************************************************/
/************************************************************/
\\ INCOMPLETE - This function should allow you to select a signature, and then produce an irreducible poly whose
\\ 							discriminant is close to
\\ OUTPUT:
random_poly(r,s, magnitude)={

    my(r1,s1, candidate, flag, intervalre, offset, sgn, range);
    flag = 0;

		\\ counting terms in the discriminant formula
		\\ separated as real products, z*conj(z), and real x imag products + z*w where z, w are complex non-conjugate
    crossterms = (binomial(2*s,2)-s + r*s*2);																		\\ real x imag, and z*w terms
    realterms = binomial(r,2);																									\\ r1*r2, real roots
    imterms = s;																																\\ z * conj(z)


		\\ in the totally real case, choose roots so that they are spaced apart somewhat evenly
		if( s == 0,
			my(extra_factor = 1);
			for(i=2, r,
				extra_factor*= ( i^(r-i) );
			);
			extra_factor = round(log(abs(extra_factor))/log(10));
			magnitude = (magnitude-extra_factor)/(2*realterms);
		);


		if(s != 0,
			magnitude = (magnitude/2);
			magnitude = ( (magnitude)/(realterms+imterms+crossterms) );
		);

  	\\  print(magnitude, " terms ", (realterms+imterms+crossterms));

    range = (10^(magnitude));
    intervalre = [-range,range];
    intervalim = [10^(magnitude-1), 10^(magnitude-1) + range];
    flag = 0;
    while(flag == 0,

    realroots = [];
    complexroots = [];

		\\
    if(s != 0,
	    for(i = 1, s,
				s1 = random_croot_norm(range, 1.0);
	      s1 = random_croot(intervalre, intervalim, 1.0);
	      complexroots = concat(complexroots, s1);
				\\intervalre+= [10^magnitude, 10^magnitude];
				\\intervalim+= [10^magnitude, 10^magnitude];
	      \\print(norm(s1));
	    );
    );

    if(r!=0,
    intervalre = [0, range];
    for(i = 1, r,
      r1 = random_rroot(intervalre, 1.0);
      realroots = concat(realroots, r1);
      intervalre+= [10^magnitude, 10^magnitude];
      \\print(abs(r1));
    );
    offset = random(r+1)*magnitude;

    );

	    for(i=1,r,
	      realroots[i] = realroots[i]-offset
			);

      candidate = 1;

      for(i =1, r+s,
          if(i < r+1,
          candidate*=(x-realroots[i]);,
          candidate*=(x-complexroots[i-r])*(x-conj(complexroots[i-r]));
          );
      );
      candidate = round(candidate);

      \\print(candidate, "  " , poldisc(candidate));
      if (polisirreducible(candidate), flag = 1;);
      if(flag == 1 ,
          K = nfinit(candidate);
          if (K.index == 1,
            flag = 1,
            flag = 0
          );
					if(K.r1 != r, flag = 0);
      )

    );
    return(candidate);
}

{

writefile = "stest-poly-4-1.gp";
maxreal = 8;
maxcomplex = 2;

\\ testing new complex root choosing function
\\for(i=1, 500, croot = random_croot_norm(10^10, 1.0); print(precision(croot,10), "  " round(norm(croot))) );


print("Gathering polynomials");
disc_cap = 41;
for(r =3, 3,
	for(s = 0, 0,
		writefile = concat(["test-poly-", Str(r),"-",Str(s)]);
		\\if(r+2*s < 7 || r+2*s > 14 || r+s-1 > 6, ,
		if(0, ,
				write(writefile, "\\\\ Signatures ", r, " ",s);
				write(writefile, data,Str(r),"_"Str(s), " = [\\" );
				discsize = 40;
				while(discsize<disc_cap,
						tally = 0;
						while(tally < 10,
							pol1 = random_poly(r,s, discsize);
							K1 = bnfinit(pol1, 1);
							if( (  abs( log(abs( poldisc(pol1) ) )/log(10)-discsize )<1.1   )&& K1.clgp.no == 1,
								tally+=1;
								write(writefile, "[" , pol1, ", " , poldisc(pol1), ", \\\n",  K1[3]  , "], \\" ); );
						);

						discsize+=5;
				);
				write(writefile, "];");
		);
	);
);
}


/*
{for(i=1, 10,
	pol1 = random_poly(0,2, 29);
	print(pol1, ",    \\\\", poldisc(pol1), "      ", precision(log(abs(poldisc(pol1)))/log(10),10) );
)
}
{for(i=1, 10,
	pol1 = random_poly(0,3, 30);
	print(pol1, ",    \\\\ ", poldisc(pol1), "      ", precision(log(abs(poldisc(pol1)))/log(10),10) );
)
}

{for(i=1, 10,
	pol1 = random_poly(0,4, 30);
	print(pol1, ",    \\\\ ", poldisc(pol1), "      ", precision(log(abs(poldisc(pol1)))/log(10),10) );
)
}
{for(i=1, 10,
	pol1 = random_poly(0,5, 30);
	print(pol1, ",    \\\\ ", poldisc(pol1), "      ", precision(log(abs(poldisc(pol1)))/log(10),10) );
)
}
*/
\\for(i=1, 10, f = generate_real_cubic(24);print("[",polcoeff(f,3),",",polcoeff(f,2),",",polcoeff(f,1),",",polcoeff(f,0) "],  \\\\ ", poldisc(f)););
\\print("];")
\\testpoly = mixed_quartic(20);
\\{for(i=1,50,
\\f = mixed_quartic(30);
\\print("[",polcoeff(f,4), ",",  polcoeff(f,3),",",polcoeff(f,2),",",polcoeff(f,1),",",polcoeff(f,0) "], \\\\ ", poldisc(f), " " exponent(poldisc(f)) );
\\);
\\}

\\ f = generate_real_cubic(); print("[",polcoeff(f,3),",",polcoeff(f,2),",",polcoeff(f,1),",",polcoeff(f,0) "  \\\\ ", poldisc(f));
\\print("Double checking, is it irreducible? ", polisirreducible(f), ". And is it index one?  ", nfinit(f).index)
