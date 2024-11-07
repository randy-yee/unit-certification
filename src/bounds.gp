/*
%Computing a simpler upper bound for N(n) and T(n)
%Here for Baby Stock: the pre

1. When n \le 30:
*/
/*
n=30;
logdetlamp=20;
q=max(4*n^2*logdel+2*n^4, 2^n/n^4) + n^2 +logdetlamp);
exfun(x)= x^(x-1/2);
N(n)=exfun(n +30*sqrt(n)-1)/(exfun(n-1)*exfun(30*sqrt(n)))
T(n)=n^(n/1+2)*(3/2)^((n^2-n)/2)
expT(n)=n/2*(log(n)+6/n*log(n)+log(3/2)*n-log(3/2))
T(n)=2^(expT(n))

for(i=1,30,print(i, ", ", T(i)*1.0, ",  ", N(i)*1.0))
for(i=1,30,print(i, ", ", expT(i)*1.0, ",  ", i^2*1.0))
*/

Nfunc(n) ={
  my(output =0);
  output = ((n +30*sqrt(n)-1)^(n +30*sqrt(n)-1));
  output = output/( (n-1)^(n-1/2) );
  output = output/( (30*sqrt(n))^(30*sqrt(n)+1/2 ) );
  return(output);
}
Tfunc(n)={
  my(output =0);
  output = n^(n/2 + 2);
  output *= (3/2)^((n^2-n)/2);
  return(output);
}

idealPrecision(G, ideal, maxnorm)=
{
  my(log2 = log(2), rho);
  rho = poldegree(G.pol)^2;
  rho = rho*log(4*(denominator(ideal)^2)*maxnorm)/log2 +2;
  return(log(denominator(ideal))/log2 + 2*log(maxnorm+1)/log2 + 3 +rho);
}



B = 1;
p = B;
indexLambda = 1;


prec_reduction(n, logdisc)={
  return( max(1, ceil(4*n^2*logdisc +2*n^4 +(24*n^2-n^3)*log(n)))    );
};

reduction_complexity(n,q)={
  return(n^5*q^2+q*2^n);
}



prec_jump(n, logdisc, qout, tv)={
  my(precision_bound = 0);
  precision_bound = ceil(max( prec_reduction(n, logdisc)+logdisc+log(n), qout));
  return( precision_bound + 2*tv+2);
}

\\ q should be prec_jump
jump_complexity(n,q, loginfv) = {
  (n^5*q^2+q*2^n)*loginfv;
}


prec_compact(degree, logdisc,loginfv)={

  expression = ceil( (4*(degree^2))*logdisc + 2*(degree^4) -(degree^3 -24*(degree^2) -6)*log(degree) + 7 + (degree^2+2)*loginfv);
  return(expression);
}

prec_reduce(G)=
{
  my(n = poldegree(G.pol),
    ldisc = log(abs(G.disc))/log(2),
    m = G.r1 +G.r2,
    prec = 0
  );
  prec = (n^2+2*n+1+1/n)*ldisc + ((3*n^3 - 3*n^2 +3*n)/4)+ ((2*n^2+6*n+3)/2)*log(n)/log(2);
  prec += (m/2)*log(n)/log(2) + 3*m + 2*log(m)/log(2) -m*log(m-1)/log(2) + 8;

  oldprec = 2*ldisc + (n+1)^2*(10*log(n)/log(2) + 3*n*(n-1)/4 -(n-1/2)*log(n)/log(2)+ldisc + (n^2+1)*ldisc+ log(abs(G.disc)+2)/log(2) );
  print(ceil(oldprec), "  ", ceil(prec));
  return(ceil(prec));
}
prec_baby(n,log_disc, infsumt)={
  expression = (4*(n^2)*log_disc +2*n^4+(-n^3+24*n^2+6)*log(n) +7 +(n^2+4)*log(infsumt + (sqrt(n)/4)*log_disc) );
}
prec_giant(n, logdisc, logdetlamp, infsumu)={
  expression = (4*n^2 +1)*logdisc +2*n^4 -(n^3-24*n^2)*log(n);
  expression += (n^2+4)*log(infsumu +3);
}

prec_bsgs(fdegree, logdisc, inf_sumv)=
{
  return(ceil(max(prec_baby(fdegree,logdisc, inf_sumv), prec_giant(fdegree, logdisc, 1, inf_sumv) )));
}
\\ membership test complexity
TOtest(n, logdisc)={
  my(expression);
  expression = (n^5)*(n^2+logdisc)^2 + Tfunc(n)*(n^2 + logdisc);
  return(expression);
}

heuristic_complexity(n,p,logdetlamp, logdisc)={
    expression = n^2*p^(0.5) + n*(n^2+logdetlamp +log(logdisc)*(n*log(p)) );
}

giant_n(n,logdisc,q,logdetlamp)={
  return((n^2)*(n^2 +logdisc+logdetlamp )*(n^7*(n^2+logdisc+logdetlamp+2^n) )*(logdetlamp+n^2) );
}

baby_n(n,logdisc,q,logdetlamp)={
  return( (n^4*(logdisc + logdetlamp +n^2)^2) * (n^7*(logdisc+logdetlamp+n^2)+(2^n) )+Nfunc(n)*TOtest(n,logdisc) );
}

bsgs_complexity(detlamp, B, b_n, g_n) = {
  return(sqrt(detlamp/B)*sqrt(b_n*g_n));
}

volumeB(detlamp, B, g_n, b_n) ={
  return( sqrt((detlamp*g_n)/(B*b_n) ) );
}

prec_pthroot(n, B, logdisc, loginfv, logpsi_eta)={
  expression = (4*n^2+1)*logdisc +2*n^4+(-n^3+24*n^2)*log(n)+log(n) +(2/p)*logpsi_eta*loginfv+2;
  return(expression);
}

prec_rigorous(n, logdisc, logsumv, logdetlamp)={
  expression = (4*n^2+1)*logdisc +2*n^4 -(n^3-24*n^2-4)*log(n) +(n^2+2)*logsumv+2;
  print("REQ_RIG ",ceil(expression));
  return(expression);
}

pmax_p1(n, logdisc, logdetlamp)={
  return(n^3*(n^2+log(logdisc)+logdetlamp +1)+n^2 );
}

pmax_p2(n, q, logdisc, logdetlamp)={
  return( logdetlamp*( (n^6*q^2 +n*q*2^n)*(n^2+logdetlamp) +n^3 + (3*logdisc +n^2 +n)^n ));
}


{
/*
forstep(logdisc = 20, 50, 5,
  print("logdisc = ", logdisc);
  logdetlamp = sqrt(logdisc);
  volB = logdetlamp/B;
  volG = logdetlamp/B;

  for(n =3, 20,
    print("n =  ", n);
    print(prec_rigorous(n, 10, logdetlamp));
    \\q=max(4*n^2*logdisc+2*n^4, 2^n/n^4) + n^2 +logdetlamp;
    \\baby_complexity = volB*(n^5*q^3 + Nfunc(n)*(n^5*q^2 + Tfunc(n)*q )  );

    \\giant_complexity = volG*n^5*q^2*(logdetlamp+n^2);

    \\pohst_complexity = indexLambda * r*n^5*q^2*(n^2 + logdetlamp) + r*n^2*p*log(p) + (3*logdisc+n^2+n);
    \\print("- baby stock : ",baby_complexity);
    \\print("- giant step : ",giant_complexity);
    \\print("-     pohst : ",pohst_complexity);
  );
);
*/
}













n=50;
logdel=20;
q=4*n^2*logdel+2*n^4
[2^n*n*q, n^5*q^2, 2^n*n*q/(n^5*q^2*1.0)]
