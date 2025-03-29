/*
%Computing a simpler upper bound for N(n) and T(n)
%Here for Baby Stock: the pre

1. When n \le 30:
*/


Nfunc(n, rho) ={
  my(output = 0);
  \\output = ((n +30*sqrt(n)-1)^(n +30*sqrt(n)-1));
  \\output = output/( (n-1)^(n-1/2) );
  \\output = output/( (30*sqrt(n))^(30*sqrt(n)+1/2 ) );
  output = (1+2*sqrt(n)*exp(2*rho))^n;
  return(output);
}
Tfunc(n)={
  my(output = 0);
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
  return( max(10, ceil(4*n^2*logdisc +2*n^4 +(24*n^2-n^3)*log(n)))    );
};
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

REQ_REDUCTION(G)=
{
  my(n = poldegree(G.pol),
    ldisc = log(abs(G.disc))/log(2),
    m = G.r1 +G.r2,
    prec = 0
  );
  prec = (n^2+2*n+1+1/n)*ldisc + ((3*n^3 - 3*n^2 +3*n)/4)+ ((2*n^2+6*n+3)/2)*log(n)/log(2);
  prec += (m/2)*log(n)/log(2) + 3*m + 2*log(m)/log(2) -m*log(m-1)/log(2) + 8;

  \\oldprec = 2*ldisc + (n+1)^2*(10*log(n)/log(2) + 3*n*(n-1)/4 -(n-1/2)*log(n)/log(2)+ldisc + (n^2+1)*ldisc+ log(abs(G.disc)+2)/log(2) );
  \\print(ceil(oldprec), "  ", ceil(prec));
  return(ceil(prec));
}

base_2_log(val) =
{
  return (log(val)/log(2))
}

REQ_JUMP(G, v)=
{
  my(
      n = poldegree(G.pol),
      log_n = ceil(base_2_log(n)),
      ldisc = base_2_log(abs(G.disc)),
      m = G.r1 +G.r2,
      q_jump,small_delta, tv,
      prec = 0
  );
  small_delta = ((2/Pi)^G.r2)*sqrt(abs(G.disc));
  tv = floor(log(n*normlp(v)/base_2_log(small_delta)) / log(2))+1;
  prec = REQ_REDUCTION(G) + 2*tv + log_n + 2;

  return(prec)
}

REQ_BABY(G, reg1, v)=
{
  return(ceil(REQ_JUMP(G,v) + base_2_log(reg1)/2 ));
}

REQ_GIANT(G, reg1, v)=
{
  my(m = G.r1+G.r2);
  return (ceil(REQ_JUMP(G,v) + base_2_log(m-1)-1+ m*base_2_log(3)-m/2+base_2_log(reg1)));
}

REQ_RIGOROUS(G,v, p)=
{
  return( ceil(REQ_JUMP(G, v)+log(p)+log(poldegree(G.pol))) );
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

/********************************************/
/********************************************/
/*
New as of March 1, 2025
*/
gfunc(n, delta_K, detLambda)={
  my(output, logdeltaK, log_detLambda);
  logdeltaK = log(abs(delta_K))/log(2);
  log_detLambda = log(detLambda)/log(2);
  prec_q = n^2*logdeltaK + log_detLambda;
  output = (n^5)*log_detLambda^2;
  output += (n^2 + log_detLambda)*((n^5)*(prec_q^2)+ (2^n)*prec_q);
}

bfunc(n, delta_K, detLambda, rho)={
  my(output, test_exp, logdeltaK, log_detLambda);
  logdeltaK = log(abs(delta_K))/log(2);
  log_detLambda = log(detLambda)/log(2);
  prec_q = n^2*logdeltaK + log_detLambda;

  test_exp = n^9+logdeltaK^2 +Tfunc(n)*(n^3+n*logdeltaK);
  output = log_detLambda+n^2;
  output *= ((n^5)*(prec_q^2)+ (2^n)*prec_q);
  output += test_exp*Nfunc(n,rho);
  return(output);
}

pfunc1(n, delta_K, detLambda)=
{
  my(output, logdeltaK, log_detLambda,q_sat);
  log_detLambda = log(detLambda)/log(2);
  logdeltaK = log(abs(delta_K))/log(2);
  q_sat = n^2*logdeltaK + 2*log_detLambda+n^3;
  output = (q_sat^2)*n^5 + q_sat*(2^n);
  output *= (n^3+ n*log_detLambda);
}

pfunc2(n, delta_K, detLambda)=
{
  my(output, logdeltaK, log_detLambda,q_sat);
  log_detLambda = log(detLambda)/log(2);
  output= n^5+(n^3)*log_detLambda;
}
/********************************************/
/********************************************/

baby_n(n,logdisc,q,logdetlamp)={
  return( (n^4*(logdisc + logdetlamp +n^2)^2) * (n^7*(logdisc+logdetlamp+n^2)+(2^n) )+Nfunc(n)*TOtest(n,logdisc) );
}

bsgs_complexity(det_lambda, B, b_n, g_n) = {
  return(sqrt(det_lambda/B)*sqrt(b_n*g_n));
}

volumeB(det_lambda, B, g_n, b_n) ={
  return( sqrt((det_lambda*g_n)/(B*b_n) ) );
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
