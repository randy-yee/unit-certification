
// Input:
// a number field object,
// a vector beta consisting of the torsion unit and an independent set of units
// prim a prime number (t_INT)
// bad = 1, gen_1


// latest version as of November 13
// If the cpct flagi s set to 1, assume that a compact representation is of the form
// [List(), vector]
//
GEN
pari_prime_check(GEN nf, GEN beta, GEN prim, GEN bad, GEN cpct)
{
  pari_sp av = avma;
  ulong p = itou(prim);
  long lb = lg(beta), rmax = lb - 1, lcpct;
  GEN M, vQ, L;
  ulong q, idealctr;
  forprime_t T;

  if (DEBUGLEVEL>3) {pari_printf("Primecert begins with p = %lu,   bad = %Ps .", p,bad);output(beta);}
  idealctr = 0;
  if (p == 2)
    L = cgetg(1,t_VECSMALL);
  else
    L = mkvecsmall(p);
  (void)u_forprime_arith_init(&T, 1, ULONG_MAX, 1, p);
  M = cgetg(lb,t_MAT); setlg(M,1);

  while ((q = u_forprime_next(&T)))
  {
    GEN qq, gg, og;
    long lQ, i, j, k;
    ulong g, m;
    if (!umodiu(bad,q)) continue;

    qq = utoipos(q);
    if(DEBUGLEVEL > 2) pari_printf("Current q: %Ps \n", qq);
    vQ = idealprimedec_limit_f(nf,qq,1);
    lQ = lg(vQ); if (lQ == 1) continue;

    /* cf rootsof1_Fl */
    g = pgener_Fl_local(q, L);
    m = (q-1) / p;
    gg = utoipos( Fl_powu(g, m, q) ); /* order p in (Z/q)^* */
    og = mkmat2(mkcol(utoi(p)), mkcol(gen_1)); /* order of g */

    if (DEBUGLEVEL>3) err_printf(" Generator of (Zk/q)^*: %lu\n", g);
    for (i = 1; i < lQ; i++)
    {
      idealctr++;
      GEN C = cgetg(lb, t_VECSMALL);
      GEN Q = gel(vQ,i); /* degree 1 */
      GEN modpr = zkmodprinit(nf, Q);

      long r;
      GEN t;
      for (j = 1; j < lb; j++)
      {
        if(itou(cpct) == 1){

          GEN cpctdenom = gel(gel(beta,j),2);
          GEN cpctnum = gel(gel(gel(beta,j),1),2);
          if(DEBUGLEVEL > 2) {
            pari_printf("Cpct Rep: %Ps   %Ps \n", cpctnum, cpctdenom);
            pari_printf("Q = %Ps \n", Q);
          }

          lcpct = lg(cpctnum);
          t = nf_to_Fp_coprime(nf, gen_1, modpr);
          GEN intermediate;
          if(DEBUGLEVEL > 2) pari_printf("t: %Ps \n",t);
          for(k=1; k < lcpct; k++){

            intermediate = nf_to_Fp_coprime(nf, gel(cpctnum,k), modpr);
            //pari_printf("Numerator: %lu   denominator: %lu\n", intermediate[2], itos(gel(cpctdenom,k)));
            intermediate = utoipos(Fl_div(intermediate[2], (itos(gel(cpctdenom,k))%q) , q));
            //pari_printf("done inverting");
            t = utoipos(Fl_mul(t[2],t[2], q));
            t= utoipos(Fl_mul(t[2], intermediate[2], q));


          }
          if (DEBUGLEVEL>2) {pari_printf("t =  %Ps .\n", t);}
        }
        else{
          t = nf_to_Fp_coprime(nf, gel(beta,j), modpr);
        }
        t = utoipos( Fl_powu(t[2], m, q) );
        C[j] = itou( Fp_log(t, gg, og, qq) ) % p;
      }
      r = lg(M);
      gel(M,r) = C; setlg(M, r+1);
      if (DEBUGLEVEL>3){pari_printf("Ideal ctr: %lu ,  r = %ld,    rmax = %ld   2lb = %ld\n", idealctr, r, rmax, 2*lb);}

      if (Flm_rank(M, p) != r) {
        setlg(M,r);
        if(idealctr > 2*lb){avma = av; return gen_0;}
        continue;
      }

      if (DEBUGLEVEL>2)
      {
        if (DEBUGLEVEL>3)
        {
          err_printf("       prime ideal Q: %Ps\n",Q);
          err_printf("       matrix log(b_j mod Q_i): %Ps\n", M);
        }
        err_printf("       new rank: %ld\n",r);
      }
      if (r == rmax) {
        if(DEBUGLEVEL>3){pari_printf("\nNumber of ideals considered: %lu \n\n", idealctr);}
        avma = av;
        return gen_1;
      }
    }
  }
  if(DEBUGLEVEL>3){pari_err_BUG("primecertify");}
  return gen_2;
}



// this second function has an additional parameter paramN which controls
// the number of ideal the check before quitting
GEN
pari_prime_check_N(GEN nf, GEN beta, GEN prim, GEN bad, GEN paramN, GEN cpct)
{
  pari_sp av = avma;
  ulong p = itou(prim);
  ulong parN = itou(paramN);
  long lb = lg(beta), rmax = lb - 1, lcpct;
  GEN M, vQ, L;
  ulong q, idealctr;
  forprime_t T;

  if (DEBUGLEVEL>3) {pari_printf("Primecert begins with p = %lu,   bad = %Ps .", p,bad);output(beta);}
  idealctr = 0;
  if (p == 2)
    L = cgetg(1,t_VECSMALL);
  else
    L = mkvecsmall(p);
  (void)u_forprime_arith_init(&T, 1, ULONG_MAX, 1, p);
  M = cgetg(lb,t_MAT); setlg(M,1);

  while ((q = u_forprime_next(&T)))
  {
    GEN qq, gg, og;
    long lQ, i, j, k;
    ulong g, m;
    if (!umodiu(bad,q)) continue;

    qq = utoipos(q);
    if(DEBUGLEVEL > 2) pari_printf("Current q: %Ps \n", qq);
    vQ = idealprimedec_limit_f(nf,qq,1);
    lQ = lg(vQ); if (lQ == 1) continue;

    /* cf rootsof1_Fl */
    g = pgener_Fl_local(q, L);
    m = (q-1) / p;
    gg = utoipos( Fl_powu(g, m, q) ); /* order p in (Z/q)^* */
    og = mkmat2(mkcol(utoi(p)), mkcol(gen_1)); /* order of g */

    if (DEBUGLEVEL>3) err_printf(" Generator of (Zk/q)^*: %lu\n", g);
    for (i = 1; i < lQ; i++)
    {
      idealctr++;
      GEN C = cgetg(lb, t_VECSMALL);
      GEN Q = gel(vQ,i); /* degree 1 */
      GEN modpr = zkmodprinit(nf, Q);

      long r;
      GEN t;
      for (j = 1; j < lb; j++)
      {
        if(itou(cpct) == 1){

          GEN cpctdenom = gel(gel(beta,j),2);
          GEN cpctnum = gel(gel(gel(beta,j),1),2);
          if(DEBUGLEVEL > 2) {
            pari_printf("Cpct Rep: %Ps   %Ps \n", cpctnum, cpctdenom);
            pari_printf("Q = %Ps \n", Q);
          }

          lcpct = lg(cpctnum);
          t = nf_to_Fp_coprime(nf, gen_1, modpr);
          GEN intermediate;
          if(DEBUGLEVEL > 2) pari_printf("t: %Ps \n",t);
          for(k=1; k < lcpct; k++){

            intermediate = nf_to_Fp_coprime(nf, gel(cpctnum,k), modpr);
            //pari_printf("Numerator: %lu   denominator: %lu\n", intermediate[2], itos(gel(cpctdenom,k)));
            intermediate = utoipos(Fl_div(intermediate[2], (itos(gel(cpctdenom,k))%q) , q));
            //pari_printf("done inverting");
            t = utoipos(Fl_mul(t[2],t[2], q));
            t= utoipos(Fl_mul(t[2], intermediate[2], q));


          }
          if (DEBUGLEVEL>2) {pari_printf("t =  %Ps .\n", t);}
        }
        else{
          t = nf_to_Fp_coprime(nf, gel(beta,j), modpr);
        }
        t = utoipos( Fl_powu(t[2], m, q) );
        C[j] = itou( Fp_log(t, gg, og, qq) ) % p;
      }
      r = lg(M);
      gel(M,r) = C; setlg(M, r+1);
      if (DEBUGLEVEL>3){pari_printf("Ideal ctr: %lu ,  r = %ld,    rmax = %ld   2lb = %ld\n", idealctr, r, rmax, 2*lb);}

      if (Flm_rank(M, p) != r) {
        setlg(M,r);
        if(idealctr > parN){avma = av; return gen_0;}
        continue;
      }

      if (DEBUGLEVEL>2)
      {
        if (DEBUGLEVEL>3)
        {
          err_printf("       prime ideal Q: %Ps\n",Q);
          err_printf("       matrix log(b_j mod Q_i): %Ps\n", M);
        }
        err_printf("       new rank: %ld\n",r);
      }
      if (r == rmax) {
        if(DEBUGLEVEL>3){pari_printf("\nNumber of ideals considered: %lu \n\n", idealctr);}
        avma = av;
        return gen_1;
      }
    }
  }
  if(DEBUGLEVEL>3){pari_err_BUG("primecertify");}
  return gen_2;
}
