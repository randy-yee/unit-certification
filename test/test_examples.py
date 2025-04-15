{
read("src/BabyStepGiantStep.py");
read("src/CompactRepresentation.py");
read("src/bounds.gp")

}

{
 \\K1 = bnfinit(x^{3} + 30* x - 360,1);nfalgtobasis(K1,K1.fu[1])
\\[7865950839478943661434907500541947030760361114794579294314873617, 1005351500390959543134043279054876965345304035388484825744600045, 1052770579104607884326123629972673238114829493466401750074041298]~
\\break> K1.zk
\\[1, x, 1/6*x^2 + 3]
    \\K1= nfinit(x^3 - 67*x^2 + 2032*x -2053);
    \\K2= bnfinit(x^3 - 67*x^2 + 2032*x -2053);
    \\K1 = nfinit(x^3 + 35*x - 92,1);
    \\K2 = bnfinit(x^3 + 35*x - 92,1);
    K1 = nfinit(x^3 - 161*x^2 + 7503*x - 99998);
    K2 = bnfinit(K1.pol,1);

    lglat = get_log_lattice_bnf(K2);


    reg1 = get_abs_determinant(lglat);
    B = 1;
    scanRadius =0.5;

    mat1 = embed_real(K1, K1[5][1]);
    short = qfminim(mat1~*mat1,,,2);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    print(K1.zk);
    print(K2.fu[1]);
    print(cpct_units[1]);
    breakpoint();
    \\bsgs_output= bsgs(K1,cpct_units, B, 18, scanRadius, eps,20,"alltest.txt");
    \\real matrix example and also check  LLL reduced basis
    rMat = [30.076166,
          -630.220923,
         27154.200950;
             0.257871,
             2.830679,
          -116.280533;
0,
             7.630279,
           129.365016];

}
{
/*
Error!: x^4 + 144*x^3 + 8399*x^2 + 147103*x + 2274843
[1, 0, 1/3, 11144165/43325466; 0, 1/3, 0, 30764695/158860042; 0, 0, 1/3, 7614331/43325466; 0, 0, 0, 1/476580126] [1, 0, 0, 40185149/147960486; 0, 1, 0, 1941326491/4290854094; 0, 0, 1, 725867059/4290854094; 0, 0, 0, 1/4290854094]
[-356008888.695695209390619870481335332775686667907412453123438061887602590349644060552236624062061309814453125, 356008910.8754468489735166473139049701749044453096252829299233110382781269441362692030782839857433763009451342074396148206228473261265865534852537166443653404712677001953125]
Assert failed. 0 does not evaluate to true
  ***   at top-level: ...mials-",if(manual_select,run_bsgs_experiment(s
  ***                                             ^---------------------
  ***   in function run_bsgs_experiment: ...getabstime();bsgs_output=
  ***   bsgs(K,cpct_units,sca
  ***   ^---------------------
  ***   in function bsgs: ...ge=="LOG",lattice_lambda=incremental_giant_ste
  ***                                                 ^---------------------
  ***   in function incremental_giant_steps: ...TE==0)&&ctr>0,base_value=
  ***   compact_rep_full_inpu
  ***   ^---------------------
  ***   in function compact_rep_full_input: ..." ",alphaOK,"\n",alpha));
  ***   GP_ASSERT_TRUE(idealB
  ***   ^---------------------
  ***   in function GP_ASSERT_TRUE: ...s not evaluate to true");breakpoint();)
  ***                                                           ^--------------
... skipping file 'experiment-babystock02.py'
*/
}

{
/*
Error!: x^3 - 361*x^2 + 327940*x - 52326638
[1, 1/2, 3932957/3984422; 0, 1/2, 527147/1992211; 0, 0, 1/3984422] [1, 0, 87195/8026591; 0, 1, 5779773/8026591; 0, 0, 1/8026591]
[-16102295.652427667985009004305442938277408356531785199289539066940099218405713088486663764342665672302246093750000000000000000000,
16102311.5506981307675074678186497565286121068044332677232658064393139116462431424107524782511263114636440912060806899003582551806244664571732272406734409742057323455810546875000000]
Assert failed. 0 does not evaluate to true
  ***   at top-level: ...mials-",if(manual_select,run_bsgs_experiment(s
  ***                                             ^---------------------
  ***   in function run_bsgs_experiment: ...getabstime();bsgs_output=
  ***   bsgs(K,cpct_units,sca
  ***   ^---------------------
  ***   in function bsgs: ...ge=="LOG",lattice_lambda=incremental_giant_ste
  ***                                                 ^---------------------
  ***   in function incremental_giant_steps: ...TE==0)&&ctr>0,base_value=
  ***   compact_rep_full_inpu
  ***   ^---------------------
  ***   in function compact_rep_full_input: ..." ",alphaOK,"\n",alpha));
  ***   GP_ASSERT_TRUE(idealB
  ***   ^---------------------
  ***   in function GP_ASSERT_TRUE: ...s not evaluate to true");breakpoint();)
  ***                                                           ^--------------
... skipping file 'experiment-babystock11.py'

  ***   Break loop: <Return> to continue; 'break' to go back to GP prompt
Goodbye!
*/
}
