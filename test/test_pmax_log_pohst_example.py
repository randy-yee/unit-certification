{
read("src/PmaxLog.py");

default(realprecision, 1000);
}

{
    my(K1, K2, O_K, n, r, cpct_units, delta_K,B,
        output_regulator,pmax_log_output,
        lglat, eps = 10^(-20),
        f, Gv, pohst_k
    );
    \\Pohst example (degree 19)
    f = x^19 + 2;
    K1 = nfinit(f); K2 = bnfinit(f);
    unit1 = [-1,2,-1,-2,-6,2, -1, 1, -2, -3, -2, 2, 1, 0, -4, 2, 0, 1, 1]~;
    unit2 = [-15, 6, 7, -15, 7, 5, -13, 8, 2,-11, 9, 0, -9, 9, -1, -7,8,-2, -5]~;
    unit3 = [-45,44,-41,41,-38, 38, -37, 33, -35,33, -29,32,-29, 26,-29,26,-24,25,-23]~;
    unit4 = [-3, -6, -5, -1,8,8,1,-5,-5,-2,-2,1,4,6,0,-5,-5,-1,2]~;
    unit5 = [-7,4,-3,-1,4,-4,4,-1,-1,3,-5,2,0,-1,4,-3,1,-1,-2]~;
    unit6 = [17,-38,0,31,-18,-21,26,5, -29,8,23,-19,-13,24,1,-23,10,18,-16]~;
    unit7 = [9,2,-2,-2,-2,-2,-5,-5,-5,1,4,6,3,1,1,2,1,-3,-5]~;
    unit8 = [-19,15,9,-10,-3,-4,15,-2,-13,5,3,7,-10,-6,13,-1,-4,-4,2]~;
    unit9 = [-91,-147,-84,21,44,-32,-109,-91,-2,58,28,-45,-67,-9,60,67,11,-34,-15]~;
    ind_unitmat = matconcat([unit1, unit2,unit3,unit4, unit5,unit6,unit7,unit8,unit9]);
    \\ CONSTRUCT K TO BE USED TO COMPUTE LOWER BOUND AND INDEX BOUND
    Gv = K1[5][2]*nfalgtobasis(K1, K1.zk[length(K1.zk)]); Gv = Gv~*Gv;
    \\print("Gv = ", precision(Gv,10));
    pohst_k = 2*Gv;

    lbound = lower_regbound(K1,  pohst_k,eps);
    index_bound = ceil(log_determinant(K1, ind_unitmat)/(lbound));
    \\print("pohst - k is ", precision(pohst_k,10));
    \\print("Lower bound is: ", precision(lbound,10), " index_bound is ", index_bound);

    unit_embedding = K1[5][1]*ind_unitmat;
    lglat = matrix( length(unit_embedding), length(unit_embedding), i,j, unit_embedding[i,j] ) ;

    lglat = log(abs(lglat));

    GP_ASSERT_NEAR(unscaled_determinant(K1,lglat), 4*27*7*K2.reg,eps);
    sumv = lglat[,1];
    for(j=2, length(lglat), sumv+=lglat[,j]);
    X1 = prec_baby(poldegree(K1.pol), log(abs(K1.disc)), infinity_norm(sumv));
    X2 = prec_giant(poldegree(K1.pol), log(abs(K1.disc)),abs(reg1),infinity_norm(sumv) );
    \\print(ceil(X1), "   ", ceil(X2), "   ", max(ceil(X1),ceil(X2)));
    REQ_BSGS = ceil(max(ceil(X1),ceil(X2)));
    default(realprecision, ceil(REQ_BSGS));
    REQ_COMPARE = ceil((poldegree(K1.pol)^2 +2)*log(infinity_norm(sumv))+2*poldegree(K1.pol)^2 +5);
    eps = 2^(-REQ_COMPARE);

    cpct_units = cpct_from_loglattice(K1, lglat, eps);
    \\print(precision(lglat,10)); breakpoint();
    index_bound = 83718;
    pmax_log_output = log_pohst_pari(K1,lglat,cpct_units, index_bound, eps);
    output_regulator = unscaled_determinant(K1, pmax_log_output);
    GP_ASSERT_NEAR(K2.reg,output_regulator, eps );

}
print("pmax log for pohst example finished");
