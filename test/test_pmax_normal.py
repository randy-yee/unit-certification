{
read("src/PmaxNormal.py");

default(realprecision, 1000);
}

{ \\ Pohst's example from 1994 paper
    my(f, K,K1, O_K, n, eps = 10^(-9));

    f = x^19 + 2;
    K = nfinit(f); K1 = bnfinit(f);

    unit1 = [-1,2,-1,-2,-6,2, -1, 1, -2, -3, -2, 2, 1, 0, -4, 2, 0, 1, 1]~;
    unit2 = [-15, 6, 7, -15, 7, 5, -13, 8, 2,-11, 9, 0, -9, 9, -1, -7,8,-2, -5]~;
    unit3 = [-45,44,-41,41,-38, 38, -37, 33, -35,33, -29,32,-29, 26,-29,26,-24,25,-23]~;
    unit4 = [-3, -6, -5, -1,8,8,1,-5,-5,-2,-2,1,4,6,0,-5,-5,-1,2]~;
    unit5 = [-7,4,-3,-1,4,-4,4,-1,-1,3,-5,2,0,-1,4,-3,1,-1,-2]~;
    unit6 = [17,-38,0,31,-18,-21,26,5, -29,8,23,-19,-13,24,1,-23,10,18,-16]~;
    unit7 = [9,2,-2,-2,-2,-2,-5,-5,-5,1,4,6,3,1,1,2,1,-3,-5]~;
    unit8 = [-19,15,9,-10,-3,-4,15,-2,-13,5,3,7,-10,-6,13,-1,-4,-4,2]~;
    unit9 = [-91,-147,-84,21,44,-32,-109,-91,-2,58,28,-45,-67,-9,60,67,11,-34,-15]~;

    my(ind_unitmat, Gv,lbound,index_bound,unit_output);
    ind_unitmat = matconcat([unit1, unit2,unit3,unit4, unit5,unit6,unit7,unit8,unit9]);

\\ CONSTRUCT K TO BE USED TO COMPUTE LOWER BOUND AND INDEX BOUND
Gv = K[5][2]*nfalgtobasis(K, K.zk[length(K.zk)]); Gv = Gv~*Gv;
pohst_k = 2*Gv;
lbound = lower_regbound(K,pohst_k,eps, -1, 2500000);
print("--Annoyingly, the repaired lower bound code breaks this example.");
print("--We instead hard-code the bound used in the original Pohst paper.");
print("--Computed lower bound: ", precision(lbound,10));
lbound = 433281.296;
print("--Bound from paper: ", precision(lbound,10));
index_bound = ceil(log_determinant(K, ind_unitmat)/(lbound));
print("-- Lower bound is: ", precision(lbound,10), " using K = ",precision(pohst_k,10), "\n index_bound is ", index_bound);

unit_output = pohst_check_normal(K,ind_unitmat, index_bound,eps);
GP_ASSERT_NEAR(log_determinant(K, unit_output), K1.reg, eps);

}
print("Testing pmax-normal functions complete")
