{
read("src/Neighbours.py");
\\ If using Alltest.gp, the file reads below are not needed, but the tests
\\ depend on these files
\\ read("src/VectorMethods.py");
\\ read("src/Neighbours.py");

}

{ \\ testing cubescan_equal, cubescan_unequal, is_minimum
    default(realprecision, 100);
    G = nfinit(x^4 - 2*x^3 + 2*x^2 + 2*x - 1); eps = 10^(-10);
    G1 =  bnfinit(x^4 - 2*x^3 + 2*x^2 + 2*x - 1);
    n = poldegree(G.pol);
    integral_basis_real = embed_real(G,G[5][1]);
    gram_zk = integral_basis_real~*integral_basis_real;
    short_list = qfminim(gram_zk, 500, , 2)[3];
    unit1 = nfalgtobasis(G,G1.fu[1]);
    testunit = unit1;
    ring_of_integers = matid(n);

    GP_ASSERT_TRUE(length(cubescan_equal(ring_of_integers, valuationvec(G,nfalgtobasis(G,1)),G, eps )) ==1);

    for(i=1, 10,

        GP_ASSERT_TRUE(cubescan_unequal(ring_of_integers, valuationvec(G,testunit),G, eps ) ==[]);
        GP_ASSERT_TRUE(is_minimum(ring_of_integers, testunit,G, eps ));
        testunit = nfeltmul(G, testunit, unit1);

    );
}

{ \\ testing NEIGHBORS
    default(realprecision, 500);
    G = nfinit(x^3 - x^2 - 10*x + 13); eps = 10^(-20);
    G1 =  bnfinit(x^3 - x^2 - 10*x + 13);
    n = poldegree(G.pol);
    unit1 = nfalgtobasis(G,G1.fu[1]);
    testunit = unit1;
    ring_of_integers = matid(n);

    GP_ASSERT_TRUE(length(cubescan_equal(ring_of_integers, valuationvec(G,nfalgtobasis(G,1)),G, eps )) ==1);
    GP_ASSERT_TRUE(cubescan_unequal(ring_of_integers, valuationvec(G,testunit),G, eps ) ==[]);
    GP_ASSERT_TRUE(is_minimum(ring_of_integers, testunit,G, eps ));

    one_neighbors = compute_one_neighbours(G, ring_of_integers, eps);
    for(i=1, length(one_neighbors),
        GP_ASSERT_TRUE(is_minimum(ring_of_integers, one_neighbors[i],G, eps ) );
        GP_ASSERT_TRUE(is_neighbour(ring_of_integers, nfalgtobasis(G,1), one_neighbors[i], G, eps ));
        nextneighbours = NEIGHBORS(G, ring_of_integers, one_neighbors[i], eps);
        for(j=1, length(nextneighbours),
            GP_ASSERT_TRUE(is_minimum(ring_of_integers, one_neighbors[j],G, eps););
        );
    );
}

{ \\ testing COLLECT
    default(realprecision, 100);
    G = nfinit(x^3 - x^2 - 10*x + 13); eps = 10^(-20);
    G1 =  bnfinit(x^3 - x^2 - 10*x + 13);
    n = poldegree(G.pol);
    integral_basis_real = embed_real(G,G[5][1]);
    gram_zk = integral_basis_real~*integral_basis_real;
    short_list = qfminim(gram_zk, 90, , 2)[3];
    unit1 = nfalgtobasis(G,G1.fu[1]);
    testunit = unit1;
    ring_of_integers = matid(n);
    boundary = [10,10,10];
    my_collection = COLLECT(G, ring_of_integers, boundary, eps);
    one_neighbors = compute_one_neighbours(G, ring_of_integers, eps);

    \\for(i=1, length(my_collection), print(norml2(valuationvec(G, my_collection[i]))));
    for(i=1, length(short_list),
        if(is_minimum(ring_of_integers,short_list[,i],G,eps ),
            GP_ASSERT_TRUE(setsearch(my_collection, vec_flip_positive(short_list[,i]))!=0 );
        );
    );
}
print("Testing neighbours complete");
