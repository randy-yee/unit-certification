{
 read("src/VectorMethods.py");
\\ If using Alltest.gp, the file reads below are not needed
\\ read("test/TestingUtility.gp");
}


{ \\ test cases for make_poly
    my(coeffs = [1,2,3,4]);
    GP_ASSERT_EQ(make_poly(coeffs), x^3 + 2*x^2 + 3*x + 4);

    coeffs = [9];
    GP_ASSERT_EQ(make_poly(coeffs),9);

}

{ \\ test cases for invert_coordinates
    my(vector1, vector2, eps,output_vector,output_vector2);

    vector1 = [1,2,3,4,5]~;
    output_vector = invert_coordinates(vector1);
    GP_ASSERT_EQ(output_vector, [1, 1/2, 1/3, 1/4, 1/5]~);


    eps = 10^(-20);
    vector2 = [1.1, 2.2, 4.5];
    output_vector2 = invert_coordinates(vector2);
    GP_ASSERT_NEAR(output_vector2[1], 1/(1.1),eps);
    GP_ASSERT_NEAR(output_vector2[2], 1/(2.2),eps);
    GP_ASSERT_NEAR(output_vector2[3], 1/(4.5),eps);

}

{ \\ test cases for mulvec
    my(v1,v2);

    M1 = Mat([1,2;3,4]);
    v1 = [10,5];
    GP_ASSERT_EQ(mulvec(M1,v1), Mat([10,20;15,20]) );

    eps = 10^(-20);
    M2 = Mat([1.1113,2.2222, 4.5678;1111000000,400.1212121212, 99999999]);
    v2 = [1.55,2.2222];
    Mv = mulvec(M2, v2);
    GP_ASSERT_NEAR(Mv[1,1], 1.1113*1.55, eps );
    GP_ASSERT_NEAR(Mv[1,2], 2.2222*1.55, eps );
    GP_ASSERT_NEAR(Mv[1,3], 4.5678*1.55, eps );
    GP_ASSERT_NEAR(Mv[2,1], 1111000000*2.2222, eps );
    GP_ASSERT_NEAR(Mv[2,2], 400.1212121212*2.2222, eps );
    GP_ASSERT_NEAR(Mv[2,3], 99999999*2.2222, eps );
}

{ \\ test cases for pointwise_vector_mul
    my(v1,v2);

    v1 = [1,2,3,4,5]~;
    v2 = [10,11,5,7,8]~;
    GP_ASSERT_EQ(pointwise_vector_mul(v1,v2), [10,22, 15, 28,40]~);

}


{ \\ test cases for dividevector

    my(v1,M1);

    M1 = Mat([1,2;3,4]);
    v1 = [10,11]~;
    GP_ASSERT_EQ(dividevector(M1,v1), Mat([1/10,2/10; 3/11, 4/11]));

}

{ \\ test cases for pointwise_vector_div

    my(v1,M1);

    M1 = Mat([1,2;3,4]);
    v1 = [10,11]~;
    GP_ASSERT_EQ(dividevector(M1,v1), Mat([1/10,2/10; 3/11, 4/11]));

}

{ \\ test cases for concat_negative_sum
    my(v);
    v=[1,2,3,4, 1.1, 6.7];
    GP_ASSERT_EQ(sumvec(v), 17.8);
}

{ \\ test cases for concat_negative_sum
    my(v1);
    v1 = [1,2,3,4,5,6,7,8,9];
    GP_ASSERT_EQ(concat_negative_sum(v1), [1,2,3,4,5,6,7,8,9,-45]);
}
{ \\ test cases for remove_last_coordinate
    my(v1);
    v1 = [1,2,3,4,5,6,7,8,9];
    GP_ASSERT_EQ(remove_last_coordinate(v1), [1,2,3,4,5,6,7,8]);
}

{ \\ test cases for expvec and inverse_expvec
    my(v1, v2, eps = 10^(-20));
    v1 = [1,2,3,4,5];
    GP_ASSERT_VEC_NEAR(expvec(v1), [exp(-1),exp(-2), exp(-3),
    exp(-4), exp(-5), exp(15)], eps );

    GP_ASSERT_VEC_NEAR(inverse_expvec(v1), [exp(1),exp(2), exp(3),
    exp(4), exp(5), exp(-15)], eps );
}

{ \\ test cases for vector_approximate
    my(v1, v2, eps = 10^(-9));

    GP_ASSERT_VEC_NEAR(, , eps);
}

{ \\ test cases for logvector
    my(v1, v2, eps = 10^(-9));
    v1 = [1.111111111, 2.121212121212, 3.3313131313131313];
    GP_ASSERT_VEC_NEAR(vector_approximate(v1, eps), [1.111111111,2.121212121,3.331313131 ], eps*10^(-10));
}
{ \\ test cases for is_trace_zero
    my(v1, v2,v3, eps = 10^(-9));
    v1 = [1,2,3,4, -10];
    v2 = [11.11, 12121.987243, 100];
    v3 = concat(v2, -sumvec(v2));
    GP_ASSERT_TRUE(is_trace_zero(v1, eps));
    GP_ASSERT_FALSE(is_trace_zero(v2, eps));
    GP_ASSERT_TRUE(is_trace_zero(v3, eps));
}

{ \\ test cases for samevecs
    my(v1, v2, eps = 10^(-9));
    v1= [1.1111111,2.23455678, 3.333555577];
    v2= [1.1111111003,2.23455678002, 3.333555577111];
    GP_ASSERT_TRUE(samevecs(v1,v2, eps));

    my(v1, v2, eps = 10^(-9));
    v1 = [1,1,1,2,2.234556780021, 3.3335555771111];
    v2= [0,0,0, -1.1111111003,2.23455678002, 3.333555577111];
    v3 = [0.1,0.1,0.1, -1.11111110029999,2.23455678002, 3.3335555771110001];

    GP_ASSERT_TRUE(samevecs(v2,v3,0.1 ) );
}

{ \\ test cases for vec_flip_positive

    my(v1, v2, v3, eps = 10^(-9));
    v1= [1.1111111,2.23455678, 3.333555577];
    v2= [0,0,0, -1.1111111003,2.23455678002, 3.333555577111];
    v3= [-0.00001111111,2.23455678, 3.333555577];
    GP_ASSERT_VEC_NEAR(vec_flip_positive(v1),v1, eps);
    GP_ASSERT_VEC_NEAR(vec_flip_positive(v2),-v2, eps);
    GP_ASSERT_VEC_NEAR(vec_flip_positive(v3),-v3, eps);
}

{ \\ test cases for vec_less_than

    my(v1, v2, eps = 10^(-9));
    v1 = [1,1,1,2,2.234556780021, 3.3335555771111];
    v2= [0,0,0, -1.1111111003,2.23455678002, 3.333555577111];
    v3 = [0.1,0.1,0.1, -1.11111110029999,-100, 3.3335555771110001];

    GP_ASSERT_TRUE(vec_less_than(v2,v1 ) );
    GP_ASSERT_FALSE(vec_less_than(v1,v2 ) );
    GP_ASSERT_TRUE(vec_less_than(v2,v3,5 ) );

}


print("Testing vector methods finished");
