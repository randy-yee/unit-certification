\\ test-vector_methods
read("src/notextendedlll.gp")
read("src/VectorMethods.py");
M = [1,-1,3;1,0,5;1,2,6];

M = [201, 1648;37,297];M = matconcat([M, [201,37]~]);

nrow = 5;
ncol = 6;
{
for(i = 1, 100,
    M = matrix(nrow,ncol, i, j, random(100));
    m1 = qflll(M);
    m2 = notextendedlll(M );

    m3 = my_mlll(M,10^(-10));

    \\print("qflll \n",M*m1);
    \\print("mlll \n",m2);
    \\print("my_mlll \n",m3);
    GP_ASSERT_EQ(abs(matdet(m2)), abs(matdet(m3)));
    \\print(abs(matdet(m2)) == abs(matdet(m3)));
);
print("Testing MLLL complete");

}

/* to do, add test(s) with real number inputs

L1 = [1.5859297750472545312900000000000000000000000000000, -27.414172950875466109285957156672501886085064126689;\
4.8353724194323870539200000000000000000000000000000, 6.9458669914539513610446200082691350486301820089770];

new_vec = [-1.7730975039190487418850343924828014226493954157399, -8.4445902564683182992020515929259551927226062356424];
my_mlll( matconcat([L1, new_vec~]),10^-10);

*/
