
\\ mlll implementation as described in Cohen I

DEBUG_MLLL = 0;


\\ mat is the original basis
\\ H is the transformation matrix
\\ mu_k is a vector consisting of the mu_{k,i} values for all i
\\ l indicates the position of interest
\\ See page 88 of Cohen I
cohen_subalg_red(mat, Hmat, mu_mat,k,l)={
    my(scale);
    if(abs(mu_mat[k,l])<= 0.5,
        return([mat, Hmat, mu_mat]);

    , \\else
        scale = round(mu_mat[k,l]);

        mat[,k] = mat[,k] - scale*mat[,l];
        Hmat[,k] = Hmat[,k] -scale*Hmat[,l];
        mu_mat[k,l] = mu_mat[k,l] - scale;

        for(i = 1, l-1,
            mu_mat[k,i] = mu_mat[k,i] - scale*mu_mat[l,i];
        );

        return([mat, Hmat, mu_mat]);
    )

}
\\ Swap and GS subroutine
cohen_subalg_swapg(mat, new_mat, Hmat, mu_mat, B_vec, k, kmax,eps)={

    \\print("SWAP: k=",k);
    my(vec_holder,mu_holder = 0);

    \\startdata = [mat, new_mat, Hmat, mu_mat,B_vec];

    vec_holder = mat[,k];
    mat[,k] = mat[,k-1]; mat[,k-1] = vec_holder;

    vec_holder = Hmat[,k];
    Hmat[,k] = Hmat[,k-1]; Hmat[,k-1] = vec_holder;

    if(k >2,
        for(j=1, k-2,
            mu_holder = mu_mat[k,j];
            mu_mat[k,j] = mu_mat[k-1,j];
            mu_mat[k-1,j] = mu_holder;
        );
    );
    mu = mu_mat[k,k-1];
    B = B_vec[k]+(mu^2)*B_vec[k-1];

    \\print("SWAPG ", precision(B,15), "  ", precision(B_vec[k],10),"  ",precision(abs(mu), 10));
    if(abs(B) < eps,
        \\print(" - if 1");
        mu_holder = B_vec[k];
        B_vec[k] = B_vec[k-1];
        B_vec[k-1] = mu_holder;

        vec_holder= new_mat[,k];
        new_mat[,k] = new_mat[,k-1];
        new_mat[,k-1] = vec_holder;
        for(i=k+1, kmax,
            mu_holder = mu_mat[i,k];
            mu_mat[i,k]= mu_mat[i,k-1];
            mu_mat[i,k-1] = mu_holder;
        );
    ,
    \\elif
    abs(B_vec[k]) < eps && abs(mu) > eps,
        \\print(" - if 2");
        B_vec[k-1] = B;
        new_mat[,k-1] = new_mat[,k-1]*mu;
        mu_mat[k,k-1] = 1/mu;
        for(i=k+1, kmax,
            mu_mat[i,k-1]=mu_mat[i,k-1]/mu;
        );
    ,
    \\elif
    abs(B_vec[k]) > eps,
        \\print(" - if 3");
        t= B_vec[k-1]/B;
        mu_mat[k,k-1] = mu*t;
        vec_holder = new_mat[,k-1];
        new_mat[,k-1] = new_mat[,k]+mu*vec_holder;
        new_mat[,k] = -mu_mat[k,k-1]*new_mat[,k] + (B_vec[k]/B)*vec_holder;
        B_vec[k] = B_vec[k]*t;
        B_vec[k-1] = B;

        for(i=k+1, kmax,
            t = mu_mat[i,k];
            mu_mat[i,k] = mu_mat[i,k-1]-mu*t;
            mu_mat[i,k-1] = t + mu_mat[k,k-1]*mu_mat[i,k];
        )
    );

    return([mat, new_mat, Hmat, mu_mat, B_vec]);

}

\\ Attempt at implementing MLLL (for real matrices with dependent columns)
my_mlll(mat,eps)={
    if(DEBUG_MLLL, print("Start"));
    my(k, kmax, GS_mat, B_vec,Hmat, mu_mat, new_b, zeroctr = 0);

    \\1. initialize
    k = 2;
    kmax = 1;
    GS_mat = matrix(matsize(mat)[1],matsize(mat)[2], i, j, 0);
    GS_mat[,1] = mat[,1];
    B_vec = vector(length(mat), i,0); B_vec[1] = mat[,1]~*mat[,1];
    Hmat = matid(length(mat));
    mu_mat = Hmat;

    while(k <= length(mat),
        if(DEBUG_MLLL,
            print("While loop top: k = ",k);
            print(precision(mat,10), "  \n", precision(GS_mat,10), "\n", precision(B_vec,10), "\n");
        );
        \\ incremental G-S
        if(k <= kmax,
        ,                                                       \\ comma indicates do nothing
            kmax = k;                                           \\ different from Cohen kmax = k
            new_b = mat[,k];
            for(j=1, k-1,
                if(abs(B_vec[j]) < eps,
                    mu_mat[k,j] = 0;
                , \\else
                    mu_mat[k,j] = (mat[,k]~*GS_mat[,j])/B_vec[j];
                );
                new_b = new_b - mu_mat[k,j]*GS_mat[,j];

            );
            \\print("orthocheck ", precision(new_b~*mat[,1],10));
            if (norml2(new_b) < eps,
                B_vec[k] = 0;
                GS_mat[,k] = vector(length(new_b), t, 0)~;
            ,
                B_vec[k] = new_b~*new_b;
                GS_mat[,k] = new_b;
            );

        );

        if(DEBUG_MLLL, print("\nafter G-S: " precision(GS_mat,10), "\n ", precision(B_vec,10)););
        \\\\\\\\\\\\\\\\\\\\
        \\ [Test LLL condition]
        flag = 1;
        while(flag,
            [mat, Hmat , mu_mat] = cohen_subalg_red(mat, Hmat, mu_mat, k, k-1);

            if(DEBUG_MLLL, print("LOVASZ: ", precision(B_vec[k],10), "   ", precision((0.75 - mu_mat[k,k-1]^2)*B_vec[k-1],10)));
            if(B_vec[k] < (0.75 - mu_mat[k,k-1]^2)*B_vec[k-1],
                if(DEBUG_MLLL,print("Lovasz condition fails"););

                [mat, GS_mat, Hmat, mu_mat, B_vec] = cohen_subalg_swapg(mat, GS_mat, Hmat, mu_mat, B_vec, k, kmax, eps);
                k = max(2, k-1);
                \\print("after LLL condition: " GS_mat, "\n ", B_vec);
            , \\else
                if(DEBUG_MLLL,print("Lovasz condition met"););
                for(l =2, k-1,
                    [mat, Hmat , mu_mat] = cohen_subalg_red(mat, Hmat, mu_mat, k, k-l);
                );
                k = k+1;
                flag = 0;
            );
        );
        if(DEBUG_MLLL, print("k incremented, returning to GS"));

    );


    \\ 4. [Finished?]
    for(i =1, length(GS_mat),
        if(DEBUG_MLLL, print("End Norms: ", precision(norml2(mat[,i]),10)));
        if(norml2(mat[,i]) < eps, zeroctr += 1;);
    );
    return(mat[,zeroctr+1..length(mat)]);
}
