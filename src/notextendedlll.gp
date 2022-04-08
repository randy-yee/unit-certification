\\                                                                         \\
\\  Lenstra Lenstra Lovasz algorithm.         WvdK, July 1998              \\
\\                                                                         \\
\\  NAME notextendedlll                                                    \\
\\                                                                         \\
\\  DESCRIPTION                                                            \\
\\                                                                         \\
\\  By default notextendedlll computes an LLL reduced basis of the lattice \\
\\  L spanned by the generating set. It is not an extended algorithm in    \\
\\  that it does not compute the transformation which expresses the        \\
\\  answer in terms of the given generating set. This makes it faster.     \\
\\  The code is derived from the code for an extended algorithm. It        \\
\\  has been cleaned up further than the C version notextendedlll.c        \\
\\                                                                         \\
\\                                                                         \\
\\  SYNOPSIS                                                               \\
\\                                                                         \\
\\  notextendedlll(m[,y[,reduceoff]])                                      \\
\\                                                                         \\
\\  Here m contains the generators as columns,                             \\
\\  y is a rational number, 1/4<y<1 and 
\\  reduceoff is 0 or 1.                                                   \\
{

\\  PROGRAM  \\

notextendedlll(X,\
\\                              X is matrix with integer entries           \\

\\         OPTIONS                                                         \\

           y,\
\\                              y is the ratio that determines the quality \\
\\                              of the LLL reduction. Its default is 3/4.  \\
           reduceoff,\
\\                              option to indicate if nonreduced basis of  \\
\\                              the lattice L is acceptable. Default is 0. \\

\\         LOCAL VARIABLES                                                 \\

           B,\
\\            table of denominators                                        \\
           H,\
\\            matrix containing latest generators                          \\
           Lam,\
\\              table of numerators                                        \\
           R,
\\           table of volume ratios                                        \\

           A,C,D,Matr11,Matr12,Matr21,Matr22,Q,Q2,T1,T3,Tmp,Tmp2,V1,V3,\
           fr,isodim,k,kk,kmax,l,ll,n,nxt,notfinished,rnk,s,s1,thr)\

\\         LOCAL LOOP VARIABLES ii,j                                       \\

=\
\\  INITIALIZATION  \\

\\ Check matrix argument. \\
if(type(X)!="t_MAT",print("notextendedlll: argument is no matrix");1/0,);
Q2=matsize(X);
n=Q2[2];
for(ii=1,Q2[1],\
    for(j=1,n,\
        if(type(X[ii,j])!="t_INT",\
          print("notextendedlll: matrix entry is no integer");1/0\
          ,\
          )\
       )\
   );
\\ Get numerator and denominator of the LLL quality ratio y. \\
if(y<=1/4,y=3/4,);
if(y>=1,y=3/4,);
if(type(y)!="t_FRAC",print("notextendedlll: argument y is no rational number");1/0,);
thr=numerator(y);
fr=denominator(y);
\\ Switch. \\
s1=if((1-reduceoff),1,0);
notfinished=1;
H=X;
B=vector(n+1,ii,1);
Lam=matid(n);;
rnk=0;
isodim=0;
kmax=0;
k=1;
kmax=k;
\\ Gram Schmidt coefficients of first generator. \\
B[k+1]=(H[,k]~)*H[,k];
s=sign(B[k+1]);
if(s<=0,\
\\ It is a zero vector. The rank of the lattice R of relations increases. \\
  isodim=isodim+1;
  B[isodim+1]=1;
  ,\
\\ It is not zero, and rank of lattice L increases. \\
  rnk=rnk+1\
  );
nxt=0;
k=2;
if((k<=n),nxt=1,notfinished=0);
\\  MAIN LOOP  \\
while(notfinished,\
      if(nxt,\
        if(k>kmax,\
\\        Add generator.     \\
          kmax=k;
\\        Gram Schmidt coefficients of new generator.        \\
          for(j=isodim+1,k,\
              Q=(H[,j]~)*H[,k];
              for(ii=isodim+1,j-1,Q=(B[ii+1]*Q-Lam[k,ii]*Lam[j,ii])/B[ii]);
              if(j<k,Lam[k,j]=Q,B[k+1]=Q)\
             );
          s=sign(B[k+1]);
          if(s<=0,\
\\          New relation expected.     \\
            if(s,\
              print("notextendedlll: This should not happen.");1/0\
              ,\
              );
            R=vector(kmax,ii,1);
\\          Push generator down into subspaces.          \\
            for(ii=1,rnk,\
                k=k-1;
\\              Extended Euclid.              \\
                A=B[k+1];C=Lam[k+1,k];Matr21=1;D=A;
                if(C,\
                  V1=0;V3=C;
                  while(V3,\
                        Q2=divrem(D,V3);Q=Q2[1];T3=Q2[2];
                        T1=Matr21-Q*V1;Matr21=V1;D=V3;V1=T1;V3=T3\
                       );
                  Matr22=(D-A*Matr21)/C;
                  ,\
                  Matr22=0\
                 );
                Matr11=-C/D;Matr12=A/D;
\\              Matrix which pushes generator one down is            \\
\\              Matr=[Matr11,Matr12;Matr21,Matr22];                  \\
\\              To estimate its 11 entry, note                       \\
\\              Matr11 == - (C/A)*Matr12                             \\

\\              Volume ratio between old and new sublattice.              \\
                R[k+1]=Matr12;
\\              Push one down.              \\
                Tmp=Matr11*H[,k]+Matr12*H[,k+1];
                H[,k+1]=Matr21*H[,k]+Matr22*H[,k+1];H[,k]=Tmp;
\\              Update part of Lam.              \\
                for(j=1,k-1,\
                    Q=Matr11*Lam[k,j]+Matr12*Lam[k+1,j];
                    Lam[k+1,j]=Matr21*Lam[k,j]+Matr22*Lam[k+1,j];Lam[k,j]=Q\
                   );\
\\              Reduce generator k with generators k-1 down to 1.         \\
                    forstep(l=k-1,1,-1,\
                            A=B[l+1];
                            if(abs(2*Lam[k,l])>A,\
                              Q=((2*Lam[k,l]+A)\(2*A));
                              H[,k]=H[,k]-Q*H[,l];
                              Lam[k,l]=Lam[k,l]-Q*A;
                              for(j=1,l-1,Lam[k,j]=Lam[k,j]-Q*Lam[l,j])\
                              ,\
                              )\
                           );\
               );
\\          Add new relation.          \\
            isodim=isodim+1;
\\          Cumulative volume ratios.          \\
            for(ii=2,kmax,R[ii]=R[ii]*R[ii-1]);
\\          Only now do we update B and Lam.          \\
\\          Thus we avoided a lot of updating.       \\
            forstep(k=kmax,isodim+1,-1,\
                    B[k+1]=(B[k])/(R[k]^2);
\\                  These are integers, implying a bound on R[k].          \\
                    Q=R[k]*R[k-1];
                    for(j=k+1,kmax,Lam[j,k]=(Lam[j,k-1])/(Q));
                   );
            B[isodim+1]=1;
\\          Gram Schmidt coefficients of new zero vector to be made up.    \\
\\          Pretend there is an extra row in which this vector has entry 1 \\
\\          and other entries vanish.                                      \\
            for(k=1,isodim-1,Lam[isodim,k]=0);
            for(k=isodim+1,kmax,Lam[k,isodim]=0);
\\          Now reduce all the vectors.                                    \\
            for(k=isodim,kmax,\
\\          Reduce generator k with generators k-1 down to isodim+1.       \\
                forstep(l=k-1,isodim+1,-1,\
                        A=B[l+1];
                        if(abs(2*Lam[k,l])>A,\
                          Q=((2*Lam[k,l]+A)\(2*A));
                          H[,k]=H[,k]-Q*H[,l];
                          Lam[k,l]=Lam[k,l]-Q*A;
                          for(ii=1,l-1,Lam[k,ii]=Lam[k,ii]-Q*Lam[l,ii])\
                          ,\
                          )\
                       );\
                );
            k=max(isodim,2)\
            ,\
\\          No new relation, and rank of lattice L increases.         \\
            rnk=rnk+1\
            )\
          ,\
          );
        nxt=0\
        ,\
\\      Test for swap.      \\
\\      First reduce generator k with generator k-1.      \\
        l=(k-1);A=B[l+1];
        if(abs(2*Lam[k,l])>A,\
          Q=((2*Lam[k,l]+A)\(2*A));
          H[,k]=H[,k]-Q*H[,l];
          Lam[k,l]=Lam[k,l]-Q*A;
          for(ii=1,l-1,Lam[k,ii]=Lam[k,ii]-Q*Lam[l,ii])\
          ,\
          );
\\      Compute if condition for swap is satisfied.      \\
        cond=(k>isodim+1)&&\
             s1&&(fr*B[k+1]*B[k-1]<(thr*B[k]^2-fr*Lam[k,k-1]^2));
        if(cond,\
\\        Swap.        \\
          Tmp=H[,k-1];H[,k-1]=H[,k];H[,k]=Tmp;
          for(j=1,k-2,Q=Lam[k,j];Lam[k,j]=Lam[k-1,j];Lam[k-1,j]=Q);
          Q=Lam[k,k-1];
\\        More Gram Schmidt coefficients to be updated.          \\
          u=(B[k-1]*B[k+1]+Q^2)/B[k];
          for(ii=k+1,kmax,\
              A=Lam[ii,k];
              Lam[ii,k]=(B[k+1]*Lam[ii,k-1]-Q*A)/B[k];
              Lam[ii,k-1]=(u*A+Q*Lam[ii,k])/B[k+1]\
             );
          B[k]=u;
          k=max(2,k-1)\
          ,\
\\        No swap required. Reduce generator k with         \\
\\        generators k-2 down to isodim+1.                  \\
          forstep(l=k-2,isodim+1,-1,\
                  A=B[l+1];
                  if(abs(2*Lam[k,l])>A,\
                    Q=((2*Lam[k,l]+A)\(2*A));
                    H[,k]=H[,k]-Q*H[,l];
                    Lam[k,l]=Lam[k,l]-Q*A;
                    for(ii=1,l-1,Lam[k,ii]=Lam[k,ii]-Q*Lam[l,ii])\
                    ,\
                    )\
                 );
\\        Go up.        \\
          k=k+1;
          if(k<=n,nxt=1,notfinished=0)\
          )\
        )\
   );
\\  The later columns of H form a basis of L.                   \\
\\  This basis is reduced unless reduceoff is 1.                \\
    vecextract(H,\
               vector(matsize(H)[1],ii,ii),\
               vector(rnk,ii,n-rnk+ii)\
              )\
}

