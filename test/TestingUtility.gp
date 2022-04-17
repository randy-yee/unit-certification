/*
Created on Mar 24, 2022
Author: Randy Yee


This file contains some helpful functions for testing to emulate
assert statements.
*/

GP_ASSERT_TRUE(arg1)={
  if(arg1,
  ,
  print("Assert failed. ", arg1, " does not evaluate to true");
  breakpoint();
  );
}
GP_ASSERT_FALSE(arg1)={
  if(!arg1,
  ,
  print("Assert failed. ", arg1, " does not evaluate to false");
  breakpoint();
  );
}
GP_ASSERT_EQ(arg1, arg2)={
  if(arg1 == arg2,
  ,
  print("Assert failed. ", arg1, " is not equal to ", arg2);
  breakpoint();
  );
}
GP_ASSERT_NEQ(arg1, arg2)={
  if(arg1 != arg2,
  ,
  print("Assert failed. ", arg1, " is equal to ", arg2);
  breakpoint();
  );
}

GP_ASSERT_NEAR(arg1, arg2, eps)={
  if(abs(arg1 - arg2) < eps,
  ,
    print("Assert failed. ", arg1, " is not within ", eps, " range of ", arg2);
    breakpoint();
  );
}

GP_ASSERT_VEC_NEAR(arg1, arg2, eps)={
  if(length(arg1) != length(arg2), print("Invalid Assert,
    vectors not the same length");breakpoint(); );
  for(i=1, length(arg1),
    if(abs(arg1[i] - arg2[i]) < eps,
    ,
      print("Assert failed. ", arg1[i], " is not within ", eps, " range of ", arg2[i]);
      breakpoint();
    );
  );
}

GP_ASSERT_MAT_NEAR(arg1, arg2, eps)={
  if(matsize(arg1) != matsize(arg2),
    print("Invalid Assert,
    matrices not the same dimension");breakpoint();
  );
  my(n_rows = matsize(arg1)[1], n_cols= matsize(arg1)[2]);
  for(i=1, n_rows,
    for(j=1, n_cols,
      if(abs(arg1[i,j] - arg2[i,j]) < eps,
      ,
      print("Assert failed. ", arg1[i], " is not within ", eps, " range of ", arg2[i]);
      breakpoint();
      );
    );
  );
}

GP_ASSERT_GT(arg1, arg2)={
  if( arg1 >arg2,
  ,
    print("Assert failed. ", arg1, " is not greater than ", arg2);
    breakpoint();
  );
}

GP_ASSERT_LT(arg1, arg2)={
  if( arg1 <arg2,
  ,
    print("Assert failed. ", arg1, " is not less than ", arg2);
    breakpoint();
  );
}

GP_ASSERT_VEC_LT(arg1, arg2,eps)={
  if(length(arg1) != length(arg2), print("Invalid Assert,
    vectors not the same length");breakpoint(); );
  for(i=1, length(arg1),
    if((arg1[i] - arg2[i]) < eps,
    ,
      print("Assert failed. ", arg1[i], " is not less than ", arg2[i]);
      breakpoint();
    );
  );
}


\\ Shortcut function way to print values at a particular precision.
debug_print(string, value, p_level = 10)={
    print(string, " ", precision(value,10));
}

\\ shortcut function for debugging. Prints two values at a given precision and breakpoints
debug_compare(value1, value2, p_level = 10)={
    debug_print("A ", value1, p_level);
    debug_print("B ", value2, p_level);
    breakpoint();
}
