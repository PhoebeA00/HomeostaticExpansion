     syms x y z w
     eq1 = x - 2*y + z + 2* w = 0;
     eq2 = 2*y - 8*z + w = 8;
     eq3 = -4*x + 5*y + 9*z - w = -9;
     [x,y,z] = solve(eq1, eq2, eq3, x,y,z)