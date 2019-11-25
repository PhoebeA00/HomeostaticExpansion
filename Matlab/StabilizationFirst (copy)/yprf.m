function yprf =yprf(t,y)
    a = 0.4;
    b = 0.1; 
    d = 0.01;
    e = 0.01;
    f = 1; 
    
    %T(t) = y(1); 
    %I(t) = y(2);
    
    %dT/dt
    yprf(1) = a*y(1)*y(2) - b*y(1);
    
    %dI/dt 
    yprf(2) = d*y(1) - e*y(1)*y(2) - f*y(2);
    
    yprf = [yprf(1) yprf(2)]';