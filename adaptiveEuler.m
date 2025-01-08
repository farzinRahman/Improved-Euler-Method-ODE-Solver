function [t,y] = adaptiveEuler(t0, tN, y0, h, f)

    % input argument description:
        % t0 = left bound of interval on which to solve ODE
        % tN = right bound of interval on which to solve ODE
        % y0 = initial y value
        % h = initial step size for t
        % f = right hand side of the given ODE = what dy/dt is equal to
    
    % output argument description:
        % t = the array of independent variables
        % y = the array of dependent variables
    
    % initial conditions
    t = t0;
    y = y0;

    % declaring the error tolerance
    tol = 1e-8;

    % starting the while loop
    i=1;
    while t(size(t)) < tN
        % (a) 1st estimate of y(i+1) = Y (one Euler step of size h)
        s_pre = f(t(i), y(i));                      % s_pre = slope at i
        Y = y(i) + (s_pre*h);                       % Y = 1st estimate using regular step size
        s_post = f(t(i) + h, Y);                    % s_post = slope at i+1
        Y = y(i) + (0.5*(s_pre+s_post)*h);          % recalculating Y using s_pre and s_post
        
        % (b) 2nd estimate of y(i+1) = Z (two succesive euler step of size h/2)
        Z_mid = y(i) + (s_pre * (h/2));                 % Z_mid = estimate at midpoint using half step size
        s_mid = f(t(i) + (h/2), Z_mid);                 % s_mid = slope at midpoint
        Z_mid = y(i) + (0.5* (s_pre + s_mid)*(h/2));    % recalculating Z_mid using s_pre and s_mid

        Z = Z_mid + (s_mid*(h/2));                      % Z = estimate of y(i+1) using Z_mid
        s_post = f(t(i) + h, Z);                        % s_post = slope at i+1
        Z = Z_mid + (0.5*(s_mid + s_post)*(h/2));       % recalculating Z using s_mid and s_post

        D = abs(Z-Y);                       % calculating estimate of the error

        if D < tol
            y_next = Z+D;                   % accepting the calculated values
            t_next = t(length(t)) + h;
            
            y = [y y_next];                 % concatenating the new values to the arrays
            t = [t t_next];

            i = i+1;                        % advancing counter i
        end
        % (c) updating the step size h
        h = 0.9*h*min(max(tol/abs(D),0.3),2);
    end
end

