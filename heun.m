function [t,y] = heun(t0, tN, y0, h, f)

    % input argument description:
        % t0 = left bound of interval on which to solve ODE
        % tN = right bound of interval on which to solve ODE
        % y0 = initial condition of ODE
        % h = step size for t
        % f = right hand side of the given ODE = what dy/dt is equal to
    
    % output argument description:
        % t = the array of independent variables
        % y = the array of dependent variables
    

    N = (tN-t0)/h;                      % N = number of steps = number of discrete t points
    t = linspace(t0, tN, N);            % t = array of independent variable t
    y = t.*0;                           % y = 1xN array of dependent variable y
    y(1) = y0;                          % setting the 1st element of y to be the initial value

    for i=1:(N-1)
        s1 = f(t(i), y(i));                     % s1 = slope at i
        y(i+1) = y(i) + (s1*h);                 % calculating y(i+1) using s1
        s2 = f(t(i+1), y(i+1));                 % s2 = slope at i+1
        y(i+1) = y(i) + (0.5*(s1+s2)*h);        % recalculating y(i+1) using s1 and s2
    end

end

