function [Aineq, bineq] = limits(a, Nx, L, x)
    %creates linear constraints for the upper and lower bounds
    % a is the matrix of coefficents
    % Nx is the amount of subdivisions in the Xdirection
    % L is the length of the section
    % x is the matrix of evenly spaced points along the x-axis
    % Ax < b
    % Aineq is A in the inequality
    % bineq is b in the inequality
    
    nvar = numel(a); % number of design variables
    Aineq = zeros(2*(Nx+1),nvar);
    bineq = zeros(2*(Nx+1),1);
    for i = 1:(Nx+1)
      % first, the upper bound
      Aineq(i,1) = 1.0; % this coefficient corresponds to variable a_1
      for k = 2:nvar
        Aineq(i,k) = sin(2*pi*(k-1)*x(i)/L); % this coefficient corresponds to variable a_k
      end
      bineq(i) = 0.05; % the upper bound value
      % next, the lower bound; we use ptr to shift the index in Aineq and bineq
      ptr = Nx+1;
      Aineq(ptr+i,1) = -1.0; % note the negative here!!! fmincon expects inequality in a form A*x < b
      for k = 2:nvar
        Aineq(ptr+i,k) = -sin(2*pi*(k-1)*x(i)/L); % again, a negative
      end
      bineq(ptr+i) = -0.01; % negative lower bound
    end