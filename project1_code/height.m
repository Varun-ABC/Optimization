function [h] = height(a, L, x)
    % height calculates the y value at each x value using a sine series
    % parametrization
    % a is the matrix of coefficents 
    % L is the length of the section
    % x is the matrix of evenly spaced points along the x-axis
    % h (output) is a matrix of the y value at each x value
    nvar = numel(a);
    h = zeros(size(x,1),1); %initilizes a matrix of zeros that is the same size as x
    h = h+ a(1); % in the parametrization used, the first values, a1, is added to all elements 
    for k = 2:nvar 
        % loops over the number of design varibales
        for i = 1:size(x,1)
            %loops over the number of x elements
            h(i) = h(i) + a(k)* sin(2*pi*(k-1)*x(i)/L); 
            % the height of each element is iterativly added to get the
            % sine series
        end
    end