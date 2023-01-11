function [Aineq, B, lb, ub] = bounds(Nnode, rmin, rmax, tmin)
% generates linear equality constraints and lower and upper bounds
% Inputs: 
% Nnode: number of nodes
% rmin: minimum inner radius
% rmax: Max outter radius
% tmin: min thickness

    Aineq = eye(Nnode, 2* Nnode);
    B = zeros(Nnode,1)+rmax;
    j = 1;
    for i = Nnode+1:2*Nnode
        % populates the second half of the equality constraint matrix
        % essentially creates an identity matrix for the secong half
        Aineq(j,i) = 1;
        j = j +1;
    end
    %lower bound on inner rad and thichness
        % 1cm = rmin m ; 2.5mm = tmin 
    %upperbound on inner rad and thichness
        % 5cm - 2.5mm; 5cm - 1cm
    lb1 = ones(Nnode,1)*rmin;
    lb2 = ones(Nnode,1)*tmin;
    lb = [lb1; lb2];
    ub1 = ones(Nnode,1)*(rmax - tmin);
    ub2 = ones(Nnode,1)*(rmax - rmin);
    ub = [ub1; ub2];
end