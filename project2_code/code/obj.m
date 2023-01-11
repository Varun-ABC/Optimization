function [m, dm] = obj(rad_thic, Nnode , L, rho)
% finds the mass of the system, is our objective function that we want to
% minimize
% rad_thic: Column vector of inner radii followed by thickness at each node
% Nnode: Number of nodes
% L: Length of spar
% rho: Density of material spar is made of
% Outputs: 
% m: mass of spar
% dm: the gradient of the mass of the spar

    m = subobj(rad_thic); % finds mass of spar at
    dm = zeros(2 * Nnode,1);
    h = 1e-60;
    for j = 1:2 * Nnode
        % populates the gradient matrix using complex step
        xc = rad_thic; 
        xc(j) = xc(j) + complex(0,h);
        dm(j) = imag(subobj(xc))/h;
    end
    function [mass] = subobj(RT)
        % Finds the mass of the spar
        % input: 
        % RT: Column vector of inner radii followed by thickness at each node
        % Output:
        % mass: Mass of spar
        
        area = zeros(Nnode ,1);
        s = L/ (Nnode-1); % Spacing of the nodes 
        for i = 1:Nnode
            % finds the area of the anulus at each node
            ri = RT(i);
            ro = ri + RT(Nnode + i);
            area(i) = pi * (ro^2 - ri^2);
        end
        % volume can be found by taking the inetgral over each annulus and
        % multiplying by the spacing of each node
        % mass can be found by multiplying the volume by density
        mass = rho * s * trapz(area); 
    end
end