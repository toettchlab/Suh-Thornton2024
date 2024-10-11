function model_fcn = construct_2D_model_flux_rescaled(N,L,pattern)
% function model_fcn = construct_2D_model_flux_rescaled(N,L,pattern)
% 
% This function builds a differential equation model of a 2D layer of
% cells that expands by diffusion and also experiences a flux at boundaries
% of a light pattern. 
%
% Input arguments:
% N       - the number of compartments to discretize each dimension
% L       - the length and width of the square of cells
% pattern - the illumination pattern driving localized cell movement
% 
% Output
% model_fcn - a function handle to an ODE system that can be run using
%             an ODE solver like ode15s().

% Get size of each box for scaling diffusion/flux in the rate equations
dxy = L/N; % size of each compartment

% Get diffusion parameters
D = 50; % in um^2/min

% Parameter associated with flux:
f = 0.3; % in min-1

% This is a set of coordinates representing every compartment in the 2D
% field of cells. For every index in "center", the index in "right" is the
% compartment to its right, "left" is the compartment to its left, etc.

% Assume matrix positions:
% y(center,1) = first COLUMN
% y(1,center) = first ROW

center = 1:N; % simulating N compartments

right = circshift(center, -1); % position 2 is RIGHT compared to 1
left  = circshift(center, 1);  % position 1 is LEFT compared to 2
down  = circshift(center, -1); % position 2 is DOWN compared to 1
up    = circshift(center, 1);  % 1 is UP compared to 2

% If there is no light pattern, initialize with a pattern of zeros.
if nargin < 3 || isempty(pattern)
    pattern = zeros(N,N);
end

% Get pattern indices
[idf, idt, ...
 iuf, iut, ...
 ilf, ilt, ...
 irf, irt] = get_pattern_inds(pattern);

model_fcn = @ode_system; % return function handle to model


    % THIS IS THE ODE SYSTEM RETURNED BY CONSTRUCT_2D_MODEL
    function dydt = ode_system(t,y)
    
        % reshape the ODE system to a square pattern of cells!
        y = reshape(y,N,N);
        
        % Assume matrix positions:
        % y(center,1) = first COLUMN
        % y(1,center) = first ROW
        % diffusion terms:
        
        % This equation defines diffusion within the 2D geometry
        dydt(center,center) = y(left,center) + y(right,center) + ... % IN from left & right
                              y(center,up) + y(center,down) - ...    % IN from up & down
                              4*y(center,center);                    % OUT from center
        
        % These equations define diffusion at the edges
        dydt(1,center)      = y(2,center) + ...                 % IN from down
                              y(1,right) + y(1,left) - ...      % IN from right & left
                              3*y(1,center);                    % OUT from center
        dydt(end,center)    = y(end-1,center) + ...             % IN from up
                              y(end,right) + y(end,left) - ...  % IN from right & left
                              3*y(end,center);                  % OUT from center
        dydt(center,1)      = y(up,1) + y(down,1) + ...         % IN from up & down
                              y(center,2) - ...                 % IN from right
                              3*y(center,1);                    % OUT from center
        dydt(center,end)    = y(up,end) + y(down,end) + ...     % IN from up & down
                              y(center,end-1) - ...             % IN from left
                              3*y(center,end);                  % OUT from center
        
        % These equations define diffusion at corner squares
        dydt(1,1)           = y(1,2)+y(2,1)-2*y(1,1);
        dydt(end,end)       = y(end,end-1)+y(end-1,end)-2*y(end,end);
        dydt(1,end)         = y(1,end-1)+y(2,end)-2*y(1,end);
        dydt(end,1)         = y(end-1,1)+y(end,2)-2*y(end,1);
        
        % All diffusion terms must be scaled by the box size "dxy" as below:
        dydt = D/dxy^2*dydt;
        
        % Now we further add cell movement based on light boundary fluxes. 
        % Flux upwards:
        dydt(iuf) = dydt(iuf) - f/dxy*y(iuf); % flux up FROM thes indices
        dydt(iut) = dydt(iut) + f/dxy*y(iuf); % flux up TO these indices
        % Flux downwards:
        dydt(idf) = dydt(idf) - f/dxy*y(idf); % flux up FROM thes indices
        dydt(idt) = dydt(idt) + f/dxy*y(idf); % flux up TO these indices
        % Flux leftwards:
        dydt(ilf) = dydt(ilf) - f/dxy*y(ilf); % flux up FROM thes indices
        dydt(ilt) = dydt(ilt) + f/dxy*y(ilf); % flux up TO these indices
        % Flux rightwards:
        dydt(irf) = dydt(irf) - f/dxy*y(irf); % flux up FROM thes indices
        dydt(irt) = dydt(irt) + f/dxy*y(irf); % flux up TO these indices
        
        % At the end, we re-pack the equations into a column vector!
        dydt = reshape(dydt,N^2,1);
    end
end
