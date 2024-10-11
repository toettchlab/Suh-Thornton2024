%% Initialize a clean workspace
clear all, close all

%% PART 1: SIMULATE THE PATTERN
load illuminati_pattern 

% this is just some simple manipulations to make sure the pattern is
% well-behaved (e.g., it is binary, and the boundaries are dark).
pattern = flipud(pattern);
pattern(:,1:2) = 0;
pattern(:,end-1:end) = 0;
pattern(1:2,:) = 0;
pattern(end-1:end,:) = 0;
pattern = pattern > 0;

% Simulate the expanding tissue under illumination
N = 101;
L = 5400;  % simulation width in um (5.4 mm = 5,400 um)
tF = 2880; % simulation time in min (48 h = 2,880 min)

f = construct_2D_model(N,L,pattern);

ic = reshape(ic,N^2,1);
tic
[t,y] = ode15s(f, [0 tF], ic);
toc

tv = linspace(0,tF,290);
yv = interp1(t, y, tv);

save Fig6_simulation_results L N tF tv yv

%% PART 2: Plot the results!
% This code plots the model simulation of Figure 6.
load Fig6_simulation_results 

xv = linspace(0,L,N);
[X Y] = meshgrid(xv, xv);
for i = 1:length(tv)
    Z = reshape(yv(i,:),N,N);
    h = pcolor(X,Y,Z);
    set(h, 'edgecolor', 'none');
    title(sprintf('t = %g', tv(i)/60))
    set(gca, 'clim', [0 1.2])
    axis square
    drawnow
end
