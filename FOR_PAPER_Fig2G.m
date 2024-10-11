%% Initialize a clean workspace
clear all, close all

%% PART 1: SIMULATE THE PATTERN
% This code simulates migration in response to the circular pattern from
% Figure 2E-G.

% 1.1 Make a uniform initial pattern of cells
ic = ones(101,101); % a 101 x 101 discretized grid of cells

imshow(ic)

% 1.2 Construct an illumination circle centered at the origin

R = 9.35; % 500 um / 5400 um * 101 units = 9.35 unit diameter

x = -50:50;
y = -50:50;
[xx, yy] = meshgrid(x,y);

u = zeros(size(xx));
u((xx.^2+yy.^2) < R^2)=1; 

figure(1)
imagesc(u)
axis equal
pattern = u;

% 1.3 Simulate this model
N = 101;
L = 5400;  % simulation width in um (6 mm = 6,000 um)
tF = 1920; % simulation time in min (32 h = 1,920 min)

f = construct_2D_model(N,L,pattern);

ic = reshape(ic,N^2,1);
tic
[t,y] = ode15s(f, [0 tF], ic);
toc

tv = linspace(0,tF,100);
yv = interp1(t, y, tv);

save Fig2G_simulation_results L N tF tv yv

%% PLOT RESULTS
% This code plots the simulation related to the experimental data of Figure 2G.

load Fig2G_simulation_results

xv = linspace(0,L,N);
[X Y] = meshgrid(xv, xv);

% Get 20 hour data
i20 = find(tv >= 20*60, 1, 'first');
Z20 = reshape(yv(i20,:),N,N);

% Get 32 hour data
i32 = find(tv >= 32*60, 1, 'first');
Z32 = reshape(yv(i32,:),N,N);

% Plot the 20 and 32 hour data as in Fig. 2G
xx = X(51,51:end)-L/2;
plot(xx, Z20(51,51:end), xx, Z32(51,51:end))
xlabel('distance from center')
ylabel('relative cell density')
legend({'20 h' '32 h'}), legend boxoff
set(gca, 'ylim', [0 2])
set(gca, 'xlim', [0 1000])
set(gcf, 'position', [360 290 360 310])
