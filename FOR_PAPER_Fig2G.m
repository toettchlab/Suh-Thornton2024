%% Make a uniform initial pattern of cells
ic = ones(101,101); % a 101 x 101 discretized grid of cells

imshow(ic)

%%
% Construct illumination circle centered at the origin

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

%% Simulate with "outgrowth" initial conditions
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

return

%%
load Fig2G_simulation_results

xv = linspace(0,L,N);
[X Y] = meshgrid(xv, xv);
i20 = find(tv >= 20*60, 1, 'first');
Z20 = reshape(yv(i20,:),N,N);
i32 = find(tv >= 32*60, 1, 'first');
Z32 = reshape(yv(i32,:),N,N);
xx = X(51,51:end)-L/2;
plot(xx, Z20(51,51:end), xx, Z32(51,51:end))
xlabel('distance from center')
ylabel('relative cell density')
set(gca, 'ylim', [0 2])
set(gca, 'xlim', [0 1000])