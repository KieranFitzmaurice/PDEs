% 2-Dimensional Cahn-Hilliard Equation

W = 1;
M = 1;
epsilon = 0.1;

t_0 = 0;
t_f = 5;

h = 0.1;
dt = 0.0003;

% Centered Difference
central_diff = @(x,h) (circshift(x,[0 1]) + circshift(x,[0 -1]) ...
    + circshift(x,[1 0]) + circshift(x,[-1 0]) - 4*x)/h^2;

% Build Cahn-Hilliard Equation
df_dc = @(c) 0.5*W*c.*(1-c).*(1-2*c);
dF_dc = @(c) df_dc(c) - epsilon^2*central_diff(c,h);

Nx = 100;
Ny = 100;

time = t_0:dt:t_f;
concentration = zeros([length(time),Nx,Ny]);

% Initial Conditions

% Specify average concentration and amount of initial noise
c_ave = 0.5;
fluct = 0.01;

c = c_ave*ones([Nx,Ny]) + (2*rand([Nx,Ny]) - 1)*fluct;

for n = 1:length(time)
    c = c + dt*M*central_diff(dF_dc(c),h);
    concentration(n,:,:) = c;
end

for i = 1:200:length(time)
hold on;
pcolor(reshape(concentration(i,:,:),[Nx,Ny]));
colormap jet;
shading flat
axis equal
hold off;
pause(0.001);
end
