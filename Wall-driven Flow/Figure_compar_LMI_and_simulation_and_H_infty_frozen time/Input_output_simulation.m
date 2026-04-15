% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))

% clear; clc;
% FD parameters
N  =  32;       % total FD grid points including boundaries
L  = 2;        % domain length [-1,1]
Re = 500;

% Spectral ranges (as in paper)
kx_list = logspace(-4, 0.48, 3);
kx = 1.2;
kz_list = logspace(-2, 1.2, 3);
kz = 0;
t0=[0, 20, 40, 60, 80, 100];
t0=[0];

P_option = 'full'; % P_option={'full','band','band_12','chord'};
bandwidth = 1; %bandwidth=n-1 will be equivalent to full P matrix. bandwidth=1 will be tridiagonal matrix



% results = cell(length(kx_list), length(kz_list));
% delete(gcp('nocreate'));
% parpool(8);
[operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0, P_option, bandwidth);

%% Main loop
% for i_kx = 1:length(kx_list)
%     [local_results] = Inside_loop(N, L, Re, kz_list, i_kx, kx_list);
%     results(i_kx,:) = local_results;
% end

disp('successful')
save('I_O_sims_and_SVD.mat')

% function [local_results]= Inside_loop(N, L, Re, kz_list, i_kx, kx_list)


%      for j_kz = 1:length(kz_list)

% Build operator for this (kx,kz)
%         kx = kx_list(i_kx);
%         kz = kz_list(j_kz);
%         operator = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0);


%     end
% end


function [operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0, P_option, bandwidth)
local_results = cell(1, length(t0));
syms y t t2 n;
ttt1=0; %TIME START POINT
T=200;
dt=1;
% t0=[0, 20, 40, 60, 80, 100];
% t0=20;
[x, D, D2, D4, w] = finitediff(N, L);
Ny = length(x);   % interior grid points

h=eye(size(D));
h1=eye(size(D));
k2=(kx.^2+kz.^2)*h;
k=(kx.^2+kz.^2).^0.5*h1;
w1=w.^0.5;
W=diag(w1);
d=D+k;
WW=W*d;
WWW=blkdiag(W,W,W);

%v
V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126
c=D2-k2;


%U
Cn=0;
NN=100;
k0 = 0.1;
% Re= 500;
g=1-exp(-k0.*t); %Acc
% g=exp(-k0.*t); %Dec
an = (pi*n).^2./Re;  %lately
m=diff(g,t);
a1=subs(m,t,0);
m2=diff(m,t);
a2=subs(m2,t,t2);
e=exp(an.*t2);
a3=simplify(e*a2);
a4=simplify(int(a3, t2, 0, t),'Steps',200);
a5=a4+a1;
aa=-2.*Re.*(-1).^n./(pi.*n).^3.*(a5)+Cn;
a6=simplify(aa.*exp(-an.*t).*sin(n.*y.*pi),'Steps',100);
b1=Re./6.*m.*(y.^3-y);
b2=g.*y;
bb=b1+b2;
a7 = 0;
for ii = 1:NN
    a7 = a7 + subs(a6, n, ii);
end
U=a7+bb;

mm=diff(U);
mm2=diff(mm);
% U_yi = (subs(U, y, yi(2:N)));
% m_yi = (subs(mm, y, yi(2:N)));
% m2_yi = (subs(mm2, y, yi(2:N)));

U_yi = (subs(U, y, x));
m_yi = (subs(mm, y, x));
m2_yi = (subs(mm2, y, x));
U_yi_function = matlabFunction(U_yi);
m_yi_function = matlabFunction(m_yi);
m2_yi_function = matlabFunction(m2_yi);
o=1/(1i*Re)*(D4-2*D2*k2+k2*k2);
oo=-inv(c);
ooo=(1/(1i*Re))*c;

%% Construct Nonlinear Time-varying System

%     A = zeros(2*(N-2), 2*(N-2), t_steps);%
%U_value
%     G = zeros(1, t_steps);
AA = eye(size(V));
%         U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
%         m_yi_value = diag(m_yi_function(t_value));
%         m2_yi_value = diag(m2_yi_function(t_value));
k2 = kx^2 + kz^2;
% Laplacians
Lap  = D2 - k2*eye(Ny);                     % ∇^2
Lap2 = D4 - 2*k2*D2 + (k2^2)*eye(Ny);       % ∇^4
LHS = blkdiag(Lap, eye(Ny));

t_steps = (T-t0)/dt;
A = zeros(2*(N-2), 2*(N-2), t_steps);%
%U_value
G = zeros(1, t_steps);
AA = eye(size(V));
for i = 1:t_steps  %step
    t_value=t0+i*dt;
    k2 = kx^2 + kz^2;
    % Laplacians
    Lap  = D2 - k2*eye(Ny);                     % ∇^2
    Lap2 = D4 - 2*k2*D2 + (k2^2)*eye(Ny);       % ∇^4
    LHS = blkdiag(Lap, eye(Ny));
    U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
    m_yi_value = diag(m_yi_function(t_value));
    m2_yi_value = diag(m2_yi_function(t_value));

    %L
    Los=oo*(o-kx*U_yi_value*c+kx*m2_yi_value);%!
    Lsq=kx*U_yi_value-ooo;
    Lc=kz*m_yi_value;%diag
    L_value{i}=[Los,zeros(N-2,N-2);Lc,Lsq];
    %         L=V*(-1i*L_value{i})/V;%new added


    %         A(:, :, i) = L;
    %         AA = expm(-1i*L_value{i}*dt)*AA;
    %         gh1=double(V*AA/V);
    %         G(i)=norm(gh1)^2;
end


Bmat = [ -1i*kx*D     -(k2)*eye(Ny)   -1i*kz*D;
    1i*kz*eye(Ny)  zeros(Ny)     -1i*kx*eye(Ny) ];
operator.B = LHS \ Bmat;
operator.B_tilde = Bmat;

% C operator
Cmat = 1/k2 * [ 1i*kx*D   -1i*kz*eye(Ny);
    k2*eye(Ny) zeros(Ny);
    1i*kz*D   1i*kx*eye(Ny) ];
operator.C = Cmat;

%
L_function = @(t)  [ ...
    oo*( o - kx*diag(U_yi_function(t))*c + kx*diag(m2_yi_function(t)) ),  zeros(Ny,Ny); ...
    kz*diag(m_yi_function(t)),                                              kx*diag(U_yi_function(t)) - ooo ...
    ];




nnn = logspace(-4, -0.3, 128);
rand_value = 1;

%% H_infty

% t_s = 0;
% A = -1i * L_value{t_s};
% sys = ss(A, operator.B, operator.C, zeros(size(operator.C,1)));
% Hinf_frozen = norm(sys, inf);        % single scalar
% max_SVD = Hinf_frozen;   % store for this t_s (or just a single value)
% [sv, w] = sigma(sys, nnn);
% gain = max(sv, [], 1);



%% Input-output Simulation

for omega_ind = 1:length(nnn)

    nnn1 = nnn(omega_ind);     % current frequency
    random_lin = rand(3*Ny, rand_value);

    % storage for G over random trials
    G = zeros(rand_value,1);

    for ra = 1:rand_value

%         random_vec = random_lin(:, ra);   % (3*Ny x 1)
% %         random_vec = random_vec / norm(random_vec);
%         % Build symbolic dd if t is symbolic, then convert
%         dd = @(tt) sin(nnn1 * tt) * random_vec;
%         dd = matlabFunction(dd_sym, 'Vars', t);   % dd(t) returns column vector

%         spatial_shape = ones(size(operator.B, 2), 1);         
%         dd = @(tt) sin(nnn1 * tt) * spatial_shape;



%         jwI_minus_A = 1i * nnn1 * eye(size(A)) - A;
%         H_jw = operator.C * (jwI_minus_A \ operator.B);
%         [U_svd, S_svd, V_svd] = svd(H_jw);
%         optimal_vec = V_svd(:, 1);         
%         dd = @(tt) sin(nnn1 * tt) * optimal_vec;
        

%         tspan = linspace(0, 20000, 1000);
        tspan = linspace(0, 200, 1000);
        u0 = zeros(2*Ny, 1);

        % Preallocate for norms
        fai_norm = zeros(length(tspan),1);
        dd_norm  = zeros(length(tspan),1);

        % Define the ODE
%         result = @(tt, uu) -1i * L_function(tt) * uu + operator.B * dd(tt);
        % Use the operator frozen at t_s 
%         result = @(tt, uu) A * uu + operator.B * dd(tt);
        result = @(tt, uu) time_varying_system(tt, uu, nnn1, operator, L_function);
        [t_values3, u_values3] = ode45(result, tspan, u0);


        for j4 = 1:length(t_values3)
            u = u_values3(j4,:)';              % convert row → column
            fai = operator.C * u;

%             fai_norm(j4) = norm(fai .* WWW)^2;
            fai_norm(j4) = norm(fai)^2;

%             dd_val = dd(t_values3(j4));
            [~, dd_val] = get_optimal_input(t_values3(j4), nnn1, operator, L_function);
            dd_norm(j4) = norm(dd_val)^2;
        end

        % L2-like total amplification
        fai_total = sum(fai_norm);
        dd_total  = sum(dd_norm);
        G(ra)     = sqrt(fai_total / dd_total);

    end
    G_norm(omega_ind) = max(G);
end



[G_max_omega,index]=max(G_norm);
figure
plot(nnn,G_norm);
hold on;
plot(nnn, gain);
set(gca, 'XScale', 'log','YScale', 'log');    % log scale
xlabel('ω');
ylabel('gain');
legend('G (simulation)','$\|H\|_\infty$', 'Interpreter', 'latex');
grid on;
saveas(gcf, 'Gnorm_SVD_plot.png');
nnn_max=nnn(index);
disp(['lambda^2=' num2str(G_max_omega, '%e')]);
%% Without SVD
% local_results{1} = struct( ...
%     'kx', kx, ...
%     'kz', kz, ...
%     'G_norm', G, ...
%     'G_max_omega', G_max_omega, ...
%     'nnn_max', nnn_max, ...
%     'nnn', nnn);
%% With SVD
local_results{1} = struct( ...
    'kx', kx, ...
    'kz', kz, ...
    'G_norm', G_norm, ...
    'G_max_omega', G_max_omega, ...
    'nnn_max', nnn_max, ...
    'nnn', nnn, ...
    'max_SVD', max_SVD, ...
    'gain', gain);
end


%     G_storage{j3} = G;
%     operator.max_trans=max(G);


function du = time_varying_system(t, u, omega, operator, L_func)

    A_t = -1i * L_func(t); 
    
    %instantaneous Transfer Function
    jwI_minus_A = 1i * omega * eye(size(A_t)) - A_t;
    H_t = operator.C * (jwI_minus_A \ operator.B);
    
    [~, ~, V_svd] = svd(H_t);
    v_optimal = V_svd(:, 1); 
    d_t = sin(omega * t) * v_optimal;
 
    du = A_t * u + operator.B * d_t;
end

function [v_opt, d_t] = get_optimal_input(t, omega, operator, L_func)
    % Instantaneous Transfer Function
    A_t = -1i * L_func(t);
    jwI_minus_A = 1i * omega * eye(size(A_t)) - A_t;
    H_t = operator.C * (jwI_minus_A \ operator.B);
    
    % SVD for optimal shape
    [~, ~, V_svd] = svd(H_t);
    v_opt = V_svd(:, 1);
    
    % Apply harmonic forcing
    d_t = sin(omega * t) * v_opt;
end




function [x, D1, D2, D4,w] =finitediff(N,L)

%Differentiation matrix using finite difference scheme.
%This is suitable for Dirichlet boundary condition v(x=L/2)=v(x=-L/2)=0 at
%the boundary and Neumann boundary condition v'(x=L/2)=v'(x=-L/2)=0.

dx=L/N;%get grid spacing.

x=linspace(-L/2,L/2,N);%get the grid point location

x=x(2:end-1); %delete the first and last points that are zero.

w=dx*ones(1,N-2);%integration weighting

N_diff=N-2;%The size of differentiation matrices

%First order derivative based on central difference
%f'(x_i)=(f(x_{i+1})-f(x_{i-1})/(2*dx)
%We also use the boundary condition that x_0=x_N=0
D1_stencil=diag(-1*ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1);
D1=D1_stencil/(2*dx);

%Second order derivative based on central difference
%f''(x_i)=(f(x_{i+1})-2f(x_i)+f(x_{i-1})/(dx^2)
%We also use the boundary condition that x_0=x_N=0
D2_stencil=(diag(-2*ones(1,N_diff))+diag(ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1));
D2=D2_stencil/dx^2;

%Forth order derivative based on central difference
%f''''(x_i)=(f(x_{i+2})-4f(x_{i+1})+6f(x_i)-4f(x_{i-1})+f(x_{i-2})/(dx^4)
%This differentiation matrix only go through x_1 up to x_{N-1}
D4_stencil=(diag(6*ones(1,N_diff))+diag(-4*ones(1,N_diff-1),1)+diag(-4*ones(1,N_diff-1),-1)...
    + diag(ones(1,N_diff-2),2)+diag(ones(1,N_diff-2),-2));

%Here, we use the Neumann boundary condition that x_0'=x_N'=0
%such that x_1=x_{-1} and x_{N-1}=x_{N+1} for the ghost points. Then we
%also use the condition x_0 and x_N=0 to express all values based on x_1 up
%to x_{N-1}
D4_stencil(1,1)=D4_stencil(1,1)+1;
D4_stencil(end,end)=D4_stencil(end,end)+1;

D4=D4_stencil/dx^4;
end


function [P,F] = sparse_method(P_option, n, bandwidth)

tic;
if strcmp(P_option,'full')
    %full n*n P matrix
    P = sdpvar(n,n,'hermitian','complex');

    %constraint that P is PSD.
    F = P - 0.01*eye(n) >= 0;

elseif strcmp(P_option,'band')
    %
    P=diag(sdpvar(1,n));
    if bandwidth>=1
        for ind=1:bandwidth
            off_diag=diag(sdpvar(1,n-ind,'full','complex'),ind);
            P=P+off_diag+off_diag';
        end
    end

    %constraint that P is PSD.
    F = P - 0.01*eye(n) >= 0;
elseif strcmp(P_option,'band_12') && bandwidth<n/2
    %
    P=diag(sdpvar(1,n));
    P12=diag(sdpvar(1,n/2,'full','complex'));

    if bandwidth>=1
        for ind=1:bandwidth
            off_diag=diag(sdpvar(1,n-ind,'full','complex'),ind);
            P=P+off_diag+off_diag';

            off_diag_P12_upper=diag(sdpvar(1,n/2-ind,'full','complex'),ind);
            off_diag_P12_lower=diag(sdpvar(1,n/2-ind,'full','complex'),-ind);
            P12=P12+off_diag_P12_upper+off_diag_P12_lower;
        end
    end

    O=zeros(n/2,n/2);
    P=P+[O,P12;
        P12',O];


    %constraint that P is PSD.
    F = P - 0.01*eye(n) >= 0;

elseif strcmp(P_option,'chord')
    %chord decomposition that decompose P into several smaller
    %P PSD in the size of bandwidth.
    F=[];
    P=zeros(n,n);

    for ind=1:n-bandwidth
        %create the E matrix that represent the interaction of
        %sparse matrix
        E_C{ind}=zeros(bandwidth+1,n);
        E_C{ind}(1:bandwidth+1,ind:ind+bandwidth)=eye(bandwidth+1);

        P_C{ind}=sdpvar(bandwidth+1,bandwidth+1,'hermitian','complex');

        F=[F,P_C{ind}-0.01*eye(bandwidth+1)>=0];
        P=P+E_C{ind}'*P_C{ind}*E_C{ind};
    end

else
    error('Wrong P_option');
end



end


function [t, x] = rk5_solver(ode_func, tspan, u0)
% RK5 Solver for ODEs
% Inputs:
%   ode_func - Function handle for the ODE system (f(t, u))
%   tspan - Time vector [t0, t1, ..., tn]
%   u0 - Initial condition (row vector)
% Outputs:
%   t - Time vector (same as tspan)
%   x - Solution matrix (size: length(tspan) x length(u0))
%   forcing - Forcing term evaluated at each step (optional)

% Number of time steps
N = length(tspan);

% Preallocate solution and forcing arrays
x = zeros(N, length(u0));
forcing = zeros(N, length(u0));
x(1, :) = u0; % Set initial condition

% Coefficients for RK5 (Dormand-Prince example)
c = [0, 1/4, 3/8, 12/13, 1, 1/2];
b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
a = [
    0,      0,      0,       0,       0,      0;
    1/4,    0,      0,       0,       0,      0;
    3/32,   9/32,   0,       0,       0,      0;
    1932/2197, -7200/2197, 7296/2197, 0,       0,  0;
    439/216, -8,    3680/513, -845/4104, 0,    0;
    -8/27,   2,    -3544/2565, 1859/4104, -11/40, 0
    ];

% Time-stepping loop
for i = 1:N-1
    h = tspan(i+1) - tspan(i); % Time step size
    u_current = x(i, :)'; % Current state

    % Compute RK5 coefficients (k1 to k6)
    k1 = h * ode_func(tspan(i), u_current);
    k2 = h * ode_func(tspan(i) + c(2)*h, u_current + a(2,1)*k1);
    k3 = h * ode_func(tspan(i) + c(3)*h, u_current + a(3,1)*k1 + a(3,2)*k2);
    k4 = h * ode_func(tspan(i) + c(4)*h, u_current + a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3);
    k5 = h * ode_func(tspan(i) + c(5)*h, u_current + a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4);
    k6 = h * ode_func(tspan(i) + c(6)*h, u_current + a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5);

    % Update state
    x(i+1, :) = (u_current + b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4 + b(5)*k5 + b(6)*k6)';

%     % Optionally, store forcing term for analysis
%     [~, forcing(i+1, :)] = ode_func(tspan(i), u_current);
end

t = tspan; % Output time vector
end
