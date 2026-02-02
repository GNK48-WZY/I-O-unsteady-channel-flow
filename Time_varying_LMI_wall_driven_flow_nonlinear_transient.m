% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))

% clear; clc;
% FD parameters
N  = 4;       % total FD grid points including boundaries
L  = 2;        % domain length [-1,1]
Re = 500;

% Spectral ranges (as in paper)
kx_list = logspace(-4, 0.48, 3);
kx = 1.2;
kz_list = logspace(-2, 1.2, 3);
kz = 0;
t0=[0, 20, 40, 60, 80, 100];
t0=[0];

% delta_list = logspace(-6, 0, 1);
% T = 200;
% t0 = 0;
% dt = 0.1;
% t_steps = (T-t0)/dt;
% results = cell(length(kx_list), length(kz_list));
delete(gcp('nocreate'));
% parpool(8);
[operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0);

%% Main loop
% for i_kx = 1:length(kx_list)
%     [local_results] = Inside_loop(N, L, Re, kz_list, i_kx, kx_list);
%     results(i_kx,:) = local_results;
% end

disp('successful')
save('mu_full_scan_results.mat')

% function [local_results]= Inside_loop(N, L, Re, kz_list, i_kx, kx_list)


%      for j_kz = 1:length(kz_list)

        % Build operator for this (kx,kz)
%         kx = kx_list(i_kx);
%         kz = kz_list(j_kz);
%         operator = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0);

            
%     end
% end


function [operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0)
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

%v
V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126
c=D2-k2;


%U
Cn=0;
NN=100;
k0 = 0.1;
% Re= 500;
g=1-exp(-k0.*t); %Acc
g=exp(-k0.*t); %Dec
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
for j3 = 1:length(t0)
    %A

    t_steps = (T-t0(j3))/dt;
    A = zeros(2*(N-2), 2*(N-2), t_steps);%
    %U_value
    G = zeros(1, t_steps);
    AA = eye(size(V));
    for i = 1:t_steps  %step
        t_value=t0(j3)+i*dt;
%         U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
%         m_yi_value = diag(m_yi_function(t_value));
%         m2_yi_value = diag(m2_yi_function(t_value));
        k2 = kx^2 + kz^2;
        % Laplacians
         Lap  = D2 - k2*eye(Ny);                     % ∇^2
         Lap2 = D4 - 2*k2*D2 + (k2^2)*eye(Ny);       % ∇^4
%         
%         % A operator
%         A11 = -1i.*kx.*diag(U_yi_function(t_value)).*Lap + 1i.*kx.*diag(m2_yi_function(t_value)) + (1/Re).*Lap2;
%         A12 = zeros(Ny);
%         A21 = -1i.*kz.*diag(m_yi_function(t_value));
%         A22 = -1i.*kx.*diag(U_yi_function(t_value)) + (1/Re).*Lap;
%         
         LHS = blkdiag(Lap, eye(Ny));
%         RHS = [A11 A12; A21 A22];
%         operator.A = LHS \ RHS;
%         A(:, :, i) = L; 
%         AA = expm(-1i*operator.A*dt)*AA;
%         gh1=double(V*AA/V);
%         G(i)=norm(gh1)^2;
        U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
        m_yi_value = diag(m_yi_function(t_value));
        m2_yi_value = diag(m2_yi_function(t_value));

        %L
        Los=oo*(o-kx*U_yi_value*c+kx*m2_yi_value);%!
        Lsq=kx*U_yi_value-ooo;
        Lc=kz*m_yi_value;%diag
        L_value=[Los,zeros(N-2,N-2);Lc,Lsq];
        L=V*(-1i*L_value)/V;%new added



        A(:, :, i) = L; 
        AA = expm(-1i*L_value*dt)*AA;
        gh1=double(V*AA/V);
        G(i)=norm(gh1)^2;
    end

%     G_storage{j3} = G;
%     operator.max_trans=max(G);


    Bmat = [ -1i*kx*D     -(k2)*eye(Ny)   -1i*kz*D;
          1i*kz*eye(Ny)  zeros(Ny)     -1i*kx*eye(Ny) ];
    operator.B = LHS \ Bmat;
    
    % C operator
    Cmat = 1/k2 * [ 1i*kx*D   -1i*kz*eye(Ny);
                     k2*eye(Ny) zeros(Ny);
                     1i*kz*D   1i*kx*eye(Ny) ];
    operator.C = Cmat;
    
    % ---- build gradient-output operators ----
    Grad_block = [1i*kx*eye(Ny); D; 1i*kz*eye(Ny)];
    Grad_big   = blkdiag(Grad_block, Grad_block, Grad_block);
    
    C_grad = Grad_big * operator.C;  % (9*Ny) × (2*Ny)
    
    n = size(operator.C, 2); % = 2*Ny
    operator.C_grad_u = C_grad(1:3*Ny, 1:n);
    operator.C_grad_v = C_grad(3*Ny+1:6*Ny, 1:n);
    operator.C_grad_w = C_grad(6*Ny+1:9*Ny, 1:n);

    Ny = size(operator.B,2)/3;
    operator.E = LHS;
%     n  = size(A,1);

%% LMI computation
    yalmip('clear');

   % time-varying
   for j1=1:t_steps
       P{j1}=sdpvar(n,n,'hermitian','complex');
       sx{j1} = sdpvar(1,1); sy{j1} = sdpvar(1,1); sz{j1} = sdpvar(1,1);
   end
%    sx = sdpvar(1,1); sy = sdpvar(1,1); sz = sdpvar(1,1);
%    gamma_H_inf_complex2 = sdpvar(1,1);
%    tt = sdpvar(1);
   G_bar = sdpvar(1,1);
   I = eye(n);
%    scaling = 1;
%    E=I*scaling;
   delta_2 = 1e-16;
%    operator.A_new=operator.A*scaling;
   Gamma = sdpvar(1,1);
   F = [];
%    F = [];
for j1 = 1:t_steps
    % Always include positivity constraints for every time step
    F = [F, ...
        I <= P{j1} <= G_bar * I, ...
        sx{j1} >= 0, ...
        sy{j1} >= 0, ...
        sz{j1} >= 0];

    % Only form and include dV_ineq if j1 < t_steps
    if j1 < t_steps
        M = blkdiag(sx{j1}*eye(Ny), sy{j1}*eye(Ny), sz{j1}*eye(Ny));
        dP_dt = (P{j1 + 1} - P{j1}) / dt;

%         dV_ineq = [ ...
%             dP_dt + operator.A_new' * P{j1} * E + E' * P{j1} * operator.A_new ...
%             + sx{j1} * operator.C_grad_u' * operator.C_grad_u ...
%             + sy{j1} * operator.C_grad_v' * operator.C_grad_v ...
%             + sz{j1} * operator.C_grad_w' * operator.C_grad_w, ...
%             P{j1} * operator.B;
%             operator.B' * P{j1}, -Gamma * M ...
%         ];
%         dV_ineq = [ ...
%             dP_dt + A(:,:,j1)' * P{j1} * E + E' * P{j1} * A(:,:,j1) ...
%             + sx{j1} * operator.C_grad_u' * operator.C_grad_u ...
%             + sy{j1} * operator.C_grad_v' * operator.C_grad_v ...
%             + sz{j1} * operator.C_grad_w' * operator.C_grad_w, ...
%             P{j1} * operator.B;
%             operator.B' * P{j1}, -Gamma * M ...
%         ];
        dV_ineq = [ ...
            dP_dt + A(:,:,j1)' * P{j1} * operator.E + operator.E' * P{j1} * A(:,:,j1) ...
            + sx{j1} * delta_2 *operator.C_grad_u' * operator.C_grad_u ...
            + sy{j1}* delta_2*operator.C_grad_v' * operator.C_grad_v ...
            + sz{j1} * delta_2*operator.C_grad_w' * operator.C_grad_w, ...
            P{j1} * operator.B;
            operator.B' * P{j1}, -M ...
        ];

        % Add inequality constraint
        F = [F, dV_ineq <= 0];
    end
end

%    for j1=1:t_steps
%        M = blkdiag(sx{j1}*eye(Ny), sy{j1}*eye(Ny), sz{j1}*eye(Ny));
%        dP_dt = (P{j1 + 1}-P{j1})/dt;
%        dV_ineq{j1} = [dP_dt + operator.A_new'*P{j1}*E + E'*P{j1}*operator.A_new ...
%                         + sx{j1}*operator.C_grad_u'*operator.C_grad_u ...
%                         + sy{j1}*operator.C_grad_v'*operator.C_grad_v ...
%                         + sz{j1}*operator.C_grad_w'*operator.C_grad_w, ...
%                         P{j1}*operator.B;
%                         operator.B'*P{j1}, -Gamma*M];     
%        F = [F, P{j1} - tt*eye(n) >= 0, dV_ineq{j1} <= 0, sx{j1} >= 0, sy{j1} >= 0, sz{j1} >= 0];
% 
% %         F = [I <= P{j1} <= tt*I,  dV_ineq <= 0, sx >= 0, sy >= 0, sz >= 0];
%    end
% %    for j1 = 1:t_steps
% %         F = [F, P{j1} - tt*eye(n) >= 0, dV_ineq{j1} <= 0, sx{j1} >= 0, sy{j1} >= 0, sz{j1} >= 0];
% %    end
%     F=[F];

    sdp_options = sdpsettings('solver','mosek','verbose',1);
%     sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
     diagnostics = bisection(F,G_bar,sdp_options);

    diagnostics.info

     if diagnostics.problem == 0
        operator.G_bar = value(G_bar);
        operator.P_optimal = value(P);
        operator.P = P;
%         operator.sx =sx;
%         operator.sy = sy;
%         operator.sz = sz;
        sx_numeric = value(sx);
        sy_numeric = value(sy);
        sz_numeric = value(sz);
        eigvals = cell(t_steps, 1);
        eigvecs = cell(t_steps, 1);
        for j1 = 1:t_steps
            P_numeric = double(operator.P_optimal{j1});
            P_cell{j1} = P_numeric;
            [V_eig, D_eig] = eig(P_numeric);
            eigvals{j1} = diag(D_eig);
            eigvecs{j1} = V_eig;

        end
        
        all_eigenvalues{j3} = eigvals;
        all_eigenvectors{j3} = eigvecs;
       

     end
     local_results{1, j3} = struct( ...
                'kx', kx, ...
                'kz', kz, ...
                'G_bar', operator.G_bar, ...
                'P', P_cell, ...
                'sx', sx_numeric, ...
                'sy', sy_numeric, ...
                'sz', sz_numeric);
  
end

operator.local_results = local_results;  % embed optional

% figure;
% hold on;
% cmap = parula(length(t0));  % Colormap with 9 colors (one per t0 value)
% 
% % Plot each curve with color based on t0 value
% for j3 = 1:length(t0)
%     if ~isempty(all_eigenvalues{j3})
%         time_points = t0(j3) + (1:length(all_eigenvalues{j3}))*dt;
%         min_eig = cellfun(@min, all_eigenvalues{j3});
%         semilogy(time_points, min_eig, 'Color', cmap(j3,:), 'LineWidth', 1.5);
%     end
% end
% 
% % Configure colorbar with exact t0 values
% c = colorbar;
% colormap(cmap);
% 
% % Set discrete ticks at centers of each color segment
% c.Ticks = linspace(0, 1, length(t0));  % 9 ticks from 0 to 1
% c.TickLabels = arrayfun(@(x) num2str(x), t0, 'UniformOutput', false);  % Exact t0 values
% 
% legend({'min \lambda(P)'}, 'Location', 'best');
% 
% % Axis labels and formatting
% xlabel('t');
% ylabel('Minimum Eigenvalue of \bf{P}');
% grid on;
% box on;
% hold off;
% 
% 
% figure;
% hold on;
% cmap = parula(length(t0));  % Colormap with 9 colors (one per t0 value)
% 
% % Plot each curve with color based on t0 value
% for j3 = 1:length(t0)
%     if ~isempty(all_eigenvalues{j3})
%         time_points = t0(j3) + (1:length(all_eigenvalues{j3}))*dt;
%         max_eig = cellfun(@max, all_eigenvalues{j3});
%         semilogy(time_points, max_eig, '--', 'Color', cmap(j3,:), 'LineWidth', 1.5);
%     end
% end
% 
% % Configure colorbar with exact t0 values
% c = colorbar;
% colormap(cmap);
% 
% % Set discrete ticks at centers of each color segment
% c.Ticks = linspace(0, 1, length(t0));  % 9 ticks from 0 to 1
% c.TickLabels = arrayfun(@(x) num2str(x), t0, 'UniformOutput', false);  % Exact t0 values
% 
% legend({'max \lambda(P)'}, 'Location', 'best');
% 
% % Axis labels and formatting
% xlabel('t');
% ylabel('Maximum Eigenvalue of \bf{P}');
% grid on;
% box on;
% hold off;


% 
% k2 = kx^2 + kz^2;
% % Laplacians
% Lap  = Dyy - k2*eye(Ny);                     % ∇^2
% Lap2 = D4 - 2*k2*Dyy + (k2^2)*eye(Ny);       % ∇^4
% 
% % A operator
% A11 = -1i*kx*diag(U_yi)*Lap + 1i*kx*diag(m2_yi) + (1/Re)*Lap2;
% A12 = zeros(Ny);
% A21 = -1i*kz*diag(m_yi);
% A22 = -1i*kx*diag(U_yi) + (1/Re)*Lap;
% 
% LHS = blkdiag(Lap, eye(Ny));
% RHS = [A11 A12; A21 A22];
% operator.A = LHS \ RHS;

% B operator


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
