% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))

clear;
% FD parameters
N  = 4;       % total FD grid points including boundaries
L  = 2;        % domain length [-1,1]
Re = 500;
flowType = 'couette';

% Spectral ranges (as in paper)
kx_list = logspace(-4, 0.48, 1);
kx_list = 1.2;
kz_list = logspace(-2, 1.2, 1);
kz_list = 0;

P_option = 'full'; % P_option={'full','band','band_12','chord'};
bandwidth = 10; %bandwidth=n-1 will be equivalent to full P matrix. bandwidth=1 will be tridiagonal matrix


% Delta list (can be more values if needed)
% delta_list = logspace(-6, 0, 1);
% T = 200;
% t0 = 0;
% dt = 0.1;
% t_steps = (T-t0)/dt;
results = cell(length(kx_list), length(kz_list));
delete(gcp('nocreate'));
% parpool(8);

%% Main loop
% time=zeros(i_kx,1);
for i_kx = 1:length(kx_list)
    [local_results] = Inside_loop(N, L, Re, kz_list, flowType, i_kx, kx_list,P_option,bandwidth);
    results(i_kx,:) = local_results;
end

%% Save everything
disp('successful')
% save('mu_full_scan_results.mat','results','kx_list','kz_list')

function [local_results]= Inside_loop(N, L, Re, kz_list, flowType,i_kx, kx_list,P_option,bandwidth)
     local_results = cell(1, length(kz_list));
     for j_kz = 1:length(kz_list)

        % Build operator for this (kx,kz)
        kx = kx_list(i_kx);
        kz = kz_list(j_kz);
        operator = build_operators_fd(N, L, Re, kx, kz, flowType);

        Ny = size(operator.B,2)/3;
        n  = size(operator.A,1);

%         for ind_delta = 1:length(delta_list)
%             delta = delta_list(ind_delta);

            %% LMI computation
            yalmip('clear');

            type = "time_varying_Hinf_sparse";  
          

  switch type
      case "time_varying"
            t_steps = 200;
            dt = 1;
            for j1=1:t_steps
                   P{j1}=sdpvar(n,n,'hermitian','complex');
                   sx{j1} = sdpvar(1,1); sy{j1} = sdpvar(1,1); sz{j1} = sdpvar(1,1);
            end
%             sx = sdpvar(1,1); sy = sdpvar(1,1); sz = sdpvar(1,1);
            gamma_H_inf_complex2 = sdpvar(1,1);
            F = [];
            for j1=1:t_steps
                F = [F, P{j1} - 0.01*eye(n) >= 0, sx{j1} >= 0, sy{j1} >= 0, sz{j1} >= 0];   
                if j1 < t_steps
                    M = blkdiag(sx{j1}*eye(Ny), sy{j1}*eye(Ny), sz{j1}*eye(Ny));
                    dP_dt = (P{j1 + 1}-P{j1})/dt;
    %                 dV_ineq = [dP_dt + operator.A'*P{j1} + P{j1}*operator.A ...
    %                         + sx{j1}*operator.C_grad_u'*operator.C_grad_u ...
    %                         + sy{j1}*operator.C_grad_v'*operator.C_grad_v ...
    %                         + sz{j1}*operator.C_grad_w'*operator.C_grad_w, ...
    %                         P{j1}*operator.B;
    %                         operator.B'*P{j1}, -gamma_H_inf_complex2*M];   
                    dV_ineq = [dP_dt + operator.A_tilde'*P{j1}*operator.E + operator.E'*P{j1}*operator.A_tilde ...
                            + sx{j1}*operator.C_grad_u'*operator.C_grad_u ...
                            + sy{j1}*operator.C_grad_v'*operator.C_grad_v ...
                            + sz{j1}*operator.C_grad_w'*operator.C_grad_w, ...
                            operator.E'*P{j1}*operator.B_tilde;
                            operator.B_tilde'*P{j1}*operator.E, -gamma_H_inf_complex2*M];   
                      
                    F = [F, dV_ineq <= 0];
                end
            end
            F=[F];
            sdp_options = sdpsettings('solver','mosek','verbose',1);
            yalmip_output = bisection(F,gamma_H_inf_complex2,sdp_options); 
            mu_LMI = sqrt(value(gamma_H_inf_complex2));

      case "time_varying_sparse"
            t_steps = 200;
            dt = 1;
            F = [];
            for j1=1:t_steps
                   [P{j1},F_local] = sparse_method(P_option, n, bandwidth);
                   F = [F, F_local];
                   sx{j1} = sdpvar(1,1); sy{j1} = sdpvar(1,1); sz{j1} = sdpvar(1,1);
            end
            gamma_H_inf_complex2 = sdpvar(1,1);
            for j1=1:t_steps-1
                M = blkdiag(sx{j1}*eye(Ny), sy{j1}*eye(Ny), sz{j1}*eye(Ny));
                dP_dt = (P{j1 + 1}-P{j1})/dt;
                dV_ineq = [dP_dt + operator.A_tilde'*P{j1}*operator.E + operator.E*P{j1}*operator.A_tilde ...
                        + sx{j1}*operator.C_grad_u'*operator.C_grad_u ...
                        + sy{j1}*operator.C_grad_v'*operator.C_grad_v ...
                        + sz{j1}*operator.C_grad_w'*operator.C_grad_w, ...
                        operator.E'*P{j1}*operator.B_tilde;
                        operator.B_tilde'*P{j1}*operator.E, -gamma_H_inf_complex2*M];     
                F = [F, dV_ineq <= 0, sx{j1} >= 0, sy{j1} >= 0, sz{j1} >= 0];    
            end
                       
           F=[F];
           sdp_options = sdpsettings('solver','mosek','verbose',1);
           yalmipt_output = bisection(F,gamma_H_inf_complex2,sdp_options); 
           mu_LMI = sqrt(value(gamma_H_inf_complex2));

            %% time-varying for H_inf norm
      case "time_varying_Hinf_sparse"

            t_steps = 200;
            dt = 1;
            F = [];
            for j1=1:t_steps
                   [P{j1},F_local] = sparse_method(P_option, n, bandwidth);
                   F = [F, F_local];
            end
            sx = 1; sy = 1; sz = 1;
            gamma_H_inf_complex2 = sdpvar(1,1);
            for j1 = 1:t_steps-1
                M = blkdiag(sx*eye(Ny), sy*eye(Ny), sz*eye(Ny));
                dP_dt = (P{j1 + 1}-P{j1})/dt;
                dV_ineq = [dP_dt + operator.A_tilde'*P{j1}*operator.E + operator.E'*P{j1}*operator.A_tilde ...
                        + operator.C'*operator.C, ...
                        operator.E'*P{j1}*operator.B_tilde;
                        operator.B_tilde'*P{j1}*operator.E, -gamma_H_inf_complex2*M];
                
            
                F = [F, dV_ineq <= 0];
            end
            F = [F];
            sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            yalmip_output = bisection(F,gamma_H_inf_complex2,sdp_options);
            
            mu_LMI = sqrt(value(gamma_H_inf_complex2));
            
            
          %% time-independent(Figure 4a)  
      case "time-independent_sparse"

            [P,F] = sparse_method(P_option, n, bandwidth);
            sx = sdpvar(1,1); sy = sdpvar(1,1); sz = sdpvar(1,1);
            gamma_H_inf_complex2 = sdpvar(1,1);
            
            M = blkdiag(sx*eye(Ny), sy*eye(Ny), sz*eye(Ny));
            
            dV_ineq = [operator.A_tilde'*P*operator.E + operator.E*P*operator.A_tilde ...
                        + sx*operator.C_grad_u'*operator.C_grad_u ...
                        + sy*operator.C_grad_v'*operator.C_grad_v ...
                        + sz*operator.C_grad_w'*operator.C_grad_w, ...
                        operator.E'*P*operator.B_tilde;
                        operator.B_tilde'*P*operator.E, -gamma_H_inf_complex2*M];
            
            F = [F, dV_ineq <= 0, sx >= 0, sy >= 0, sz >= 0];
            
            sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            yalmip_output = bisection(F,gamma_H_inf_complex2,sdp_options);
            
            mu_LMI = sqrt(value(gamma_H_inf_complex2));

            
             %% (Figure 4b) with Gradient C and sx,sy,sz = 1
     case "Figure_4b"
            P = sdpvar(n,n,'hermitian','complex');
            sx = 1; sy = 1; sz = 1;
            gamma_H_inf_complex2 = sdpvar(1,1);
            
            M = blkdiag(sx*eye(Ny), sy*eye(Ny), sz*eye(Ny));
            
            dV_ineq = [operator.A'*P + P*operator.A ...
                        + sx*operator.C_grad_u'*operator.C_grad_u ...
                        + sy*operator.C_grad_v'*operator.C_grad_v ...
                        + sz*operator.C_grad_w'*operator.C_grad_w, ...
                        P*operator.B;
                        operator.B'*P, -gamma_H_inf_complex2*M];
            
            F = [P - 1e-3*eye(n) >= 0, dV_ineq <= 0, sx >= 0, sy >= 0, sz >= 0];
            
            sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            yalmip_output = bisection(F,gamma_H_inf_complex2,sdp_options);
            
            mu_LMI = sqrt(value(gamma_H_inf_complex2));
            
            %% (Figure 4c) without Gradient C and sx,sy,sz = 1 
      case "Figure_4c_sparse"
             [P,F] = sparse_method(P_option, n, bandwidth);
%             P = sdpvar(n,n,'hermitian','complex');
            sx = 1; sy = 1; sz = 1;
            gamma_H_inf_complex2 = sdpvar(1,1);
%             F = [];
            M = blkdiag(sx*eye(Ny), sy*eye(Ny), sz*eye(Ny));
            % %Old version, we move Laplacian inverse to the RHS
%             dV_ineq = [operator.A'*P + P*operator.A ...
%                         + operator.C'*operator.C, ... 
%                         P*operator.B;
%                         operator.B'*P, -gamma_H_inf_complex2*M];
            %new version:
            dV_ineq = [operator.A_tilde'*P*operator.E + operator.E*P*operator.A_tilde ...
                        + operator.C'*operator.C, ...
                        operator.E'*P*operator.B_tilde;
                        operator.B_tilde'*P*operator.E, -gamma_H_inf_complex2*M];

            F = [F, dV_ineq <= 0];

%             sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            sdp_options = sdpsettings('solver','mosek','verbose',1);
        
            % sdp_options = sdpsettings('bisection.solver','sparsecolo','verbose',2);
            % sdp_options.sparsecolo.SDPsolver='mosek';
            % sdp_options.sparsecolo.domain=0;
            
            yalmip_output=bisection(F,gamma_H_inf_complex2,sdp_options);

            mu_LMI_2 = value(gamma_H_inf_complex2);
            mu_LMI = sqrt(mu_LMI_2);

            
  end
            time_LMI=toc;

            %% H_inf norm; old one, with inverse(Laplacian) on the RHS

            % sys = ss(operator.A, operator.B, operator.C, zeros(size(operator.C,1), size(operator.B,2)));
            % hinf_norm = norm(sys, inf);

            %% H_inf norm, new one, with Laplacian in the LHS as E operator. This will need to use dss (descriptor state-space model)
            tic;
            sys = dss(operator.A_tilde, operator.B_tilde, operator.C, zeros(size(operator.C,1), size(operator.B,2)),operator.E);
            hinf_norm = hinfnorm(sys);
            time_hinf = toc;
%             G = C * inv(1i * omega * I - A) * B;

            %% === Mussv computation (pos/neg freq) ===
            tic;
            C_all = [operator.C_grad_u;
             operator.C_grad_v;
             operator.C_grad_w];

            BlockStructure=[Ny,3*Ny;
            Ny,3*Ny;
            Ny,3*Ny];

            sys_grad = dss(operator.A_tilde, operator.B_tilde, C_all, zeros(size(C_all,1), size(operator.B,2)),operator.E);
            [mu_frd_H_inf_grad, mu_info] = mussv(sys_grad, BlockStructure, 'Ufs');
            [mu_upper_H_inf_grad, max_ind] = max(squeeze(mu_frd_H_inf_grad(1).ResponseData));
            mu_upper_freq = mu_frd_H_inf_grad(1).Frequency(max_ind);
            % Negative frequency μ-analysis
            K_operator_discretized = eye(size(C_all,1));   % placeholder
            sys_grad_neg_freq = dss(-operator.A_tilde, -operator.B_tilde, K_operator_discretized*C_all, ...
                              zeros(size(C_all,1), size(operator.B,2)),operator.E);
            [mu_frd_H_inf_grad_neg_freq, mu_info_neg] = mussv(sys_grad_neg_freq, BlockStructure, 'Ufs');
            [mu_upper_H_inf_grad_neg_freq, max_ind]   = max(squeeze(mu_frd_H_inf_grad_neg_freq(1).ResponseData));
            mu_upper_freq_neg_freq = -mu_frd_H_inf_grad_neg_freq(1).Frequency(max_ind);

            time_mu=toc;

            local_results{1, j_kz} = struct( ...
                'kx', kx, ...
                'kz', kz, ...
                'mu_LMI', mu_LMI, ...
                'P', value(P), ...
                'sx', value(sx), ...
                'sy', value(sy), ...
                'sz', value(sz), ...
                'mu_upper_H_inf_grad', mu_upper_H_inf_grad, ...
                'mu_upper_freq', mu_upper_freq, ...
                'mu_upper_H_inf_grad_neg_freq', mu_upper_H_inf_grad_neg_freq, ...
                'mu_upper_freq_neg_freq', mu_upper_freq_neg_freq, ...
                'yalmip_output', yalmip_output, ...
                'hinf_norm', hinf_norm, ...
                'time_LMI',time_LMI,...
                'time_hinf',time_hinf,...
                'time_mu',time_mu...
            );
%             local_results{i_kx,j_kz,ind_delta}.kx = kx;
%             local_results{i_kx,j_kz,ind_delta}.kz = kz;
%             local_results{i_kx,j_kz,ind_delta}.delta = delta;
% 
%             local_results{i_kx,j_kz,ind_delta}.mu_LMI   = mu_LMI;
%             local_results{i_kx,j_kz,ind_delta}.P        = value(P);
%             local_results{i_kx,j_kz,ind_delta}.sx       = value(sx);
%             local_results{i_kx,j_kz,ind_delta}.sy       = value(sy);
%             local_results{i_kx,j_kz,ind_delta}.sz       = value(sz);
% 
%             local_results{i_kx,j_kz,ind_delta}.mu_upper_H_inf_grad = mu_upper_H_inf_grad;
%             local_results{i_kx,j_kz,ind_delta}.mu_upper_freq    = mu_upper_freq;
%             local_results{i_kx,j_kz,ind_delta}.mu_upper_H_inf_grad_neg_freq = mu_upper_H_inf_grad_neg_freq;
%             local_results{i_kx,j_kz,ind_delta}.mu_upper_freq_neg_freq = mu_upper_freq_neg_freq;
        
    end
end
function operator = build_operators_fd(N, L, Re, kx, kz, flowType)
% Build state-space operators (A, B, C, C_grad_u/v/w)
% using finite-difference differentiation matrices.
%
% Inputs:
%   N        : total grid points (including boundaries)
%   L        : domain length in wall-normal direction
%   Re       : Reynolds number
%   kx, kz   : streamwise and spanwise wavenumbers
%   flowType : 'Couette' or 'Poiseuille'
%
% Output:
%   operator struct with fields:
%       A, B, C, C_grad_u, C_grad_v, C_grad_w

% ---- get finite-difference matrices ----
[x, Dy, Dyy, D4, ~] = finitediff(N, L);
Ny = length(x);   % interior grid points

% ---- base flow and derivatives ----
switch lower(flowType)
    case 'couette'
        U   = x;                % U(y) = y
        Uy  = ones(size(x));    % U'(y) = 1
        Uyy = zeros(size(x));   % U''(y) = 0
    case 'poiseuille'
        U   = 1 - x.^2;         % U(y) = 1 - y^2
        Uy  = -2*x;             % U'(y) = -2y
        Uyy = -2*ones(size(x)); % U''(y) = -2
    otherwise
        error('flowType must be ''Couette'' or ''Poiseuille''');
end
syms y t
omeg = 10; % Stokes Number
I_2 = eye(2);
Z_2 = zeros(2);
H = [I_2 Z_2];
G = [0 0 1 0;
    0 0 0 1;
    0 -omeg 0 0;
    omeg 0 0 0];
b = [1;0;0;0];
N_1 = [I_2 Z_2;
    Z_2 Z_2];
N_2 = [Z_2 Z_2;
    I_2 Z_2];

p = H * exp(G*(y+1)) * (N_1 + N_2*exp(2*G))^(-1) * b;
W_s = p(1);
W_c = p(2);
R_w_divide_R = 1e-4;
omega_zero = 2*pi/L;

W_yi = 2*R_w_divide_R*(W_c*cos(omeg*t/Re) + W_s*sin(omeg*t/Re));
W_m_yi = diff(W_yi, y);  
W_m2_yi = diff(W_m_yi, y);
W_yi_sub  = subs(W_yi,  y, x);
W_m_yi_sub  = subs(W_m_yi,  y, x);
W_m2_yi_sub = subs(W_m2_yi, y, x);

W_yi_function  = matlabFunction(W_yi_sub);
W_m_yi_function  = matlabFunction(W_m_yi_sub);
W_m2_yi_function = matlabFunction(W_m2_yi_sub);

k2 = kx^2 + kz^2;
% Laplacians
Lap  = Dyy - k2*eye(Ny);                     % ∇^2
Lap2 = D4 - 2*k2*Dyy + (k2^2)*eye(Ny);       % ∇^4

% A operator
A11 = -1i*kx*diag(U)*Lap + 1i*kx*diag(Uyy) + (1/Re)*Lap2;
A12 = zeros(Ny);
A21 = -1i*kz*diag(Uy);
A22 = -1i*kx*diag(U) + (1/Re)*Lap;

LHS = blkdiag(Lap, eye(Ny));
RHS = [A11 A12; A21 A22];
operator.A = LHS \ RHS;

operator.A_tilde=RHS;
operator.E=LHS;

% B operator
Bmat = [ -1i*kx*Dy     -(k2)*eye(Ny)   -1i*kz*Dy;
          1i*kz*eye(Ny)  zeros(Ny)     -1i*kx*eye(Ny) ];
operator.B = LHS \ Bmat;

operator.B_tilde=Bmat;


% C operator
Cmat = 1/k2 * [ 1i*kx*Dy   -1i*kz*eye(Ny);
                 k2*eye(Ny) zeros(Ny);
                 1i*kz*Dy   1i*kx*eye(Ny) ];
operator.C = Cmat;

% ---- build gradient-output operators ----
Grad_block = [1i*kx*eye(Ny); Dy; 1i*kz*eye(Ny)];
Grad_big   = blkdiag(Grad_block, Grad_block, Grad_block);

C_grad = Grad_big * operator.C;  % (9*Ny) × (2*Ny)

n = size(operator.C, 2); % = 2*Ny
operator.C_grad_u = C_grad(1:3*Ny, 1:n);
operator.C_grad_v = C_grad(3*Ny+1:6*Ny, 1:n);
operator.C_grad_w = C_grad(6*Ny+1:9*Ny, 1:n);

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