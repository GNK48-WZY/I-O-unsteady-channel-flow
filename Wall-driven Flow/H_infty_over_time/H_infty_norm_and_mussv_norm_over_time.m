N  = 16;       % total FD grid points including boundaries
L  = 2;        % domain length [-1,1]
Re = 500;

% Spectral ranges (as in paper)
kx_list = logspace(-4, 0.48, 3);
kx = 0;
kz_list = logspace(-2, 1.2, 3);
kz = 1.6;
t0=[0, 20, 40, 60, 80, 100];
t0=[0];

P_option = 'full'; % P_option={'full','band','band_12','chord'};
bandwidth = 4; %bandwidth=n-1 will be equivalent to full P matrix. bandwidth=1 will be tridiagonal matrix
% delete(gcp('nocreate'));
% parpool(8);
[operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0, P_option, bandwidth);

%% Main loop
% for i_kx = 1:length(kx_list)
%     [local_results] = Inside_loop(N, L, Re, kz_list, i_kx, kx_list);
%     results(i_kx,:) = local_results;
% end

disp('successful')
save('mu_full_scan_results.mat', 'operator', 'local_results', '-v7.3');

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

%v
V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126
c=D2-k2; % Laplace


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
oo=-inv(c); % inverse laplace
ooo=(1/(1i*Re))*c;

%% Construct linear Time-varying System
for j3 = 1:length(t0)
    %A

    t_steps = (T-t0(j3))/dt;
    A = zeros(2*(N-2), 2*(N-2), t_steps);%
    %U_value
    G = zeros(1, t_steps);
    AA = eye(size(V));
    for i = 1:t_steps  %step
        t_value=t0(j3)+i*dt;
        k2 = kx^2 + kz^2;
        % Laplacians
        Lap  = D2 - k2*eye(Ny);                     % ∇^2
        Lap2 = D4 - 2*k2*D2 + (k2^2)*eye(Ny);       % ∇^4
        LHS = blkdiag(Lap, eye(Ny));
        U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
        m_yi_value = diag(m_yi_function(t_value));
        m2_yi_value = diag(m2_yi_function(t_value));

        %L
%         Los=oo*(o-kx*U_yi_value*c+kx*m2_yi_value);%!
        Los = -(c \ (o - kx*U_yi_value*c + kx*m2_yi_value));
        Lsq=kx*U_yi_value-ooo;
        Lc=kz*m_yi_value;%diag
        L_value{i}=[Los,zeros(N-2,N-2);Lc,Lsq];
        L=V*(-1i*L_value{i})/V;%new added


        A(:, :, i) = L;
        AA = expm(-1i*L_value{i}*dt)*AA;
        gh1=double(V*AA/V);
        G(i)=norm(gh1)^2;
    end

    %     G_storage{j3} = G;
    %     operator.max_trans=max(G);

    
    Bmat = [ -1i*kx*D     -(k2)*eye(Ny)   -1i*kz*D;
             1i*kz*eye(Ny)  zeros(Ny)     -1i*kx*eye(Ny) ];
    
    B0 = LHS \ Bmat;
    
    C0 = 1/k2 * [ 1i*kx*D   -1i*kz*eye(Ny);
                  k2*eye(Ny) zeros(Ny);
                  1i*kz*D   1i*kx*eye(Ny) ];
    
    operator.B = V * B0;
%     operator.B = B0;
    operator.C = C0 / V;
%     operator.C = C0;
    
    Grad_block = [1i*kx*eye(Ny); D; 1i*kz*eye(Ny)];
    Grad_big   = blkdiag(Grad_block, Grad_block, Grad_block);
    
    C_grad = Grad_big * operator.C;
    
    n = size(A,1);
    operator.C_grad_u = C_grad(1:3*Ny, 1:n);
    operator.C_grad_v = C_grad(3*Ny+1:6*Ny, 1:n);
    operator.C_grad_w = C_grad(6*Ny+1:9*Ny, 1:n);
%     Ny = size(operator.B,2)/3;
%     operator.E=  LHS;

%% (H_infty norm and mussv norm)

            % Use nan so failed steps are obvious
            hinf_norm               = nan(1, t_steps);
            hinf_norm_grad_t        = nan(1, t_steps);

            mu_upper_H_inf_grad_t   = nan(1, t_steps);
            mu_upper_freq_t         = nan(1, t_steps);

            mu_upper_H_inf_grad_neg_t = nan(1, t_steps);
            mu_upper_freq_neg_t       = nan(1, t_steps);

            % Time vector consistent with your current loop:
            % since t_value = t0(j3) + i*dt, i=1,...,t_steps
            ts = t0(j3) + dt*(1:t_steps);

            % Optional: save full sigma/mu curves too
            sigma_gain_cell = cell(1, t_steps);
            sigma_freq_cell = cell(1, t_steps);
            mu_curve_cell   = cell(1, t_steps);
            mu_freq_cell    = cell(1, t_steps);

            w_range = logspace(-4, -0.3, 128);

            C_all = [operator.C_grad_u;
                operator.C_grad_v;
                operator.C_grad_w];

            BlockStructure = [Ny, 3*Ny;
                Ny, 3*Ny;
                Ny, 3*Ny];

            for j1 = 1:t_steps
                disp('good')
                A_1(:,:,j1) = V * (-1i * L_value{j1}) / V;

                % ----- Standard H-infinity norm -----
                sys = ss(A_1(:,:,j1), operator.B, operator.C, ...
                    zeros(size(operator.C,1), size(operator.B,2)));
                hinf_norm(j1) = norm(sys, inf);

                % Optional: save sigma curve
                [sv, w] = sigma(sys, w_range);
                sigma_gain_cell{j1} = sv;
                sigma_freq_cell{j1} = w;

                % ----- Gradient-output system -----
                sys_grad = ss(A_1(:,:,j1), operator.B, C_all, ...
                    zeros(size(C_all,1), size(operator.B,2)));

                % H-infinity norm for gradient output
                hinf_norm_grad_t(j1) = hinfnorm(sys_grad);

                % ----- mu upper bound over chosen frequency grid -----
                % Better to force the same frequency grid at every time step
                sys_grad_frd = frd(sys_grad, w_range);
                [mu_frd_H_inf_grad, mu_info] = mussv(sys_grad_frd, BlockStructure, 'Ufs');

                mu_curve = squeeze(mu_frd_H_inf_grad(1).ResponseData);
                mu_freq  = mu_frd_H_inf_grad(1).Frequency;

                mu_curve_cell{j1} = mu_curve;
                mu_freq_cell{j1}  = mu_freq;

                [mu_upper_H_inf_grad_t(j1), max_ind] = max(mu_curve);
                mu_upper_freq_t(j1) = mu_freq(max_ind);

%                 % ----- Negative-frequency branch, if you still want it -----
%                 sys_grad_neg = ss(-A_1(:,:,j1), -operator.B, C_all, ...
%                     zeros(size(C_all,1), size(operator.B,2)));
%                 sys_grad_neg_frd = frd(sys_grad_neg, w_range);
% 
%                 [mu_frd_H_inf_grad_neg, mu_info_neg] = mussv(sys_grad_neg_frd, BlockStructure, 'Ufs');
% 
%                 mu_curve_neg = squeeze(mu_frd_H_inf_grad_neg(1).ResponseData);
%                 mu_freq_neg  = mu_frd_H_inf_grad_neg(1).Frequency;
% 
%                 [mu_upper_H_inf_grad_neg_t(j1), max_ind_neg] = max(mu_curve_neg);
%                 mu_upper_freq_neg_t(j1) = -mu_freq_neg(max_ind_neg);
            end

            % Peaks over time
            [Hinf_max, peak_idx] = max(hinf_norm);
            t_peak = ts(peak_idx);

            [mu_max, mu_peak_idx] = max(mu_upper_H_inf_grad_t);
            t_mu_peak = ts(mu_peak_idx);
%% Plot

%             figure('Color','w');
% 
%             h1 = plot(ts, hinf_norm, '-', 'LineWidth', 1.2);
%             hold on;
%             h2 = plot(ts, mu_upper_H_inf_grad_t, '-', 'LineWidth', 1.2);
%             
%             % Peak points
%             plot(t_peak, Hinf_max, 'o', ...
%                 'MarkerSize', 8, ...
%                 'MarkerFaceColor', h1.Color, ...
%                 'MarkerEdgeColor', h1.Color);
%             
%             plot(t_mu_peak, mu_max, 'd', ...
%                 'MarkerSize', 8, ...
%                 'MarkerFaceColor', h2.Color, ...
%                 'MarkerEdgeColor', h2.Color);
%             
%             grid on;
%             xlabel('$t$', 'Interpreter', 'latex');
%             ylabel('Gain / Upper Bound', 'Interpreter', 'latex');
%             title('Temporal evolution of $\|\mathcal{H}\|_\infty$ and $\|\mathcal{H}_{\nabla}\|_{\mu}$', ...
%                 'Interpreter', 'latex');
%             
%             legend({ ...
%                 '$\|\mathcal{H}\|_\infty$', ...
%                 '$\|\mathcal{H}_{\nabla}\|_{\mu}$', ...
%                 '$\max\!\left(\|\mathcal{H}\|_\infty\right)$', ...
%                 '$\max\!\left(\|\mathcal{H}_{\nabla}\|_{\mu}\right)$'}, ...
%                 'Interpreter', 'latex', ...
%                 'Location', 'best');


            local_results{1, j3} = struct( ...
                'kx', kx, ...
                'kz', kz, ...
                't0', t0(j3), ...
                'ts', ts, ...
                'hinf_norm', hinf_norm, ...
                'Hinf_max', Hinf_max, ...
                't_peak', t_peak, ...
                'hinf_norm_grad_t', hinf_norm_grad_t, ...
                'mu_upper_H_inf_grad_t', mu_upper_H_inf_grad_t, ...
                'mu_upper_freq_t', mu_upper_freq_t, ...
                'mu_max', mu_max, ...
                't_mu_peak', t_mu_peak, ...
                'mu_upper_H_inf_grad_neg_t', mu_upper_H_inf_grad_neg_t, ...
                'mu_upper_freq_neg_t', mu_upper_freq_neg_t);
end

operator.local_results = local_results;  % embed optional
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


function [P,F] = sparse_method(P_option,n,bandwidth,epsP,eta)

tic;

floorP = epsP*eta;

if strcmp(P_option,'full')
    %full n*n P matrix
    P = sdpvar(n,n,'hermitian','complex');

    %constraint that P is PSD.
    F = P - floorP*eye(n) >= 0;

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

