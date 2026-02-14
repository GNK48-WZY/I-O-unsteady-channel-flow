% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))

% clear; clc;
% FD parameters
% test
N  = 18;       % total FD grid points including boundaries
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
delete(gcp('nocreate'));
% parpool(8);
[operator,local_results] = build_operators_fd_and_LMI(N, L, Re, kx, kz, t0, P_option, bandwidth);

%% Main loop
% for i_kx = 1:length(kx_list)
%     [local_results] = Inside_loop(N, L, Re, kz_list, i_kx, kx_list);
%     results(i_kx,:) = local_results;
% end

disp('successful')
save('mu_full_scan_results.mat')

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
    operator.B = LHS \ Bmat;
    operator.B_tilde = Bmat;

    % C operator
    Cmat = 1/k2 * [ 1i*kx*D   -1i*kz*eye(Ny);
        k2*eye(Ny) zeros(Ny);
        1i*kz*D   1i*kx*eye(Ny) ];
    operator.C = Cmat;

    % build gradient-output operators
    Grad_block = [1i*kx*eye(Ny); D; 1i*kz*eye(Ny)];
    Grad_big   = blkdiag(Grad_block, Grad_block, Grad_block);

    C_grad = Grad_big * operator.C;  % (9*Ny) × (2*Ny)
    Ny = size(operator.B,2)/3;
    n  = size(A,1);
%     n = size(operator.C, 2); % = 2*Ny
    operator.C_grad_u = C_grad(1:3*Ny, 1:n);
    operator.C_grad_v = C_grad(3*Ny+1:6*Ny, 1:n);
    operator.C_grad_w = C_grad(6*Ny+1:9*Ny, 1:n);
%     Ny = size(operator.B,2)/3;
%     operator.E=  LHS;

    %% LMI computation
    yalmip('clear');

    type = "time_varying_3";
    switch type
        %% with Gradient C and sx,sy,sz = sdpvar(1,1)
        case "time_varying"
            yalmip('clear');
            F = [];
            % time-varying
            for j1=1:t_steps
                %% Sparse
%                 [~,F_local] = sparse_method(P_option, n, bandwidth);
%                 F =[F,F_local];
                %% Full
                P{j1}=sdpvar(n,n,'hermitian','complex');
                F = [F, P{j1} >= 0.1*eye(n)];

                sx{j1} = sdpvar(1,1); sy{j1} = sdpvar(1,1); sz{j1} = sdpvar(1,1);

            end
            tt = 0.01;
            I = eye(n);
            scaling = 1;
            E=I*scaling;
            Gamma = sdpvar(1,1);


            for j1 = 1:t_steps
                % Always include positivity constraints for every time step
                F = [F, ...
                    sx{j1} >= 0, ...
                    sy{j1} >= 0, ...
                    sz{j1} >= 0];

                % Only form and include dV_ineq if j1 < t_steps
                if j1 < t_steps
                    M = blkdiag(sx{j1}*eye(Ny), sy{j1}*eye(Ny), sz{j1}*eye(Ny));
                    dP_dt = (P{j1 + 1} - P{j1}) / dt;
%                     dV_ineq = [ ...
%                         dP_dt + A(:,:,j1)' * P{j1} * operator.E + operator.E' * P{j1} * A(:,:,j1) ...
%                         + sx{j1} * operator.C_grad_u' * operator.C_grad_u ...
%                         + sy{j1} * operator.C_grad_v' * operator.C_grad_v ...
%                         + sz{j1} * operator.C_grad_w' * operator.C_grad_w, ...
%                         operator.E' * P{j1} * operator.B;
%                         operator.B' * P{j1}* operator.E, -Gamma * M ...
%                         ];
                    dV_ineq = [ ...
                        E'*dP_dt*E + (E * A(:,:,j1))' * P{j1} * E + E' * P{j1} * E* A(:,:,j1) ...
                        + sx{j1} * operator.C_grad_u' * operator.C_grad_u ...
                        + sy{j1} * operator.C_grad_v' * operator.C_grad_v ...
                        + sz{j1} * operator.C_grad_w' * operator.C_grad_w, ...
                        E' * P{j1} * E * operator.B;
                        (E * operator.B)' * P{j1}* E, -Gamma * M ...
                        ];


                    % Add inequality constraint
                    F = [F,dV_ineq <= 0];
                end
            end
            F=[F];
%             dV_ineq_f = [ A(:,:,end)' * P{end} * operator.E + operator.E' * P{end} * A(:,:,end) ...
%                         + sx{end} * operator.C_grad_u' * operator.C_grad_u ...
%                         + sy{end} * operator.C_grad_v' * operator.C_grad_v ...
%                         + sz{end} * operator.C_grad_w' * operator.C_grad_w, ...
%                         operator.E' * P{end} * operator.B;
%                         operator.B' * P{end}* operator.E, -Gamma * M ...
%                         ];
%             F=[F, dV_ineq_f <=0];

            sdp_options = sdpsettings('solver','mosek','verbose',1);
            %     sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            diagnostics = bisection(F,Gamma,sdp_options);



            %% with Gradient C and sx,sy,sz = 1
        case "time_varying_2"

            F = [];
            for j1=1:t_steps
                [P{j1},F_local] = sparse_method(P_option, n, bandwidth);
                F =[F,F_local];
            end
            sx = 1; sy = 1; sz = 1;
            %    gamma_H_inf_complex2 = sdpvar(1,1);
            %    tt = sdpvar(1);
            tt = 0.01;
            I = eye(n);
            scaling = Re;
            E=I*scaling;
            Gamma = sdpvar(1,1);

            for j1 = 1:t_steps
                % Always include positivity constraints for every time step
                F = [F];

                % Only form and include dV_ineq if j1 < t_steps
                if j1 < t_steps
                    M = blkdiag(sx*eye(Ny), sy*eye(Ny), sz*eye(Ny));
                    dP_dt = (P{j1 + 1} - P{j1}) / dt;
                    dV_ineq = [ ...
                        dP_dt + A(:,:,j1)' * P{j1} * operator.E + operator.E' * P{j1} * A(:,:,j1) ...
                        + sx * operator.C_grad_u' * operator.C_grad_u ...
                        + sy * operator.C_grad_v' * operator.C_grad_v ...
                        + sz * operator.C_grad_w' * operator.C_grad_w, ...
                        operator.E' * P{j1} * operator.B;
                        operator.B' * P{j1}* operator.E, -Gamma * M ...
                        ];

                    % Add inequality constraint
                    F = [F,dV_ineq <= 0];
                end
            end

            sdp_options = sdpsettings('solver','mosek','verbose',1);
            %     sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            diagnostics = bisection(F,Gamma,sdp_options);

            %% Without Gradient C(H_infty norm)
        case "time_varying_3"
            yalmip('clear');
            F = [];
           
            for j1=1:t_steps
                %                 [P{j1}, F_local] = sparse_method(P_option, n, bandwidth);
                P{j1}=sdpvar(n,n,'hermitian','complex');
                F = [F, P{j1} >= 0.1*eye(n)];
            end
            Gamma = sdpvar(1,1);
            I = eye(n);
            scaling = 1;
            E=I*scaling;
            for j1 = 1:t_steps
                if j1 < t_steps
                    M = blkdiag(eye(Ny), eye(Ny), eye(Ny));
                    dP_dt = (P{j1 + 1} - P{j1}) / dt;

%                     dV_ineq = [ dP_dt + A(:,:,j1)' * P{j1} * operator.E + operator.E' * P{j1} * A(:,:,j1) ...
%                         + operator.C' * operator.C, ...
%                        operator.E' * P{j1} * operator.B;
%                         operator.B' * P{j1} * operator.E, -Gamma * M ];
                dV_ineq = [ E'*dP_dt*E + (E * A(:,:,j1))' * P{j1}*E + E'*P{j1} * E* A(:,:,j1) ...
                        + operator.C' * operator.C, ...
                        E'*P{j1} * E* operator.B;
                        (E * operator.B)' * P{j1} * E, -Gamma * M ];
                    % Add inequality constraint
                    F = [F, dV_ineq <= 0];
                end
            end
            F=[F];
            sdp_options = sdpsettings('solver','mosek','verbose',1);
            %     sdp_options = sdpsettings('solver','mosek','cachesolvers',1);
            diagnostics = optimize(F,Gamma,sdp_options);


            %% Time-independent: Without Gradient C (H_infty norm)
            case "time_varying_4"   

            folder_name = 'Hinf_Results';
            if ~exist(folder_name, 'dir')
                mkdir(folder_name);
            end
            
            % Initialize storage for H-inf norms
            hinf_norm = zeros(1, t_steps);
            
            figure('Visible', 'off'); % Create a hidden figure for saving plots
%             ttt = (0:t_steps-1) * dt;
            w_range = {1e-4, 10^(-0.3)};
            for j1 = 1:t_steps    
%                 A0 = -1i * L_value{1};
%                 sys = ss(A0, operator.B, operator.C, zeros(size(operator.C,1)));
                A = -1i * L_value{j1};
                sys = ss(A, operator.B, operator.C, zeros(size(operator.C,1)));
                hinf_norm(j1) = norm(sys, inf);
                [sv, w] = sigma(sys, w_range); 
    
%                 gain_dB = 20*log10(sv(1,:));
%                 gain_dB = 20*log10( max(sv,[],1) );
                gain_dB = max(sv,[],1);

                
                hFig = figure('Visible', 'off'); 
                semilogx(w, gain_dB, 'b-', 'LineWidth', 1.5);
                
                grid on;
                xlabel('Frequency (rad/s)');
                ylabel('Gain (dB)');
                title(['Resolvent Spectrum at Step ', num2str(j1)]);
                
                drawnow; 
                
                file_path = fullfile(folder_name, sprintf('Sigma_step_%03d.png', j1));
                saveas(hFig, file_path);
                close(hFig);
                
                if mod(j1, 10) == 0, fprintf('Saved step %d\n', j1); end
            end
            ts = linspace(0, (t_steps-1)*dt, t_steps);
            figure('Name', 'Final H-infinity Trend');
            plot(ts, hinf_norm, '-b', 'LineWidth', 2);
            xlabel('t');
            ylabel('gain');
            title('Temporal Evolution of H_\{infty} norm');
            grid on;
            
            saveas(gcf, fullfile(folder_name, 'Hinf_Final_Plot.png'));
   end

    diagnostics.info

%     if diagnostics.problem == 0
        operator.mu_LMI = value(Gamma);
        operator.P_optimal = value(P);
%         operator.t_optimal = value(tt);
        operator.P = P;
        %         operator.sx =sx;
        %         operator.sy = sy;
        %         operator.sz = sz;
%         sx_numeric = value(sx);
%         sy_numeric = value(sy);
%         sz_numeric = value(sz);
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


%     end
%     local_results{1, j3} = struct( ...
%         'kx', kx, ...
%         'kz', kz, ...
%         'mu_LMI', operator.mu_LMI, ...
%         'P', P_cell, ...
%         'sx', sx_numeric, ...
%         'sy', sy_numeric, ...
%         'sz', sz_numeric);
    local_results{1, j3} = struct( ...
        'kx', kx, ...
        'kz', kz, ...
        'mu_LMI', operator.mu_LMI, ...
        'P', P_cell);

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
