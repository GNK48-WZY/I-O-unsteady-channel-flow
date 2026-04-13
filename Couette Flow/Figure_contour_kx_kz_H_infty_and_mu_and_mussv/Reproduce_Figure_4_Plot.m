%% Load results
% load('mu_full_scan_results.mat');  % contains results, kx_list, kz_list, delta_list

% Pick which delta index to visualize
% ind_delta = 1;

% Preallocate mu matrix
% mu_LMI_re = zeros(length(kz_list), length(kx_list));
% mu_upper_H = zeros(length(kz_list), length(kx_list));
% mu_upper_H_neg= zeros(length(kz_list), length(kx_list));

% Extract µ values from results struct
for i_kx = 1:length(kx_list)
    for j_kz = 1:length(kz_list)
          mu_LMI_1(j_kz,i_kx) = results{i_kx,j_kz}.mu_LMI;  
          mu_upper_H(j_kz,i_kx) = results{i_kx,j_kz}.mu_upper_H_inf_grad;
          mu_upper_H_neg(j_kz,i_kx) = results{i_kx,j_kz}.mu_upper_H_inf_grad_neg_freq;
          % hinfty_norm(j_kz,i_kx) = results{i_kx,j_kz}.hinf_norm;
          % hinf_norm_grad(j_kz,i_kx) = results{i_kx,j_kz}.hinf_norm_grad;
        % or use .mu_upper_H_inf_grad if you want mussv result
    end
end

%% Contour plot
figure;
[KKX, KKZ] = meshgrid(kx_list, kz_list);

contourf(KKZ, KKX, log10(mu_LMI_1), 1000, 'LineColor','none');  % filled contour
set(gca,'XScale','log','YScale','log');  % log scales like in the paper
colormap("jet");
colorbar;
clim([1 3])
% clim([1.5 4]) % for \gamma_\infty and \gamma_{\mu, \gamma}
xlabel('$k_z$','Interpreter','latex','FontSize',16);
ylabel('$k_x$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',16)  
% title(['$\mu$ contour (delta = ' num2str(delta_list(ind_delta)) ')'],'Interpreter','latex');


figure;
[KKX, KKZ] = meshgrid(kx_list, kz_list);

contourf(KKZ, KKX, log10(mu_upper_H), 1000, 'LineColor','none');  % filled contour
set(gca,'XScale','log','YScale','log');  % log scales like in the paper
colormap(jet);
colorbar;
clim([1 3])
xlabel('$k_z$','Interpreter','latex','FontSize',16);
ylabel('$k_x$','Interpreter','latex','FontSize',16);
% title(['$\mu$ contour (delta = ' num2str(delta_list(ind_delta)) ')'],'Interpreter','latex');


% figure;
% [KKX, KKZ] = meshgrid(kx_list, kz_list);
% 
% contourf(KKZ, KKX, log10(mu_upper_H_neg), 1000, 'LineColor','none');  % filled contour
% set(gca,'XScale','log','YScale','log');  % log scales like in the paper
% colormap(jet);
% colorbar;
% clim([1 3])
% xlabel('$k_z$','Interpreter','latex');
% ylabel('$k_x$','Interpreter','latex');
% title(['$\mu$ contour (delta = ' num2str(delta_list(ind_delta)) ')'],'Interpreter','latex');
% 
% figure;
% [KKX, KKZ] = meshgrid(kx_list, kz_list);
% 
% contourf(KKZ, KKX, log10(hinfty_norm), 1000, 'LineColor','none');  % filled contour
% set(gca,'XScale','log','YScale','log');  % log scales like in the paper
% colormap(jet);
% clim([1.5 4])
% colorbar;
% xlabel('$k_z$','Interpreter','latex');
% ylabel('$k_x$','Interpreter','latex');
% 
% figure;
% [KKX, KKZ] = meshgrid(kx_list, kz_list);
% 
% contourf(KKZ, KKX, log10(hinf_norm_grad), 1000, 'LineColor','none');  % filled contour
% set(gca,'XScale','log','YScale','log');  % log scales like in the paper
% colormap(jet);
% clim([1.5 4])
% colorbar;
% xlabel('$k_z$','Interpreter','latex');
% ylabel('$k_x$','Interpreter','latex');


function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end
