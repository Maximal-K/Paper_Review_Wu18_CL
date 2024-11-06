%% parameter
clear all;clc;

%% SEP func of exact Q func（11）
pes_eq = @(M, sigma_e) 2*qfunc(pi./(M*sigma_e));

%% SEP func of Lower Bound of Q func

% f = @(x) sqrt(x.^4+6*x.^2+1);           % 正确
 Q_our = @(x) sqrt(exp(1)/pi).*sqrt(sqrt(x.^4+6.*x.^2+1)+x.^2+1)./...
     (sqrt(x.^4+6*x.^2+1)+x.^2+3).*...
     exp(-(sqrt(x.^4+6*x.^2+1)+x.^2-1)./4);

% Q_beta = @(x) sqrt(2).*exp(-5/16.*x.^2-1/2-x.*sqrt(9.*x.^2+48)./16).*...
%     (-3.*x+sqrt(9.*x.^2+48))./(8*sqrt(pi)); 
% Q_our = @(x) sqrt(-3*sqrt(3)*pi.*log(-1/2+sqrt(1+24.*x)/2))/3;
% pes_lq = @(M, sigma_e, f) 2*Qco(pi./(M*sigma_e), f);       Q Q_c
% pes_lq = @(M, sigma_e) 2*Q_beta(pi./(M*sigma_e));          % Q_beta

pes_lq = @(M, sigma_e) 2.*Q_our(pi./(M.*sigma_e));             % Q_our

Lo = @(gamma, sigma_p) sqrt(1./2.*(1+3./(gamma.*sigma_p^2)));

sigma_e = @(gamma, L, sigma_p) sqrt((1+1./L).*1./(2.*gamma)+...
    (2*L.^2+3.*L+1)./(6.*L).*sigma_p^2);

%% function
% SEP of exact Q
L = 10;              % 固定窗口长度

%% 固定L, 8PSK
x1 = 10.^((12:25)./10);
M = 8;               % 8-PSK
sigma_p = sqrt(8.4e-4);     % 8-PSK时的噪声高斯分布方差
sigma_e_1 = sigma_e(x1, L, sigma_p);

P_EQ_L8 = pes_eq(M, sigma_e_1);
P_LQ_L8 = pes_lq(M, sigma_e_1);

% f1 = f(pi./(M*sigma_e_1));
% P_LQ_L8 = pes_lq(M, sigma_e_1, f1);

%% 最优Lo, 8PSK
Lo1 = round(Lo(x1, sigma_p));         % 最优窗口长度
sigma_e_2 = sigma_e(x1, Lo1, sigma_p);

P_EQ_Lo8 = pes_eq(M, sigma_e_2);
P_LQ_Lo8 = pes_lq(M, sigma_e_2);

% f2 = f(pi./(M*sigma_e_2));
% P_LQ_Lo8 = pes_lq(M, sigma_e_2, f2);
%% 固定L, 16PSK
x2 = 10.^((16:32)./10);
M = 16;
sigma_p = sqrt(2e-4);      % 16-PSK时的噪声高斯分布方差
sigma_e_1 = sigma_e(x2, L, sigma_p);

P_EQ_L16 = pes_eq(M, sigma_e_1);
P_LQ_L16 = pes_lq(M, sigma_e_1);

% f1 = f(pi./(M*sigma_e_1));
% P_LQ_L16 = pes_lq(M, sigma_e_1, f1);
% 最优Lo, 16PSK

Lo2 = round(Lo(x2, sigma_p));         % 最优窗口长度
sigma_e_2 = sigma_e(x2, Lo2, sigma_p);

P_EQ_Lo16 = pes_eq(M, sigma_e_2);
P_LQ_Lo16 = pes_lq(M, sigma_e_2);

% f2 = f(pi./(M*sigma_e_2));
% P_LQ_Lo16 = pes_lq(M, sigma_e_2, f2);
%% 固定L, 32PSK
x3 = 10.^((23:35)./10);
M = 32;
sigma_p = sqrt(5e-5);      % 32-PSK时的噪声高斯分布方差
sigma_e_1 = sigma_e(x3, L, sigma_p);

P_EQ_L32 = pes_eq(M, sigma_e_1);
P_LQ_L32 = pes_lq(M, sigma_e_1);

% f1 = f(pi./(M*sigma_e_1));
% P_LQ_L32 = pes_lq(M, sigma_e_1, f1);

% 最优Lo, 32PSK
Lo3 = round(Lo(x3, sigma_p));         % 最优窗口长度
sigma_e_2 = sigma_e(x3, Lo3, sigma_p);

P_EQ_Lo32 = pes_eq(M, sigma_e_2);
P_LQ_Lo32 = pes_lq(M, sigma_e_2);

% f2 = f(pi./(M*sigma_e_2));
% P_LQ_Lo32 = pes_lq(M, sigma_e_2, f2);



% 固定L, 64PSK
x4 = 10.^((29:39)./10);
M = 64;
sigma_p = sqrt(1.3e-5);      % 64-PSK时的噪声高斯分布方差
sigma_e_1 = sigma_e(x4, L, sigma_p);

P_EQ_L64 = pes_eq(M, sigma_e_1);
P_LQ_L64 = pes_lq(M, sigma_e_1);

% f1 = f(pi./(M*sigma_e_1));
% P_LQ_L64 = pes_lq(M, sigma_e_1, f1);

% 最优Lo, 64PSK
Lo4 = round(Lo(x4, sigma_p));         % 最优窗口长度
sigma_e_2 = sigma_e(x4, Lo4, sigma_p);

P_EQ_Lo64 = pes_eq(M, sigma_e_2);
P_LQ_Lo64 = pes_lq(M, sigma_e_2);

% f2 = f(pi./(M*sigma_e_2));
% P_LQ_Lo64 = pes_lq(M, sigma_e_2, f2);
%% figure
figure;
% 8-PSK
h = semilogy((12:25), P_EQ_L8, 'g', 'LineWidth', 1.2);
hold on;
h = [h, semilogy((12:25), P_LQ_L8, 'gs', 'LineWidth', 1.2, 'MarkerSize', 10, 'MarkerFaceColor','w')];

h = [h, semilogy((12:25), P_EQ_Lo8, 'r', 'LineWidth', 1.2)];
h = [h, semilogy((12:25), P_LQ_Lo8, 'ro', 'LineWidth', 1.2, 'MarkerFaceColor','r')];


% 16-PSK
semilogy((16:32), P_EQ_L16, 'g', 'LineWidth', 1.2);
semilogy((16:32), P_LQ_L16, 'gs', 'LineWidth', 1.2, 'MarkerSize', 10, 'MarkerFaceColor','w');

semilogy((16:32), P_EQ_Lo16, 'r', 'LineWidth', 1.2);
semilogy((16:32), P_LQ_Lo16, 'ro', 'LineWidth', 1.2, 'MarkerFaceColor','r');


% 32-PSK
semilogy((23:35), P_EQ_Lo32, 'r', 'LineWidth', 1.2);
semilogy((23:35), P_LQ_Lo32, 'ro', 'LineWidth', 1.2, 'MarkerFaceColor','r');

semilogy((23:35), P_EQ_L32, 'g', 'LineWidth', 1.2);
semilogy((23:35), P_LQ_L32, 'gs', 'LineWidth', 1.2, 'MarkerSize', 10, 'MarkerFaceColor','w');


% 64-PSK
semilogy((29:39), P_EQ_Lo64, 'r', 'LineWidth', 1.2);
semilogy((29:39), P_LQ_Lo64, 'ro', 'LineWidth', 1.2, 'MarkerFaceColor','r');

semilogy((29:39), P_EQ_L64, 'g', 'LineWidth', 1.2);
semilogy((29:39), P_LQ_L64, 'gs', 'LineWidth', 1.2, 'MarkerSize', 10, 'MarkerFaceColor','w');

%% config
% grid on;
xlabel('γ(dB)');
ylabel('P(e_s)');
title('Demenstration of the accuracy of the invertible Q-function bound');
t = legend(h([1,2,3,4]), ' exact Q(x), L = 10', ' Q_{our}(x), L = 10',...
    ' exact Q(x), L = L_o', ' Q_{our}(x), L = L_o');

axis([10, 40, 10^-11, 10^-1-10^-2]);
set(t,'FontSize', 20, 'box', 'off', 'FontWeight', 'bold', 'Location', 'Southwest', 'Position', ([0.157908008638902,0.22260015117158,0.174446205889122,0.231150798976016]));
set(gcf,'color', 'white');
set(gca,'FontSize',20, 'FontName', 'Times New Roman');