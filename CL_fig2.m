%% parameter
clear all;clc;
temp = 0.52:0.01:6.5;
x = -log10(qfunc(temp));
y = qfunc(temp);        %Q函数值
e = 1e-1;         %需要的精度
c = @(x) (sqrt(x^4+6*x^2+1)+x^2+1)/4;        %最优c
val = @(c, y) sqrt(-4*c/(2*c+1)*log(sqrt(pi/(exp(1)*c))*(2*c+1)*y));          %反函数值

record_times = [];
% vf = [];
% record_val = [];

%% function
for i = 1:length(temp)
    xf = sqrt(-pi/2.*log(4*y(i)));        %下界函数求x值
    % record(end+1) = xf;
    co = c(xf);
    xn = val(co, y(i));
    % record(end+1) = xn;          %记录x的值
    n = 1;                     %记录迭代次数
    while abs(xn-xf) > e
        co = c(xn);
        xf = xn;
        xn = val(co, y(i));
        n = n + 1;
    end
    record_times(i) = n;
end
plot(x, record_times,'k','LineWidth',1.2);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 1e-2;         %需要的精度
for i = 1:length(temp)
    xf = sqrt(-pi/2.*log(4*y(i)));        %下界函数求x值
    % record(end+1) = xf;
    co = c(xf);
    xn = val(co, y(i));
    % record(end+1) = xn;          %记录x的值
    n = 1;                     %记录迭代次数
    while abs(xn-xf) > e
        co = c(xn);
        xf = xn;
        xn = val(co, y(i));
    %     record(end+1) = xn;
        n = n + 1;
    end
%     vf(i) = xn;
    record_times(i) = n;
end
plot(x, record_times,'r','LineWidth',1.2);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 1e-7;         %需要的精度
for i = 1:length(temp)
    xf = sqrt(-pi/2.*log(4*y(i)));        %下界函数求x值
    % record(end+1) = xf;
    co = c(xf);
    xn = val(co, y(i));
    % record(end+1) = xn;          %记录x的值
    n = 1;                     %记录迭代次数
    while abs(xn-xf) > e
        co = c(xn);
        xf = xn;
        xn = val(co, y(i));
    %     record(end+1) = xn;
        n = n + 1;
    end
%     vf(i) = xn;
    record_times(i) = n;
end
plot(x, record_times,'b','LineWidth',1.2);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 1e-12;         %需要的精度
for i = 1:length(temp)
    xf = sqrt(-pi/2.*log(4*y(i)));        %下界函数求x值
    % record(end+1) = xf;
    co = c(xf);
    xn = val(co, y(i));
    % record(end+1) = xn;          %记录x的值
    n = 1;                     %记录迭代次数
    while abs(xn-xf) > e
        co = c(xn);
        xf = xn;
        xn = val(co, y(i));
    %     record(end+1) = xn;
        n = n + 1;
    end
%     vf(i) = xn;
    record_times(i) = n;
end
plot(x, record_times,'g','LineWidth',1.2);
hold on;

%% config
t = legend(' e = 10^{-1}',' e = 10^{-2}',' e = 10^{-7}',' e = 10^{-12}');
xlabel('-log_{10}Q(x)');
ylabel('iterations');
axis([0.52,10,0,7]);

set(t,'FontName', 'Times New Roman', 'FontSize', 16, 'fontweight','bold');
set(gcf,'color', 'white');
set(gca,'FontName', 'Times New Roman', 'FontSize',16, 'fontweight','bold');