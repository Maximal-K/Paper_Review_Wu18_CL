%% 参数
clear 
clc
close
db=linspace(0,10,30);
t=10.^(db./20);
%% 实际的q函数
exactQ=qfunc(t);
semilogy(db,exactQ,'k-','Linewidth',1.5)
hold
%% Qlbkw自己的上下界
Co = @(t)(sqrt(t.^4+6.*t.^2+1)+t.^2+1)./4;
f_x = @(t)sqrt(t.^4+6.*t.^2+1);
Q_co = @(t) sqrt(exp(1)./pi).*sqrt(f_x(t)+t.^2+1)./(f_x(t)+t.^2+3)...
    .*exp((f_x(t)+t.^2-1)./(-4));
Q_c = @(t, given)sqrt(exp(1)./pi).*sqrt(Co(given))./(2.*Co(given)+1)...
    .*exp((2.*Co(given)+1)./(-4.*Co(given)).*t.^2);
semilogy(db,Q_c(t,10.^(0./20)),'b-','Linewidth',1.5)
semilogy(db,Q_c(t,10.^(10./20)),'g-','Linewidth',1.5)
semilogy(db,Q_co(t),'r-','Linewidth',1.5)
axis([0 10 10.^(-4) 0.2])
legend( 'boxoff')
legend('ExactQ','Q_{c}(x) optimized at x = 0 dB',...
    'Q_{c}(x) optimized at x = 10 dB','Q_{Co}(x)',...
    'Location','southwest','FontSize',16,'FontName','Times New Roman')
%% 画画
% imshow(uint8(null),'border','tight','initialmagnification','fit')
% semilogy(db,Qlbkw1,'-k')
% hold on
% semilogy(db,Qlbkw2,'.-k')
% semilogy(db,Qlbkw3,'k-+')
% semilogy(db,exactQ,'r')
% imshow(uint8(data),'border','tight','initialmagnification','fit')
xlabel('x(dB)=20*log(10)x')
ylabel('Q(x)')
% legend('Q_{LB-KW-1}_(x)','Q_{LB-KW-2}_(x)','Q_{LB-KW-3}_(x)',...
%     'Q_{LB-CDS-3}_(x)','Q_{LB-A-3}_(x)','Q_{UB-CDS-3}_(x)',...
%     'exactQ','Location','southwest','FontSize',16)
% axis([-6 5 0.03 0.3])







