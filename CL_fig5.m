 clc;
clear;
% qfuncinv
   L=10;
   x_val=linspace(10.^(-6),10.^(-3),1000);
   x_c=[];
  y_real=[];
 for k=1:1
        y_real(k)=10.^(-6)./2;
        for i=0:20
            if i==0 
               x_c(i+1)=sqrt(-pi./2.*log(4.*y_real(k)));       %用来放每次xc的数值的数组
               x=x_c(i+1);                                     %数组对x进行赋值，进行后续计算
            else        
                x=x_c(i+1);                                   %数组对x进行赋值，进行后续计算
            end
             c_o=(sqrt(x.^4+6*x.^2+1)+x.^2+1)/4 ;          %用上个循环的x来计算c
              x=sqrt(-(4*c_o)/(2*c_o+1)*log(sqrt(pi/(exp(1)*c_o))*(2*c_o+1)*y_real(k)));%用上个循环的x来计算c
            if i==0
              x_c(i+2)=x;
            else
                 x_i=x_c(i+1);
              x_c(i+2)=x;          
            end
                if i~=0 && abs(x-x_i)<10^-3
                                 break;
                end
        end        
    xf(k)=x;
 end
Sigma_e_c2=@(M)pi.^2./M.^2.*(x).^(-2);

%% SEP search
vT_M8=linspace(10.^(-5),10.^(-3),1000);                         %M=8时的变量
Sigma_p=@(vT)sqrt(2.*pi.*vT);                                   %横坐标变量
Sigma_e_Exact2=@(M)pi.^2./(M.^2).*(qfuncinv(10.^(-6)./2)).^(-2);%Sigma_e的值用Q反函数来算
gamma_1=@(sigma_p,sigma_e)3./(sigma_p.^2-...                    %gamma1的式子
    4.*sqrt(3).*sigma_e.*sigma_p+6.*sigma_e.^2);
gamma_c=@(M)10*log10((qfuncinv(10.^(-6)./2)./sin(pi./M)).^2./2);%gammac的式子
gamma_dB_p_M8=10.*log10(gamma_1(Sigma_p(vT_M8),sqrt(Sigma_e_Exact2(8))))-...
    gamma_c(8);                                                 %所画曲线的表达式
gamma1_dB_p_M8=10.*log10(gamma_1(Sigma_p(vT_M8),sqrt(Sigma_e_c2(8))))-gamma_c(8);
L=@(gamma,sigma_p) sqrt((1+3./(gamma.*sigma_p.^2))./2);
gamma_2=@(l,sigma_e,sigma_p) (1+1./l)./(2.*sigma_e-...
    (2.*l.^2+3.*l+1)./3./l.*sigma_p.^2);
gamma2_dB_p_M8=10.*log10...
    (gamma_2(L(gamma_1(Sigma_p(vT_M8),sqrt(Sigma_e_c2(8))),Sigma_p(vT_M8)),...
    Sigma_e_c2(8),Sigma_p(vT_M8)))-gamma_c(8);

figure
semilogx(vT_M8,gamma_dB_p_M8, 'r-x','MarkerIndices',...
    1:75:length(vT_M8));
hold on
semilogx(vT_M8,gamma1_dB_p_M8, 'go','MarkerIndices',...
    1:75:length(vT_M8));
semilogx(vT_M8,gamma2_dB_p_M8, 'bs','MarkerIndices',...
    1:75:length(vT_M8));

%%
vT_M16=linspace(10.^(-5),2.5.*10.^(-4),1000);
gamma_dB_p_M16=10.*log10...
    (gamma_1(Sigma_p(vT_M16),sqrt(Sigma_e_Exact2(16))))-gamma_c(16);
semilogx(vT_M16,(gamma_dB_p_M16), 'r-x','MarkerIndices',...
    1:75:length(vT_M16));
gamma1_dB_p_M16=10.*log10...
    (gamma_1(Sigma_p(vT_M16),sqrt(Sigma_e_c2(16))))-gamma_c(16);
semilogx(vT_M16,gamma1_dB_p_M16, 'go','MarkerIndices',...
    1:75:length(vT_M16));
gamma2_dB_p_M16=10.*log10...
    (gamma_2(L(gamma_1(Sigma_p(vT_M16),sqrt(Sigma_e_c2(16))),Sigma_p(vT_M16)),...
    Sigma_e_c2(16),Sigma_p(vT_M16)))-gamma_c(16);
semilogx(vT_M16,gamma2_dB_p_M16, 'bs','MarkerIndices',...
    1:75:length(vT_M16));

%%
vT_M32=linspace(3.*10.^(-6),6.5.*10.^(-5),1000);
gamma_dB_p_M32=10.*log10...
    (gamma_1(Sigma_p(vT_M32),sqrt(Sigma_e_Exact2(32))))-gamma_c(32);
semilogx(vT_M32,gamma_dB_p_M32, 'r-x','MarkerIndices',...
    1:75:length(vT_M32));
gamma1_dB_p_M32=10.*log10...
    (gamma_1(Sigma_p(vT_M32),sqrt(Sigma_e_c2(32))))-gamma_c(32);
semilogx(vT_M32,gamma1_dB_p_M32, 'go','MarkerIndices',...
    1:75:length(vT_M32));
gamma2_dB_p_M32=10.*log10...
    (gamma_2(L(gamma_1(Sigma_p(vT_M32),sqrt(Sigma_e_c2(32))),Sigma_p(vT_M32)),...
    Sigma_e_c2(32),Sigma_p(vT_M32)))-gamma_c(32);
semilogx(vT_M32,gamma2_dB_p_M32, 'bs','MarkerIndices',...
    1:75:length(vT_M32));

%%
vT_M64=linspace(2.*10.^(-6),1.6.*10.^(-5),1000);
gamma_dB_p_M64=10.*log10...
    (gamma_1(Sigma_p(vT_M64),sqrt(Sigma_e_Exact2(64))))-gamma_c(64);
semilogx(vT_M64,gamma_dB_p_M64, 'r-x','MarkerIndices',...
    1:75:length(vT_M64));
gamma1_dB_p_M64=10.*log10...
    (gamma_1(Sigma_p(vT_M64),sqrt(Sigma_e_c2(64))))-gamma_c(64);
semilogx(vT_M64,gamma1_dB_p_M64, 'go','MarkerIndices',...
    1:75:length(vT_M64));
gamma2_dB_p_M64=10.*log10...
    (gamma_2(L(gamma_1(Sigma_p(vT_M64),sqrt(Sigma_e_c2(64))),Sigma_p(vT_M64)),...
    Sigma_e_c2(64),Sigma_p(vT_M64)))-gamma_c(64);
semilogx(vT_M64,gamma2_dB_p_M64, 'bs','MarkerIndices',...
    1:75:length(vT_M64));

axis([10^-6, 2.*10^-3,0, 16 ]);
xlabel('\Delta\nuT',"FontName","Times New Roman")
ylabel('\gamma_{p}(dB)',"FontName","Times New Roman")

figure
semilogy(gamma_dB_p_M8,vT_M8,'r-x', gamma_dB_p_M16,vT_M16,'r-x',...
    gamma_dB_p_M32,vT_M32, 'r-x',gamma_dB_p_M64,vT_M64, 'r-x','MarkerIndices',...
    1:75:length(vT_M8));
hold on
semilogy(gamma1_dB_p_M8,vT_M8,'g.', gamma1_dB_p_M16,vT_M16,'g.',...
    gamma1_dB_p_M32,vT_M32,'g.',gamma1_dB_p_M64,vT_M64,'g.','MarkerIndices',...
    1:75:length(vT_M8),'MarkerSize',12);
semilogy(gamma2_dB_p_M8,vT_M8,'bs', gamma2_dB_p_M16,vT_M16,'bs',...
    gamma2_dB_p_M32,vT_M32,'bs',gamma2_dB_p_M64,vT_M64,'bs','MarkerIndices',...
    1:75:length(vT_M8));
axis([0 8 10^(-6) 10^(-3)]);
