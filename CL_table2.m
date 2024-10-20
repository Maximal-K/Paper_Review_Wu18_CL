clc;
clear;
number=zeros(1,8);
xf=zeros(1,8);
y_approximate=[];
x_c=[];
for k=1:8
     
        y_real=qfunc(k);
        for i=0:100
            if i==0 
               x_c(i+1)=sqrt(-pi./2.*log(4.*y_real));       %用来放每次xc的数值的数组
               x=x_c(i+1);                                     %数组对x进行赋值，进行后续计算
            else
                 
                x=x_c(i+1);                                   %数组对x进行赋值，进行后续计算
            end
             c_o=(sqrt(x.^4+6*x.^2+1)+x.^2+1)/4;            %用上个循环的x来计算c
              x=sqrt(-(4*c_o)/(2*c_o+1)*log(sqrt(pi/(exp(1)*c_o))*(2*c_o+1)*y_real));%用上个循环的x来计算c
            if i==0
              
              x_c(i+2)=x;
            else
                 x_i=x_c(i+1);
              x_c(i+2)=x;
                    
            end
                
                if i~=0 && abs(x-x_i)<10^-4
                                 break;
                 end
            
        end
        
    xf(k)=x;

end
