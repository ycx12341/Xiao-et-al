% mean_var_ref.m
% Author: Yunchen Xiao
% This MATLAB file generates the reference summary statistics based on the
% reference dataset simulated.

% Environment settings
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

% Space
l1=0;
l2=1;
x11=linspace(l1,l2,80);
h=abs(x11(2)-x11(1));
 
% Timec
T=10;
dt = 0.001;
time = 0:dt:T;
tv=0;
 
% Parameters
dn = 0.01;
gamma = 0.05;
ita = 10;
dm = 0.01;
alpha = 0.1;
r = 5;

beta = 0;
eps = 0.01;

% Initial condition
n0 = repelem(0,length(x11));

for i = 1:length(x11)
    if x11(i)<=0.25
        n0(i)=exp(-(x11(i)^2)/eps);
    else 
        n0(i)=0;
    end 
end 

n = n0;

f0 = 1-0.5*n0;
f=f0;

m0=0.5*n0; 
m = m0;

%%%Inits.matrix%%%
inits = [n;f;m];
%%%%%%%%%%%%%%%%%%

%%%%% Initial plot %%%%%%%%%%%%%%%%%%%
%figure
%plot(x11,n,'b-',x11,f,'k-',x11,m,'g-')
%axis([min(x11) max(x11) 0 max(f)+0.1])
%axis square
%xlabel('$x_1$')
%title(['$System$ at $t=$',num2str(0)])
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = [];

p=1;

n=n0;
m=m0;
f=f0;

%%% Results matrix%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
while p*dt<=T
    f(2:length(x11)-1) = -ita*dt*m(2:length(x11)-1).*f(2:length(x11)-1)+f(2:length(x11)-1);
    m(2:length(x11)-1) = dm*(m(1:length(x11)-2)+m(3:length(x11))-2*m(2:length(x11)-1))*dt/(h^2)+alpha*n(2:length(x11)-1)*dt-beta*m(2:length(x11)-1)*dt+m(2:length(x11)-1);
    n(2:length(x11)-1) = dn*(n(1:length(x11)-2)+n(3:length(x11))-2*n(2:length(x11)-1))*dt/(h^2)...
                         -gamma*(n(3:length(x11))-n(2:length(x11)-1)).*(f(3:length(x11))-f(2:length(x11)-1))*dt/(h^2)...
                         -gamma*n(2:length(x11)-1).*(f(1:length(x11)-2)+f(3:length(x11))...
                         -2*f(2:length(x11)-1))*dt/(h^2)+n(2:length(x11)-1)+r*(1-f(2:length(x11)-1)-n(2:length(x11)-1)).*n(2:length(x11)-1)*dt;                     
    
    %Homogeneous & Zero Neumann
    n(1) = n(2);
    n(length(x11)) = n(length(x11)-1);
    
    f(1) = f(2);
    f(length(x11)) = f(length(x11)-1);
    
    m(1) = m(2);
    m(length(x11)) = m(length(x11)-1); 

     
     
     
    %% Plots
    if mod(p,1000)==0
        
        temp = time(p+1);
        
        disp(temp);
     
        figure
        plot(x11,n,'b-',x11,f,'r-',x11,m,'g-')
        axis([min(x11) max(x11) 0 1.5])
        axis square
        xlabel('$x$')
        title(['$t=$ ',num2str(round(time(p+1)*100)/100)])
        
        res = [res;n;f;m];
        
        %p;
        %tv=[tv p*dt];
        %clf
        %plot(x11,m,'g-')
        %axis([min(x11) max(x11) 0 5])
        %axis square
        %xlabel('$x_1$')
        %title(['Pattern at $t=$',num2str(round(time(p+1)*100)/100),' $r_{n} = 0$'])
        %drawnow
    end
     
    p=p+1;
     
end

%figure
%subplot(2,3,1)
%plot(x11,res(13,:),'b-',x11,res(14,:),'r-',x11,res(15,:),'g-')
%axis([min(x11) max(x11) 0 1])
%axis square
%set(gca,'Xticklabel',[])
%xlabel('$x$')
%title(['Invasion pattern at $t=$ ',num2str(1)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,3,2)
%plot(x11,res(28,:),'b-',x11,res(29,:),'r-',x11,res(30,:),'g-')
%axis square
%set(gca,'Xticklabel',[], 'Yticklabel', [])
%xlabel('$x$')
%title(['$t=$ ',num2str(2)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,3,3)
%plot(x11,res(43,:),'b-',x11,res(44,:),'r-',x11,res(45,:),'g-')
%axis square
%set(gca,'Xticklabel',[], 'Yticklabel', [])
%xlabel('$x$')
%title(['$t=$ ',num2str(3)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,3,4)
%plot(x11,res(58,:),'b-',x11,res(59,:),'r-',x11,res(60,:),'g-')
%axis([min(x11) max(x11) 0 1])
%axis square
%set(gca,'XTick',[], 'YTick', [])
%xlabel('$x$')
%title(['$t=$ ',num2str(4)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,3,5)
%plot(x11,res(73,:),'b-',x11,res(74,:),'r-',x11,res(75,:),'g-')
%axis([min(x11) max(x11) 0 1])
%axis square
%set(gca,'Yticklabel', [])
%xlabel('$x$')
%title(['$t=$ ',num2str(5)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,3,6)
%plot(x11,res(88,:),'b-',x11,res(89,:),'r-',x11,res(90,:),'g-')
%axis([min(x11) max(x11) 0 1])
%axis square
%set(gca,'Yticklabel', [])
%xlabel('$x$')
%title(['$t=$ ',num2str(6)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_arr = zeros(30,80);

res_arr(1:10,:) = res(1:3:28,:);
res_arr(11:20,:) = res(2:3:29,:);
res_arr(21:30,:) = res(3:3:30,:); % rearrange the result matrix, 1-5, TC, 6-10, ECM, 11-15, MDE.

mean_var = zeros(240,2); 


mean_var(1:80,1) = mean(res_arr(1:10,:)); % mean of the TC time series
mean_var(1:80,2) = var(res_arr(1:10,:)); % variance of the TC time series

mean_var(81:160,1) = mean(res_arr(11:20,:)); % mean of the ECM time series
mean_var(81:160,2) = var(res_arr(11:20,:)); % variance of the ECM time series

mean_var(161:240,1) = mean(res_arr(21:30,:)); % mean of the MDE time series
mean_var(161:240,2) = var(res_arr(21:30,:)); % variance of the MDE time series

% Write the means and variances into the file "mean_var_obs.txt"
fileID = fopen('mean_var_obs.txt','w');
fprintf(fileID,'%5.4f %5.4f\r\n',mean_var');
fclose(fileID);