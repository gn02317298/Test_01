% Lattice Boltzmann Method 1-D combine Neutron Transport and T eq. 
clc;clear;
close all
% set(0,'defaultaxesfontname','times');
% set(0,'defaultaxesfontsize',20);


%============== figure settings ==========================
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',16);
% Defaults for this blog post
width = 6;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
% set(figure,'visible','off')
set(figure,'visible','on')
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(gcf,'position',[80 100 800 600])  % 
%============== figure settings ==========================
tic;                            % begin timer

NGrid = 11;
% NGrid = 101; 
% NGrid = 201;                                                              
% NGrid = 801;                                                            
% NGrid = 1201;  
% NGrid = 1601;
n = NGrid;
Q = 10;
S = -Q ;
% L = 60;
% L = 20;
L = 0.2;
% L = 1;
dx = L/(n-1);

x_exact = linspace(0,L,200);
x_Num = linspace(0,L,n);


T_Num = zeros(n,1);              % Initial condition for Temp
pre_T_Num = zeros(n,1);
 
err = 1;               % give initial value to activate the while loop
e1T = 1E-12;
iter=0;
% S = -1E7;
while (err > e1T)
% while (err2g > e1T)
    
      iter=iter+1;
      %calculate dimensionless temperature theda
      pre_T_Num = T_Num;
      
      % Method 1
      %phi_0_Num(2:end-1) = ( phi_0_Num(1:end-2) + phi_0_Num(3:end) - S*dx^2)/(2+dx^2/L^2);

        
      % Method 2
      for i = 2:n-1
          T_Num(i) = 0.5*( T_Num(i-1) + T_Num(i+1) - S*dx^2);
      end
        
      T_Num(1)   = 0;   % BC 1
      T_Num(end) = 0; % BC 2
      
      err = norm( pre_T_Num - T_Num )/norm(T_Num);
 
     if mod(iter,10000) == 0
        clc
        fprintf('N:\t\t\t%i\n',NGrid)
        fprintf('Iterations:\t%i\n',iter)
        fprintf('Error Temp:  \t%5.3E\n',err)

     end  
    
     if iter>=2E7
        break
     end    
  
end
%=====================solution Finite Difference with relaxation Mehtond===
toc    % end of timer
%--------------- P-1 With Marshark BC, exact sol ---------------------
T_exact = 1/2*Q*x_exact.*(L-x_exact);
% phi_0_Marshark_exact = L^2*S*(-1+cosh( (W-2*x_exact)/(2*L))*sech(W/(2*L)));
%--------------- P-1 With Marshark BC, exact sol ---------------------

figure(1)
set(gcf,'position',[80 100 800 600])
plot(x_Num, T_Num,'-'); hold on
plot(x_exact,T_exact,'--')
legend('FDM', 'Exact')
xlabel('x (non-dimensional)')
ylabel('Temperature (non-dimensional)')
title([' TL = 1.0 , TR = 0.5, heat source = 1/N*\theta^4(x) ' 10 ' N = 0.125 , NGrid = ' num2str(n)])

