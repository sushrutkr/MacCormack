%Simulation of 1D Quasi Steady Supersonic Flow Through Nozzle
%Non-Conservative Form of governing equations
%MacCormack Method

clear all
close all
clc

%Inputs
n=31;
x=linspace(0,3,n);
dx=x(2)-x(1);
gamma=1.4;
nt=1400;    %Time Steps
error_tolerance=1e-6;

%Initial Profiles
rho=1-0.3146*x;    %Density
t=1-0.2314*x;      %temperature
v=(0.1+1.09*x).*t.^0.5;   %Velocity

a=1+2.2*(x-1.5).^2;  %area

%time step size
C=0.5; %cfl
del_t=C*dx./(sqrt(t)+v);
dt=min(del_t);

rho_initial=rho;
v_initial=v;
t_initial=t;

drhodt_p=zeros(1,n);
dvdt_p=zeros(1,n);
dtdt_p=zeros(1,n);
drhodt_c=zeros(1,n);
dvdt_c=zeros(1,n);
dtdt_c=zeros(1,n);
error=9e9;

%Time Loop
for k=1:nt
    %at beginning of time step
    rho_old=rho;
    t_old=t;
    v_old=v;   
   
    M=v./sqrt(t);
    p=rho.*t;
    
    %For Plotting variables at throat
    rho_plot(k)=rho(16);
    t_plot(k)=t(16);
    p_plot(k)=p(16);
    M_plot(k)=M(16);
    
    
   %Non-Conservative, Non- Dimensionalized form of governing equations
   
   %Continuity Equation
   %drhodt=-rho*(dvdx)-rho*v*(dlogadx)-v*(drhodx)
   
   %Momentum Equation
   %dvdy=-v*dvdx-(1/gamma)*(dtdx+t*drhodx/rho)
   
   %Energy Equation
   %dtdt=-v*dtdx-(gamma-1)*t(dvdx+v*dlogadx)
   
   for j=2:n-1 
       %Predictor Method
      dvdx=(v(j+1)-v(j))/dx;
      dlogadx=(log(a(j+1))-log(a(j)))/dx;
      drhodx=(rho(j+1)-rho(j))/dx;
      dtdx=(t(j+1)-t(j))/dx;
      
      %Forward Differencing space terms in continuity equation
      drhodt_p(j)=-rho(j)*dvdx -rho(j)*v(j)*dlogadx -v(j)*drhodx;
      
      %Forward Differencing space terms in momentum equation
      dvdt_p(j)=-v(j)*dvdx-(1/gamma)*(dtdx+t(j)*drhodx/rho(j));
      
      %Forward Differencing space terms in energy equation
      dtdt_p(j)=-v(j)*dtdx-(gamma-1)*t(j)*(dvdx+v(j)*dlogadx);
      
      
      %solution update
      rho(j)=rho(j)+drhodt_p(j)*dt;
      v(j)=v(j)+dvdt_p(j)*dt;
      t(j)=t(j)+dtdt_p(j)*dt;     
   end
   
   for j=2:n-1
   %Corrector Method
      dvdx=(v(j)-v(j-1))/dx;
      dlogadx=(log(a(j))-log(a(j-1)))/dx;
      drhodx=(rho(j)-rho(j-1))/dx;
      dtdx=(t(j)-t(j-1))/dx;
      
      %Backward Differencing space terms in continuity equation
      drhodt_c(j)=-rho(j)*dvdx -rho(j)*v(j)*dlogadx -v(j)*drhodx;
      
      %Backward Differencing space terms in momentum equation
      dvdt_c(j)=-v(j)*dvdx-(1/gamma)*(dtdx+t(j)*drhodx/rho(j));
      
      %Backward Differencing space terms in energy equation
      dtdt_c(j)=-v(j)*dtdx-(gamma-1)*t(j)*(dvdx+v(j)*dlogadx);    
   end
   
   %computing average derivative
   drhodt_av=0.5*(drhodt_p+drhodt_c);
   dvdt_av=0.5*(dvdt_p+dvdt_c);
   dtdt_av=0.5*(dtdt_p+dtdt_c);
   
   %final calculation of values at next time step
   rho=rho_old+drhodt_av*dt;
   v=v_old+dvdt_av*dt;
   t=t_old+dtdt_av*dt;
   
   %Boundary conditions
   
   %Inflow
   v(1)=2*v(2)-v(3);
   
   %Outflow
   v(n)=2*v(n-1)-v(n-2);
   rho(n)=2*rho(n-1)-rho(n-2);
   t(n)=2*t(n-1)-t(n-2);
   

   error_rho=max(abs(rho-rho_old));
   error_t=max(abs(t-t_old));
   error_v=max(abs(v-v_old));
   error=max([error_rho;error_t;error_v]);
   
   if error<error_tolerance
       break;
    end
   
end

%Convergence Plot
  figure(1)
    plot(rho_plot,'b','linewidth',1.20);
    xlabel('Time steps');
    hold on;
    plot(t_plot,'r','linewidth',1.20);
    hold on;
    plot(p_plot,'g','linewidth',1.20);
    hold on;
    plot(M_plot,'k','linewidth',1.20);
    hold on;
    title(sprintf('Variations in properties at time step %d',k));
    legend('\rho/\rho_{o}','T/T_{o}','p/p_{o}','Mach');
    pause(0.0001);

%Mach number in terms of non dimensionalized velocity and temperature
figure(2);
M=v./sqrt(t);
plot(x,M);
ylabel('\bf Mach');
xlabel('\bf Position');
title(sprintf('Time step %d',k));
grid on;


%Normalized Mass Flow Rate in terms of non dimensionalized variables
figure(3);
m=rho.*a.*v;
plot(x,m);
xlabel('\bf Position');
ylabel('\bf Mass flow rate,\rho*A*V');
axis([0 3 0.5 0.6]);
title(sprintf('Mass Flow Rate at time step %d',k));
grid on;

figure(4);
%plotting profiles    
    subplot(5,1,1)
    plot(x,a);
    ylabel('\bf A/A^{\ast}');
    title(sprintf('Non Conservative Form and time step %d',k));
  
    subplot(5,1,2)
    plot(x,rho);
    ylabel('\bf \rho/\rho_{o}');
    title('Density')
    
    subplot(5,1,3)
    plot(x,v);
    ylabel('\bf v/v_{o}');
    title('Velocity');
    
    subplot(5,1,4)
    plot(x,t);
    ylabel('\bf T/T_{o}');
    title('Temperature');
    
    subplot(5,1,5)
    plot(x,p);
    ylabel('\bf p/p_{o}');
    xlabel('\bf Position');
    title('Pressure');

