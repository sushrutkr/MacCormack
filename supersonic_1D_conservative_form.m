%Simulation of 1D Quasi Steady Supersonic Flow Through Nozzle
%Conservative Form of governing equations
%MacCormock Method

clear all
close all
clc

%Inputs
n=31;
x=linspace(0,3,n);
dx=x(2)-x(1);
gamma=1.4;
nt=1400;
error_tolerance=1e-6;

%Initial Profiles
a=1+2.2*(x-1.5).^2;


%Density and Temperature

rho=zeros(1,n);
t=zeros(1,n);

%Initial Conditions
for i=1:n
   rho(i)=1;
   t(i)=1;
   if x(i)==0.5
      break;
   end
end

for j=i:n
   rho(j)=1-0.366*(x(j)-0.5);
   t(j)=1-0.167*(x(j)-0.5);
    if x(j)==1.5
        break;
    end
end

for k=j:n
    rho(k)=0.634-0.3879*(x(k)-1.5);
    t(k)=0.833-0.3507*(x(k)-1.5);
    if x(k)==3.5
        break;
    end
end

%Velocity
v=((rho.*a)/0.59).^-1;
e=t;
%Initial Conditions for solution vectors
u1=rho.*a;
u2=rho.*a.*v;
u3=rho.*((e/(gamma-1))+(gamma/2)*v.^2).*a;

%timestep=min(C*dx/vel)
%for cfl number 0.5
dt=0.0267;
%time step size
% C=0.5; %cfl
% del_t=C*dx./(sqrt(t)+v);
% dt=min(del_t);



f1=zeros(1,n);
f2=zeros(1,n);
f3=zeros(1,n);
j2=zeros(1,n);
du1dt_p=zeros(1,n);
du2dt_p=zeros(1,n);
du3dt_p=zeros(1,n);
du1dt_c=zeros(1,n);
du2dt_c=zeros(1,n);
du3dt_c=zeros(1,n);


u1_initial=u1;
u2_initial=u2;
u3_initial=u3;

%Dependent Variables
%u1=rho*a
%u2=rho*a*v
%u3=rho*((e/(gamma-1))+(gamma/2)*v^2)*a

%e=t as non dimensionalized

%Time Loop
for l=1:nt
    %Non dimensionalized pressure and mach number
    p=rho.*t;
    M=v./sqrt(t);
    
    %for plotting variables at throat
    rho_plot(l)=rho(16);
    t_plot(l)=t(16);
    p_plot(l)=p(16);
    M_plot(l)=M(16);
    
    
    %Working with dependent variables u1,u2 and u3 instead of primitive variables rho,t,v
    u1_old=u1;
    u2_old=u2;
    u3_old=u3;
    
    %Computing the flux and source terms in pure form
    f1=u2;
    f2=((u2.^2)./u1)+((gamma-1)/gamma)*(u3-(gamma*u2.^2)./(2*u1));
    f3=(gamma.*u2.*u3./u1)-((gamma*(gamma-1)*u2.^3)./(2*u1.^2));
    
    f1_old=f1;
    f2_old=f2;
    f3_old=f3;
       
    %Predictor Step
    for m=2:n-1    
       %Modelling the governing equation
       %du1dt=-df1dx
       %du2dt=-df2dx+j2
       %du3dt=-df3dx
       
       %Forward Differncing
       df1dx=(f1(m+1)-f1(m))/dx;
       df2dx=(f2(m+1)-f2(m))/dx;
       df3dx=(f3(m+1)-f3(m))/dx;
       dadx=(a(m+1)-a(m))/dx;
       j2=rho(m)*t(m)*dadx/gamma;
       
       %Governing Equations
       du1dt_p(m)=-df1dx;        %Continuity
       du2dt_p(m)=-df2dx+j2;     %Momentum
       du3dt_p(m)=-df3dx;        %Energy
       
       %Solution Update
       %Predicted Values of dependent variables/solution vectors
       u1(m)=u1(m)+du1dt_p(m)*dt;
       u2(m)=u2(m)+du2dt_p(m)*dt;
       u3(m)=u3(m)+du3dt_p(m)*dt;      
    end
    
    %Predicted Values of flux terms using predicted values of u1,u2 and u3
    f1=u2;
    f2=((u2.^2)./u1)+((gamma-1)/gamma)*(u3-(gamma*u2.^2)./(2*u1));
    f3=(gamma.*u2.*u3./u1)-((gamma*(gamma-1)*u2.^3)./(2*u1.^2));
    
    
    %Corrector Step
    for p=2:n-1
        %Backward differencing for flux terms
        df1dx=(f1(p)-f1(p-1))/dx;
        df2dx=(f2(p)-f2(p-1))/dx;
        df3dx=(f3(p)-f3(p-1))/dx;
        dadx=(a(p)-a(p-1))/dx;
        j2=rho(p)*t(p)*dadx/gamma;
        
        %Governing Equations
        du1dt_c(p)=-df1dx;        %Continuity
        du2dt_c(p)=-df2dx+j2;     %Momentum
        du3dt_c(p)=-df3dx;        %Energy
    end
    
    %Computing average derivatives
    du1dt_av=0.5*(du1dt_p+du1dt_c);
    du2dt_av=0.5*(du2dt_p+du2dt_c);
    du3dt_av=0.5*(du3dt_p+du3dt_c);
    
    %Final Corrected Values of dependent variables u1,u2 and u3 at next time step
    u1=u1_old+(du1dt_av*dt);
    u2=u2_old+(du2dt_av*dt);
    u3=u3_old+(du3dt_av*dt);
    
    %Applying Boundary Conditions
    %Inflow Boundary
    u1(1)=rho(1)*a(1);      %Fixed Value as rho and area at inflow is fixed
    u2(1)=2*u2(2)-u2(3);
    v(1)=u2(1)/u1(1);
    u3(1)=u1(1)*((t(1)/(gamma-1))+((gamma/2)*v(1)^2));
    
    %Outflow boundary
    u1(n)=2*u1(n-1)-u1(n-2);
    u2(n)=2*u2(n-1)-u2(n-2);
    u3(n)=2*u3(n-1)-u3(n-2);
    
    %Finally obtaining the values of primitive variables by decoding u1,u2 and u3
    rho=u1./a;
    v=u2./u1;
    t=(gamma-1)*((u3./u1)-(gamma/2)*(v.^2));
    
    error_u1=max(abs(u1-u1_old));
    error_u2=max(abs(u2-u2_old));
    error_u3=max(abs(u3-u3_old));
    
    error=max([error_u1; error_u2; error_u3]);
    
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
    title(sprintf('Variations in properties at time step %d for %d gridpoints',l,n));
    legend('\rho/\rho_{o}','T/T_{o}','p/p_{o}','Mach');

    %Mach number in terms of non dimensionalized velocity and temperature
figure(2);
M=v./sqrt(t);
plot(x,M);
ylabel('\bf Mach');
xlabel('\bf Position');
title(sprintf('Time step %d',l));
grid on;
  

%Normalized Mass Flow Rate in terms of non dimensionalized variables
figure(3);
m=rho.*a.*v;   
plot(x,m);
xlabel('\bf Position');
ylabel('\bf Mass flow rate,\rho*A*V');
axis([0 3 0.5 0.6]);
title(sprintf('Mass Flow Rate at time step %d',l));
grid on;
    


%plotting profiles
figure(4);    
    subplot(5,1,1)
    plot(x,a);
    ylabel('\bf A/A^{\ast}');
    title(sprintf('Conservative Form and time step %d',l));
  
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
    p=rho.*t;
    plot(x,p);
    ylabel('\bf p/p_{o}');
    xlabel('\bf Position');
    title('Pressure');
    
    

