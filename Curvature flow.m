
% level set method to capture the interface

clear all; clc; close all;
Lx = 1.0; Ly = 1.0;                                                  % domain size 
gx = 0.0; gy = -100.0; rho1 = 2; rho2 = 1; mu = 0.01;    % parameters
unorth = 0; usouth = 0; veast = 0; vwest = 0;                        % boundary conditions
rad = 0.15; xc = 0.5; yc = 0.5;                                      % initial drop size and location
D=1; % diffusion coefficient
time = 0.0; plot_freq = 30; pi=3.415926;

nx = 256; ny = 256; dx = Lx/nx; dy = Ly/ny; dt = 0.00001;

nstep = 1200; maxit = 200; maxError = 0.001; omg = 1.5; 

u=zeros(nx+1,ny+2); ut = u   ; uplot = zeros(nx+1,ny+1);
v=zeros(nx+2,ny+1); vt = v   ; vplot = zeros(nx+1,ny+1);  % vt=u???


p=zeros(nx+2,ny+2); tmp1 = p ; tmp2  = p; r = p; chi = p; 
C=zeros(nx+2,ny+2);  gamma=1; %zeros(nx+2,ny+2); 
Cn=C; % C at n time step
phi=zeros(nx+2, ny+2); % signed distance function
phin=zeros(nx+2, ny+2);% phi at n timestep
sphi=zeros(nx+2,ny+2); % signed phi in the reinitialization process
psi=zeros(nx+2,ny+2); % psi to satisfy the steady state solution 
psin=zeros(nx+2,ny+2); % psi at n timestep
eps= 1.5*dx; % used for smooth out chi
epsilon=0.00001; % used for calculating sign of phi
Ddelta=zeros(nx+2,ny+2); % Dirac delta
%
R1=zeros(1,nstep+1); R1(1)=rad; % Radius of the circle
%xf=zeros(1,Nf+2); yf=zeros(1,Nf+2); 
%un=zeros(1,Nf+2); vn=zeros(1,Nf+2);

% stargerred grid
xh = linspace(0,Lx,nx+1)         ; yh = linspace(0,Ly,ny+1);                % velocity points
x  = linspace(-dx/2,Lx+dx/2,nx+2); y  = linspace(-dy/2,Ly+dy/2,ny+2);       % pressure points

r = zeros(nx+2,ny+2) + rho2;                                                % initial density 
fgx=zeros(nx+2,ny+2); fgy=zeros(nx+2,ny+2);              % initial surface tension

% initialization


for i=2:nx+1;
    for j=2:ny+1;
        u(i,j)=1;
    end;
end;
    

% the initial shape is a circle

for i=1:nx+2;for j=1:ny+2
   if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2);    
        phi(i,j)= -( rad-sqrt((x(i)-xc)^2 + (y(j)-yc)^2));
        %chi(i,j)= 1/2.0 + phi(i,j)/(2*eps)+1/(2*pi)*sin(pi*phi(i,j)/eps);%1.0;
        r(i,j) = rho1;
   elseif((x(i)-xc)^2+(y(j)-yc)^2 > rad^2);
        phi(i,j)=  sqrt((x(i)-xc)^2 + (y(j)-yc)^2)-rad;  
        %chi(i,j)=0.0;
   else
       phi(i,j)=0.0;
       %chi(i,j)=0.0;
   end;
   if(phi(i,j) <-eps)
       chi(i,j)=0;
       Ddelta(i,j)=0;
   elseif(phi(i,j)>eps)
       chi(i,j)=1.;
       Ddelta(i,j)=0;
   else
       chi(i,j)=1/2.0 + phi(i,j)/(2*eps)+1/(2*pi)*sin(pi*phi(i,j)/eps);
       Ddelta(i,j)=1/(2*eps)*(1+cos(pi*phi(i,j)/eps));
   end;
end;end;

figure(1); contourf(x,y,phi');  axis equal; axis([0 Lx 0 Ly]); hold on; xlabel('x', 'Fontsize', 20);
ylabel('y', 'Fontsize', 20);
title([sprintf('time t=%0.3f', time)], 'Fontsize',20);
colorbar;%caxis([-1 1])
drawnow; hold off
figure(2); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; xlabel('x', 'Fontsize', 20);
ylabel('y', 'Fontsize', 20);
title([sprintf('time t=%0.3f', time)], 'Fontsize',20);
colorbar;%caxis([-1 1])
drawnow; hold off
% for i = 2:nx+1; for j = 2:ny+1                   
%   if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2); r(i,j) = rho1; chi(i,j)=1.0; end; 
% end; end


                                
            
for is=1:nstep
    
    
%%%%%%%%%%%%%%%%%%%%%%%
% interface propagation
%%%%%%%%%%%%%%%%%%%%%%%%
phin=phi;

%%%% curvature flow  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%validation 
for it = 1:maxit
    
    phi_old=phi;
    
    for i=2:nx+1;
        for j=2: ny+1;
            
            phi(i,j)= (omg)*1.0/(1.0 + 2.0*dt/dx^2*D+2*dt/dy^2*D)* ...
                      (((phi(i+1,j)+phi(i-1,j))*dt/dx^2*D +(phi(i,j+1)+phi(i,j-1))*dt/dy^2*D) +phin(i,j) )+ ...
                      (1-omg)*phi(i,j);
        end;
    end;
    phi(1,:) = phi(2,:); phi(nx+2,:) =phi(nx+1, :); phi(:,1)=phi(:,2); phi(:,ny+2) =phi(:,ny+1) ; % no flux B.C, set ghost values
    if max(max(abs(phi_old-phi))) < maxError; break;  end
end


% figure(4); contourf(x,y,phi');  axis equal; axis([0 Lx 0 Ly]); hold on; 
% xlabel('x', 'Fontsize', 20);
% ylabel('y', 'Fontsize', 20);
% title([sprintf('time t=%0.3f', time)], 'Fontsize',20);
% colorbar;%caxis([-1 1])
% drawnow; hold off
% update phi to n+1 timestep 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interface reinitialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nx+1;for j=2:ny+1;
        
sphi(i,j)= phi(i,j)/(sqrt(phi(i,j)^2 +epsilon^2));
end;end;

figure(3); contourf(x,y,sphi');  axis equal; axis([0 Lx 0 Ly]); hold on; 
xlabel('x', 'Fontsize', 20);
ylabel('y', 'Fontsize', 20);
title([sprintf('time t=%0.3f', time)], 'Fontsize',20);
colorbar;%caxis([-1 1])
drawnow; hold off

psi=phi; % initial value of the reinitialization process


% if 

for k=1:5; % iteration for 5 time steps
    psin=psi;
    for i=2:nx+1;
        for j=2:ny+1;
            dphix1=(psi(i,j)-psi(i-1,j))/dx;
            dphix2=(psi(i+1,j)-psi(i,j))/dx;
            dphiy1=(psi(i,j)-psi(i,j-1))/dy;
            dphiy2=(psi(i,j+1)-psi(i,j))/dy;
    
            % upwind scheme
           if(psi(i,j)>0.);
             gphi=1.0-sqrt(max(max(dphix1,0.)^2, min(dphix2,0.)^2) +max(max(dphiy1, 0.)^2, min(dphiy2, 0.)^2));
            elseif(psi(i,j) <  0.);
             gphi=1.0-sqrt(max(min(dphix1,0.)^2, max(dphix2,0.)^2)+max(min(dphiy1, 0.)^2, max(dphiy2, 0.)^2));
           else
             gphi=0.;
           end;
        psi(i,j)=psin(i,j)+dt*sphi(i,j)*gphi; 
        end;
    end;
end;
% 
% phi=psi;
    
    
%calculate the characteristic function    
for i=2:nx+1;for j=2:ny+1
   if(phi(i,j) <-eps)
       chi(i,j)=0.;
       Ddelta(i,j)=0.;
   elseif(phi(i,j)>eps)
       chi(i,j)=1.;
       Ddelta(i,j)=0.;
   else
       chi(i,j)=1/2.0 + phi(i,j)/(2*eps)+1/(2*pi)*sin(pi*phi(i,j)/eps);
       Ddelta(i,j)=1/(2.0*eps)*(1+cos(pi*phi(i,j)/eps));
    end;
end;end;

%%%figure(4); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
%figure(3); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
              
  ro = r;
  r  = rho1*chi + rho2*(1-chi);  % obtain density from charact func


  % interpolate the values of the concentration to the front
%   for l=1:Nf+1
%       ip = floor(xf(l)/dx)+1; jp = floor((yf(l))/dy)+1;  % find the closest point
%       ax = xf(l)/dx-ip+1; ay = (yf(l))/dy-jp+1;	   
%       cf(l) = (1.0-ax)*(1.0-ay)*C(ip,jp) + ax*(1.0-ay)*C(ip+1,jp) + (1.0-ax)*ay*C(ip,jp+1) + ax*ay*C(ip+1,jp+1);
%   end
%   cf(Nf+2)=cf(2);
%   gammaCf=10.0 -3.0.*cf; % surface tension at the front
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Calculate the surface tension
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for i=2:nx+1; for j=2:ny+1
      fgx(i,j)=gamma * ((phi(i+1,j)-2.0 * phi(i,j)+phi(i-1,j))/(dx^2) + (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/(dy^2)) * Ddelta(i,j) * ((phi(i+1,j)-phi(i,j))/dx);
      fgy(i,j)=gamma * ((phi(i+1,j)-2.0 * phi(i,j)+phi(i-1,j))/(dx^2) + (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/(dy^2)) * Ddelta(i,j) * ((phi(i,j+1)-phi(i,j))/dy);
  end;
  end;
  
   %fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
   %fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  % boundary condition for surface tension  
  
%   for l=2:Nf+1  % distribute of Fr (surface tension in a volume)to the fixed grid
%     fglx = (gammaCf(l)*tx(l)-gammaCf(l-1)*tx(l-1)); fgly = (gammaCf(l)*ty(l)-gammaCf(l-1)*ty(l-1)); %fglx = gamma*(tx(l)-tx(l-1)); fgly = gamma*(ty(l)-ty(l-1));  % calculate Fr surface tension per volume
%     ip   = floor(xf(l)/dx)+1;     jp   = floor((yf(l)+0.5*dy)/dy)+1;  
%     % closest point as integration point
%     
%     ax   = xf(l)/dx-ip+1;         ay   = (yf(l)+0.5*dy)/dy-jp+1;
%     fgx(ip,jp)     = fgx(ip,jp)     + (1.0-ax)*(1.0-ay)*fglx/dx/dy;
%     fgx(ip+1,jp)   = fgx(ip+1,jp)   + ax*(1.0-ay)*fglx/dx/dy;
%     fgx(ip,jp+1)   = fgx(ip,jp+1)   + (1.0-ax)*ay*fglx/dx/dy;
%     fgx(ip+1,jp+1) = fgx(ip+1,jp+1) + ax*ay*fglx/dx/dy;
% 
%     ip = floor((xf(l)+0.5*dx)/dx)+1; jp = floor(yf(l)/dy)+1;
%     ax = (xf(l)+0.5*dx)/dx-ip+1;     ay = yf(l)/dy-jp+1;	  
%     fgy(ip,jp)     = fgy(ip,jp)     + (1.0-ax)*(1.0-ay)*fgly/dx/dy;
%     fgy(ip+1,jp)   = fgy(ip+1,jp)   + ax*(1.0-ay)*fgly/dx/dy;
%     fgy(ip,jp+1)   = fgy(ip,jp+1)   + (1.0-ax)*ay*fgly/dx/dy;
%     fgy(ip+1,jp+1) = fgy(ip+1,jp+1) + ax*ay*fgly/dx/dy;  
%   end
% 
%   fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
%   fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  % boundary condition for surface tension
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving  N-S equations using projection method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
%   v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC
% 
%   for i=2:nx; for j=2:ny+1     % temporary u-velocity (boundary values are not touched)
%     ut(i,j) = (2.0/(r(i+1,j)+r(i,j)))*(0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+ dt* (...
%             - (0.25/dx)*(ro(i+1,j)*(u(i+1,j)+u(i,j))^2-ro(i,j)*(u(i,j)+u(i-1,j))^2)...
%             - (0.0625/dy)*( (ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1))*(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) ...
%             - (ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))*(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))...
%             + mu*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)...
%             + 0.5*(ro(i+1,j)+ro(i,j))*gx + fgx(i,j) ) );
%   end; end
% 
%   for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
%     vt(i,j) = (2.0/(r(i,j+1)+r(i,j)))*(0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+ dt* (...     
%             - (0.0625/dx)*( (ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) ...
%             - (ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))*(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
%             - (0.25/dy)*(ro(i,j+1)*(v(i,j+1)+v(i,j))^2-ro(i,j)*(v(i,j)+v(i,j-1))^2 )...
%             + mu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)...
%             + 0.5*(ro(i,j+1)+ro(i,j))*gy + fgy(i,j) ) );    
%   end; end     
% 
%   for i = 2:nx+1; for j = 2:ny+1
%     tmp1(i,j) = (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy );
%     tmp2(i,j) =1/( (1/dx)*(1/(dx*(r(i+1,j)+r(i,j)))+ 1/(dx*(r(i-1,j)+r(i,j))) )+ ...
%                    (1/dy)*(1/(dy*(r(i,j+1)+r(i,j)))+ 1/(dy*(r(i,j-1)+r(i,j))) )   );
%   end; end
% 
%   for it = 1:maxit	               % solve for pressure by SOR
%     pold   = p;
%     p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
%     for i=2:nx+1; for j=2:ny+1
%       p(i,j) = (1.0-omg)*p(i,j) + omg*tmp2(i,j)*(        ...
%       (1/dx)*( p(i+1,j)/(dx*(r(i+1,j)+r(i,j)))+ p(i-1,j)/(dx*(r(i-1,j)+r(i,j))) )+    ...
%       (1/dy)*( p(i,j+1)/(dy*(r(i,j+1)+r(i,j)))+ p(i,j-1)/(dy*(r(i,j-1)+r(i,j))) ) - tmp1(i,j));
%     end; end
%     if max(max(abs(pold-p))) < maxError; break;  end
%   end
%                                       
%   for i=2:nx; for j=2:ny+1   % correct the u-velocity 
%     u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
%   end; end
%       
%   for i=2:nx+1; for j=2:ny   % correct the v-velocity 
%     v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
%   end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=2:nx+1;
%     for j=2:ny+1;
%         phi(i,j)=phin(i,j)
%     end;
% end;


% calculate the radius of the circle;

   for i=(round((nx+3)/2)):nx+1;
       %for j=ny/2+1:ny+2;
           if phi(i,round((ny+3)/2)) < 0 && phi(i+1, round((ny+3)/2))>0;
               x1 =x(i);
               x2=x(i+1);
               c1= phi(i, round((ny+3)/2));
               c2= phi(i+1, round((ny+3)/2));
               xr= x1-c1/(c2-c1)*(x2-x1);
               R1(is+1)=(xr-0.5);
               %R(is+1)=sqrt((x(i)-xc)^2 + (y(i)-yc)^2);
               break;
           end;
       %end;
   end;

  time = time+dt                       
  if (mod(is,plot_freq)==0) | (is==1);                         % plot solution
%  uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
%  vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
  %figure(5); contourf(x,y,chi'); axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off %contour(x,y,r'); 
  figure(6); contourf(x,y,phi'); 
  axis equal; axis([0 Lx 0 Ly]); hold on;
  xlabel('x', 'Fontsize', 20);
  ylabel('y', 'Fontsize', 20);
title([sprintf('time t=%0.3f', time)], 'Fontsize',20);
colorbar;
%caxis([-1 1]);
drawnow; hold off
%  quiver(xh,yh,uu',vv','r'); hold on; drawnow; hold off
  end

end
t1=linspace(0,0.012,nstep+1);
r_exact=sqrt(rad^2-2*t1);
figure(7);
plot(t1,r_exact); hold on;
plot(t1,R1,'--r'); 
legend('Exact solution', 'level set')
xlabel('t', 'Fontsize', 14);
ylabel('r', 'Fontsize', 14);
drawnow; hold off;
