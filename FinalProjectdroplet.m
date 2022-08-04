
% level set method to capture the interface
clear all; clc; close all;
Lx = 1.0; Ly = 1.0;                                                  % domain size 
gx = 0.0; gy = -100.0; rho1 = 2; rho2 = 1; mu = 0.01;    % parameters
unorth = 0; usouth = 0; veast = 0; vwest = 0;                        % boundary conditions
rad = 0.15; xc = 0.5; yc = 0.5;                                      % initial drop size and location
D=1; % 
time = 0.0; plot_freq = 10;

nx = 256; ny = 256; dx = Lx/nx; dy = Ly/ny; dt = 0.0001;

nstep = 2000; maxit = 200; maxError = 0.001; omg = 1.5; Nf = 100;

u=zeros(nx+1,ny+2); ut = u   ; uplot = zeros(nx+1,ny+1);
v=zeros(nx+2,ny+1); vt = u   ; vplot = zeros(nx+1,ny+1);
 

p=zeros(nx+2,ny+2); tmp1 = p ; tmp2  = p; r = p; chi = p; 
C=zeros(nx+2,ny+2);  gamma=0; %zeros(nx+2,ny+2); 
Cn=C; % C at n time step
phi=zeros(nx+2, ny+2); % signed distance function
phin=zeros(nx+2, ny+2);% phi at n timestep
sphi=zeros(nx+2,ny+2); % signed phi in the reinitialization process
psi=zeros(nx+2,ny+2); % psi to satisfy the steady state solution 
psin=zeros(nx+2,ny+2); % fi at n timestep
eps= 1.5*dx; % used for smooth out chi
epsilon=0.0000001; % used for calculating sign of phi
Ddelta=zeros(nx+2,ny+2); % Dirac delta
%

% stargerred grid
xh = linspace(0,Lx,nx+1)         ; yh = linspace(0,Ly,ny+1);                % velocity points
x  = linspace(-dx/2,Lx+dx/2,nx+2); y  = linspace(-dy/2,Ly+dy/2,ny+2);       % pressure points

r = zeros(nx+2,ny+2) + rho2;                                                % initial density 
fgx=zeros(nx+2,ny+2); fgy=zeros(nx+2,ny+2);              % initial surface tension

% initialization 
% the initial shape is a circle
 % initialize front 
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
       chi(i,j)=0.0;
       Ddelta(i,j)=0;
   elseif(phi(i,j)>eps)
       chi(i,j)=1.;
       Ddelta(i,j)=0;
   else
       chi(i,j)=1/2.0 + phi(i,j)/(2*eps)+1/(2*pi)*sin(pi*phi(i,j)/eps);
       Ddelta(i,j)=1/(2*eps)*(1+cos(pi*phi(i,j)/eps));
   end;
end;end;
%r  = rho1*chi + rho2*(1-chi); % initial density
figure(1); contourf(x,y,phi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
figure(2); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
% for i = 2:nx+1; for j = 2:ny+1                   
%   if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2); r(i,j) = rho1; chi(i,j)=1.0; end; 
% end; end


                                
            
for is=1:nstep
    
    


%%%figure(4); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
%figure(3); contourf(x,y,chi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
              
  ro = r;
  r  = rho2*chi + rho1*(1-chi);  % obtain density from charact func

  figure(5); contourf(x,y,r'); axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Calculate the surface tension
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for i=2:nx+1; for j=2:ny+1
      fgx(i,j)=gamma * ((phi(i+1,j)-2.0 * phi(i,j)+phi(i-1,j))/(dx^2) + (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/(dy^2)) * Ddelta(i,j) * ((phi(i+1,j)-phi(i,j))/dx);  % using central
      fgy(i,j)=gamma * ((phi(i+1,j)-2.0 * phi(i,j)+phi(i-1,j))/(dx^2) + (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/(dy^2)) * Ddelta(i,j) * ((phi(i,j+1)-phi(i,j))/dy);
  end;
  end;
  fgx(1:nx+2,1)=fgx(1:nx+2,2);fgx(1:nx+2,ny+2)=fgx(1:nx+2,ny+1);
  fgx(1,1:ny+2)=fgx(2,1:ny+2);fgx(ny+2,1:ny+2)=fgx(ny+1,1:ny+2); 
  fgy(1:nx+2,1)=fgy(1:nx+2,2);fgy(1:nx+2,ny+2)=fgy(1:nx+2,ny+1);
  fgy(1,1:ny+2)=fgy(2,1:ny+2);fgy(ny+2,1:ny+2)=fgy(ny+1,1:ny+2); 
  
  %fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
 % fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  % boundary condition for surface tension  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving  N-S equations using projection method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
  v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC

  for i=2:nx; for j=2:ny+1     % temporary u-velocity (boundary values are not touched)
    ut(i,j) = (2.0/(r(i+1,j)+r(i,j)))*(0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+ dt* (...
            - (0.25/dx)*(ro(i+1,j)*(u(i+1,j)+u(i,j))^2-ro(i,j)*(u(i,j)+u(i-1,j))^2)...
            - (0.0625/dy)*( (ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1))*(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) ...
            - (ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))*(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))...
            + mu*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)...
            + 0.5*(ro(i+1,j)+ro(i,j))*gx + (fgx(i,j)+fgx(i-1,j))/2.0 ));
  end; end

  for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
    vt(i,j) = (2.0/(r(i,j+1)+r(i,j)))*(0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+ dt* (...     
            - (0.0625/dx)*( (ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) ...
            - (ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))*(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
            - (0.25/dy)*(ro(i,j+1)*(v(i,j+1)+v(i,j))^2-ro(i,j)*(v(i,j)+v(i,j-1))^2 )...
            + mu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)...
            + 0.5*(ro(i,j+1)+ro(i,j))*gy + (fgy(i,j)+fgy(i,j-1))/2.0 ) );    
  end; end     

  for i = 2:nx+1; for j = 2:ny+1
    tmp1(i,j) = (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy );
    tmp2(i,j) =1/( (1/dx)*(1/(dx*(r(i+1,j)+r(i,j)))+ 1/(dx*(r(i-1,j)+r(i,j))) )+ ...
                   (1/dy)*(1/(dy*(r(i,j+1)+r(i,j)))+ 1/(dy*(r(i,j-1)+r(i,j))) )   );
  end; end

  for it = 1:maxit	               % solve for pressure by SOR
    pold   = p;
    p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
    for i=2:nx+1; for j=2:ny+1
      p(i,j) = (1.0-omg)*p(i,j) + omg*tmp2(i,j)*(        ...
      (1/dx)*( p(i+1,j)/(dx*(r(i+1,j)+r(i,j)))+ p(i-1,j)/(dx*(r(i-1,j)+r(i,j))) )+    ...
      (1/dy)*( p(i,j+1)/(dy*(r(i,j+1)+r(i,j)))+ p(i,j-1)/(dy*(r(i,j-1)+r(i,j))) ) - tmp1(i,j));
    end; end
    if max(max(abs(pold-p))) < maxError; break;  end
  end
                                      
  for i=2:nx; for j=2:ny+1   % correct the u-velocity 
    u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
  end; end
      
  for i=2:nx+1; for j=2:ny   % correct the v-velocity 
    v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
  end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
% interface propagation
%%%%%%%%%%%%%%%%%%%%%%%%
phin=phi;

for i=2:nx+1;
    for j=2:ny+1;
        if(((u(i,j)+u(i-1,j))/2.0)>=0);
            dphi_x=(phin(i,j)-phin(i-1,j))/dx;  % upwind scheme
        else;
            dphi_x=(phin(i+1,j)-phin(i,j))/dx;
        end;
        if(((v(i,j)+v(i,j-1))/2.0)>=0);
            dphi_y=(phin(i,j)-phin(i,j-1))/dy;
        else
            dphi_y=(phin(i,j+1)-phin(i,j))/dy;
        end;
%        phi(i,j)= phin(i,j)+ dt * D*  ( (phin(i+1,j)-2*phin(i,j)+phin(i-1,j))/(dx^2) + (phin(i,j+1) - 2*phin(i,j) +phin(i,j-1))/(dy^2) ); %curvature flow
       phi(i,j)= phin(i,j)+ dt * ((u(i,j)+u(i-1,j))/2.0*dphi_x+(v(i,j)+v(i,j-1))/2.0*dphi_y);  %  drop falling flow 
        % phi(i,j) =phin(i,j) - u(i,j)   
    end;
end;
figure(4); contourf(x,y,phi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
% update phi to n+1 timestep 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interface reinitialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:nx+1;for j=2:ny+1;
        
sphi(i,j)=phi(i,j)/(sqrt(phi(i,j)^2 +epsilon^2));
end;end;
figure(3); contourf(x,y,sphi');  axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off
psi=phi;

for k=1:20; % iteration for 5 times
    psin=psi;
    for i=2:nx+1;
        for j=2:ny+1;
            dphix1=(psin(i,j)-psin(i-1,j))/dx;
            dphix2=(psin(i+1,j)-psin(i,j))/dx;
            dphiy1=(psin(i,j)-psin(i,j-1))/dy;
            dphiy2=(psin(i,j+1)-psin(i,j))/dy;
    
    % reinitialization near the interface
         if(psi(i,j)>0.);
             gphi=1-sqrt(max(max(dphix1,0.)^2, min(dphix2,0.)^2) +max(max(dphiy1, 0.)^2, min(dphiy2, 0.)^2));
         elseif(psi(i,j)<0.);
             gphi=1-sqrt(max(min(dphix1,0.)^2, max(dphix2,0.)^2)+max(min(dphiy1, 0.)^2, max(dphiy2, 0.)^2));
         else
             gphi=0.;
         end;
       % gphi=((psin(i+1,j)-psin(i-1,j))/(2*dx))^2 +((psin(i,j+1)-psin(i,j-1))/(2*dy))^2;
        psi(i,j)=psin(i,j)+dt*sphi(i,j)*gphi; 
        end;
    end;
    
   
end;

phi=psi;

%%%%%%%%%%get the characteristic function at new time step%%%%%%%
for i=2:nx+1;for j=2:ny+1
   if(phi(i,j) <-eps)  % inside the drop chi=0
       chi(i,j)=0.;   
       Ddelta(i,j)=0.;
   elseif(phi(i,j)>eps)  % outside the drop chi=1
       chi(i,j)=1.;
       Ddelta(i,j)=0.;
   else
       chi(i,j)=1/2.0 + phi(i,j)/(2*eps)+1/(2*pi)*sin(pi*phi(i,j)/eps);
       Ddelta(i,j)=1/(2*eps)*(1+cos(pi*phi(i,j)/eps));
    end;
end;end;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 

  time = time+dt                       
  if (mod(is,plot_freq)==0) | (is==1);                         % plot solution
  uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
  vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
  %figure(5); contourf(x,y,chi'); axis equal; axis([0 Lx 0 Ly]); hold on; drawnow; hold off %contour(x,y,r'); 
  figure(6); contourf(x,y,phi'); 
  axis equal; axis([0 Lx 0 Ly]); %hold on;drawnow; hold off
  quiver(xh,yh,uu',vv','r'); hold on; drawnow; hold off
  end

end
