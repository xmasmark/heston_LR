% Roe s s l e r system
 clear all;
 a =0.2; b=0.2; c=5.7;

 %right hand side of ODE
 roess=@(t,y) [-y(2)-y(3);y(1)+a*y(2);b+y(3)*(y(1)-c)] ;
 %ex2=@(t,y1) [y1(2)^2-y1(2)^4;y1(1)-y1(2)^2];
 ex2=@(t,y1) [y1(2)^2-y1(2)^4;y1(1)-y1(2)^2];

 tspan =[0 200] ;
 tspan1 =[0 1] ;
 y0 = [0;0;0] ;
 y10 = [-10;10] ;

 % solve
 [t,y]=ode45(roess,tspan,y0) ;
 [t1,y1]=ode45(ex2,tspan1,y10) ;
 last=length(t);

 figure(1);
 clf;
 % plot3(y(:,1),y(:, 2),y(:,3));
 % xlabel('x') ;
 % ylabel('y') ;
 % zlabel('z') ;
 % figure(2) ;
 % clf;
 % plot(t,y(:,1));
 % xlabel('t');
 % ylabel('x');

% plot(t1,y1(:,1));
% xlabel('t1');
% ylabel('x');
% 
 plot(y1(:,1),y1(:, 2));
 xlabel('x') ;
 ylabel('y');

% 
% ODExvf = @(t,y,m,D,K,T,G) [y(2);-D/m*y(2)-K/m*y(1)-y(3)/m;-y(3)/T + G/T*y(1)];
% m=3;
% D=2;
% K=1;
% T=1;
% G=5;
% [t,F]=ode45(@(t,y) ODExvf(t,y,m,D,K,T,G),[0,1],[0.5,0,0]);
% plot(t,F);