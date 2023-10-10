% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

global alphai C alphai_min delta gamma epsi epsm Dmaxi Dmaxm pi pr k Dstar sR iR rR choice


choice =1; %we are choosing the tanh curve for alpha(D)
delta = 1;  % rate of decay of microbial diversity
gamma = 1/7;  % recovery rate
alphai = 0.1;  % maximal probability of pathogenic transmission given a contact
alphai_min = 0.02; % minimal probability of pathogenic transmission given a contact
epsi = 1; % smoothering factor for function 1
epsm = 1; % smoothering factor for function 2
Dmaxi = 3;  %maximal diversity
Dmaxm = 6;  %minimal diversity
pi = 1;
pr = 1;
s0 = 1;
i0 = 1e-6;
r0 = 0;
c0 = 0;
k = .75;
Tfin = 300;

p_max = 0.1;
p_min = 0.02;

%C1 : C closed, C2: C open.

C1 = 3;
C2 = 8;


C = C1;

d0 = 0;
sR = 1;
iR = 0;
rR = 0;
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);


Dstar_close = Yd(end,1);   
Dstar = Dstar_close;

C=C2;
i0 = 1e-6;
[T,Yclose] = ode45(@SIR_square_cumul_eq, 0:.1:Tfin, [s0,i0,r0,c0], options);



C=C1;

d0 = 0;  
iR = 0.1;
rR = 0.6;
sR = 1-iR-rR;


options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);


Dstar_close_i = Yd(end,1);   
Dstar = Dstar_close_i;

C=C2;
sR = 1;
iR = 0;
rR = 0;

i0 = 1e-6;
[T,Yclosei] = ode45(@SIR_square_cumul_eq, 0:.1:Tfin, [s0,i0,r0,c0], options);


sR = 1;
iR = 0;
rR = 0;
C = 10;
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);
 
Dstar_open = Yd(end,1);   
Dstar = Dstar_open;

i0 = 1e-5;
[T,Yopen] = ode45(@SIR_square_cumul_eq, 0:.1:Tfin, [s0,i0,r0,c0], options);


figure(2)
a3=area(T,Yclose(:,2)*p_max*100);
a3.FaceColor = [0.5 0.6 1];
alpha(a3,0.5);
hold on
f2= plot(T,Yclose(:,2)*100,'b-','Linewidth',1.5);
hold on
plot(T,Yclose(:,2)*p_max*100,'Color',[0.5 0.6 1],'Linewidth',1)
a1=area(T,Yclose(:,2)*p_min*100);
a1.FaceColor = [0.5 0.6 1];
hold on
plot(T,Yclose(:,2)*p_min*100,'b-','Linewidth',1)
hold on
f1 = plot(T,Yopen(:,2)*100,'k-','Linewidth',1.5);
hold on
plot(T,Yopen(:,2)*p_min*100,'k-','Linewidth',1)
hold on
a3=area(T,Yopen(:,2)*p_max*100);
a3.FaceColor = [0.9 0.9 0.9];
a2=area(T,Yopen(:,2)*p_min*100);
a2.FaceColor = [0.6 0.6 0.6];
alpha(a2,0.7);
axis([0 300 0 max(Yclose(:,2))*1.05*100])
legend([f1,f2],'No physical distancing', 'Prolonged physical distancing','Fontsize',12,'Location','north')
legend boxoff
xlabel('Time [arbitrary unit]')
ylabel('Infected individuals [% of population]')
set(gca,'fontsize',14)
hold on




figure(1)
a3=area(T,Yclose(:,2)*p_max*100);
a3.FaceColor = [0.5 0.6 1];
alpha(a3,0.5);
hold on
plot(T,Yclose(:,2)*100,'b-','Linewidth',1.5);
hold on
plot(T,Yclose(:,2)*p_max*100,'Color',[0.5 0.6 1],'Linewidth',1)
a1=area(T,Yclose(:,2)*p_min*100);
a1.FaceColor = [0.5 0.6 1];
hold on
plot(T,Yclose(:,2)*p_min*100,'b-','Linewidth',1)
hold on
plot(T,Yopen(:,2)*100,'k-','Linewidth',1.5);
hold on
plot(T,Yopen(:,2)*p_min*100,'k-','Linewidth',1)
hold on
plot(T,Yopen(:,2)*p_max*100,'k-','Linewidth',1)
hold on
a3=area(T,Yopen(:,2)*p_max*100);
a3.FaceColor = [0.9 0.9 0.9];
a2=area(T,Yopen(:,2)*p_min*100);
a2.FaceColor = [0.6 0.6 0.6];
alpha(a2,0.7);
axis([0 300 0 max(Yclose(:,2))*p_max*100*1.1])
xlabel('Time [arbitrary unit]')
ylabel('% hospitalization')
set(gca,'fontsize',14)
hold on



