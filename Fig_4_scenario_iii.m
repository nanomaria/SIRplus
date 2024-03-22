% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------


clear all
close all

global alphai C alphai_min delta gamma Hmaxi Hmax pi pr k h Hstar choice eta zip pd_min pd_max


choice =1; %we are choosing the sigmoid form of alpha(H) and p(H)

delta = .1;  % rate of decay of microbial diversity
h = .1;
gamma = 1/7;  % recovery rate
eta = 1/10; % recovery rate (severe cases)
alphai = 0.04;  % maximal probability of pathogenic transmission given a contact
alphai_min = 0.02; % minimal probability of pathogenic transmission given a contact
Hmaxi = 3;  %maximal diversity
Hmax = 6;  %minimal diversity
pi = 1;
pr = 1;
s0 = 1;
i0 = 1e-4;
is0 = 0;
k = 0.5;
Tfin = 400;

pd_max = 0.05;
pd_min = 0.03;

%C1 : C closed, C2: C open.

C1 = 5;
C2 = 10;


C = C1;

d0 = 0;

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd_close] = ode45(@onlyD_eq_updated, 0:.1:Tfin, d0, options);

Hstar_close = Yd_close(end,1);   
Hstar = Hstar_close;

zip = 1;
[T,Yopen] = ode45(@Eq_fig3_updated, 0:.1:Tfin, [s0,i0,is0], options);
Ic = (Yopen(:,2)+Yopen(:,3))*100;
Ic_sev = Yopen(:,3)*100;
zip=2;
[T,Yopen] = ode45(@Eq_fig3_updated, 0:.1:Tfin, [s0,i0,is0], options);
Ic_sev_plus = Yopen(:,3)*100;



C=C2;

[T,Yd_open] = ode45(@onlyD_eq_updated, 0:.1:Tfin, d0, options);
 
Hstar_open = Yd_open(end,1);   
Hstar = Hstar_open;

zip = 1;
[T,Yclose] = ode45(@Eq_fig3_updated, 0:.1:Tfin, [s0,i0,is0], options);
Io = (Yclose(:,2)+Yclose(:,3))*100;
Io_sev = Yclose(:,3)*100;
zip=2;
[T,Yclose] = ode45(@Eq_fig3_updated, 0:.1:Tfin, [s0,i0,is0], options);
Io_sev_plus = Yclose(:,3)*100;



figure(2)
a3=area(T,Ic_sev);
a3.FaceColor = [0.5 0.6 1];
alpha(a3,0.5);
hold on
f2= plot(T,Ic,'b-','Linewidth',1.5);
hold on
plot(T,Ic_sev,'Color',[0.5 0.6 1],'Linewidth',1)
a1=area(T,Ic_sev_plus);
a1.FaceColor = [0.5 0.6 1];
hold on
plot(T,Ic_sev_plus,'b-','Linewidth',1)
hold on
f1 = plot(T,Io,'k-','Linewidth',1.5);
hold on
plot(T,Io_sev_plus,'k-','Linewidth',1)
hold on
a3=area(T,Io_sev);
a3.FaceColor = [0.9 0.9 0.9];
a2=area(T,Io_sev_plus);
a2.FaceColor = [0.6 0.6 0.6];
alpha(a2,0.7);
legend([f1,f2],'Altered immune resilience', 'Unaltered immune resilience','Fontsize',14,'Location','northeast')
legend boxoff
xlabel('Time [arbitrary unit]')
ylabel({'Infected individuals';' [% of population]'})
set(gca,'fontsize',14)
hold on


figure(3)
a3=area(T,Ic_sev);
a3.FaceColor = [0.5 0.6 1];
alpha(a3,0.5);
hold on
%f2= plot(T,Ic,'b-','Linewidth',1.5);
hold on
plot(T,Ic_sev,'Color',[0.5 0.6 1],'Linewidth',1)
a1=area(T,Ic_sev_plus);
a1.FaceColor = [0.5 0.6 1];
hold on
plot(T,Ic_sev_plus,'b-','Linewidth',1)
hold on
%f1 = plot(T,Io,'k-','Linewidth',1.5);
hold on
plot(T,Io_sev_plus,'k-','Linewidth',1)
hold on
a3=area(T,Io_sev);
a3.FaceColor = [0.9 0.9 0.9];
a2=area(T,Io_sev_plus);
a2.FaceColor = [0.6 0.6 0.6];
alpha(a2,0.7);
axis([50 350 0 0.4])
xlabel('Time [arbitrary unit]')
ylabel({'Severe infections';' [% of population]'})
set(gca,'fontsize',14)
hold on



