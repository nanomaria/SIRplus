% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

global C delta Hmaxi Hmax k h Hstar

delta = .1;  % rate of decay of H
h = .1; % rate of acquisition of H
p_max = 0.05;  % maximal probability of severe disease from infection
p_min = 0.03; % minimal probability of severe disease from infection
Hmaxi = 5; %
Hmax = 10; %
Tfin = 300;

dC = .1;  % 0.1 for smooth curves


Cv = 0:dC:30;
kv = [0.25,0.5,0.75,1] ;
Nsteps = size(Cv,2);
Msteps = size(kv,2);

H0 = 0;


Hv = zeros(Msteps,Nsteps);

for j = 1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq_updated, 0:.1:Tfin, H0, options);

Hv(j,i) = Yd(end,1);

end
end


Nsteps = size(Cv,2);
Msteps = size(kv,2);
pv = zeros(Msteps,Nsteps);

Hstar = zeros(1,Nsteps);

H0 = 0;


for j=1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);


options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq_updated, 0:.1:Tfin, H0, options);
 
Hstar(j,i) = Yd(end,1);   


pv(j,i) = (p_min + (p_max-p_min)/2*(1-tanh((Hstar(j,i)-Hmaxi))));



end




end










figure(3);
%subplot(2,1,1)
fig = figure(3);
right_color =  [0 0.462 0.729];
left_color = [0.592 0.054 0.325];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




yyaxis right
plot(Cv,Hv(4,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1.8)
hold on
line(Cv,Hv(3,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1.4)
hold on
line(Cv,Hv(2,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1)
hold on
line(Cv,Hv(1,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',0.6)
hold on
axis([0 Cv(end) 0 11])
xlabel('Contact rate (C)')
ylabel({'Interaction-dependent health benefits',' at equilibrium (H*)'})
set(gca,'fontsize',14)



yyaxis left
plot(Cv,pv(4,:),'Color',[0.592 0.054 0.325],'Linestyle','-.','Linewidth',1.8)
hold on
line(Cv,pv(3,:),'Color',[0.592 0.054 0.325],'Linestyle','-.','Linewidth',1.4)
hold on
line(Cv,pv(2,:),'Color',[0.592 0.054 0.325],'Linestyle','-.','Linewidth',1)
hold on
line(Cv,pv(1,:),'Color',[0.592 0.054 0.325],'Linestyle','-.','Linewidth',0.6)
hold on
axis([0 Cv(end) -.045 0.25])
hold on
ylabel({'Probability of severe disease (p(H))'})
set(gca,'fontsize',14)


tcy = 0.5; %transparency
p1 = -0.0375;
p2 = -.025;
p3 = -0.0125;
p4 = -0.0;
x = 0.005;
thc = 4; %thickness

d1=line([ Hmaxi*delta/(kv(1)*h),Cv(end)],[p1,p1],   'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',thc);
d1.Color(4) = tcy;
hold on
d2=line([ Hmaxi*delta/(kv(2)*h),Cv(end)],[p2,p2], 'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',thc);
d2.Color(4) = tcy;
hold on
d3=line([ Hmaxi*delta/(kv(3)*h),Cv(end)],[p3,p3],   'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',thc);
d3.Color(4) = tcy;
hold on
d4=line([ Hmaxi*delta/(kv(4)*h),Cv(end)],[p4,p4], 'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',thc);
d4.Color(4) = tcy;
hold on

b1=line([ 0,Hmax*delta/(kv(1)*h)],[p1+x,p1+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b1.Color(4) = tcy;
hold on
b2=line([ 0,Hmax*delta/(kv(2)*h)],[p2+x,p2+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b2.Color(4) = tcy;
hold on
b3=line([ 0,Hmax*delta/(kv(3)*h)],[p3+x,p3+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b3.Color(4) = tcy;
hold on
b4=line([ 0,Hmax*delta/(kv(4)*h)],[p4+x,p4+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b4.Color(4) = tcy;
hold on


