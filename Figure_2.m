% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------


clear all
close all

global alphai C alphai_min delta gamma epsi epsm Dmaxi Dmaxm pi pr k Dstar sR iR rR choice

delta = 1;  % rate of decay of microbial diversity
gamma = 1/7;  % recovery rate
alpha_standard = 0.1;
alphai = alpha_standard;  % maximal probability of pathogenic transmission given a contact
alphai_min = 0.02; %0.02% minimal probability of pathogenic transmission given a contact
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
k = .5;
Tfin = 300;
dc = 2;  % choose 0.1 for smooth curves



Cv = 0:dc:30;
kv=[0.25,0.5,0.75,1];
Nsteps = size(Cv,2);
Msteps = size(kv,2);
Imax = zeros(Msteps,Nsteps);
Imax_standard = zeros(Msteps,Nsteps);
Cmax = zeros(Msteps,Nsteps);

d0 = 0;
sR= s0;
iR= i0;
rR= r0;




choice = 1;
% if
% 1 : tanh
% 2 : linear function



for j = 1:Msteps

    k=kv(j);

for i=1:Nsteps
    C = Cv(i);


options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);
 
Dstar = Yd(end,1);   


[T,Y] = ode45(@SIR_square_cumul_eq, 0:.1:Tfin, [s0,i0,r0,c0], options);

Imax(j,i) = max(max(Y(:,2)))*100;


[T,Ys] = ode45(@SIR_standard_eq, 0:.1:Tfin, [s0,i0,r0], options);
Imax_standard(j,i) = max(max(Ys(:,2)))*100;
alphai =  alphai_min;
[T,Ys_m] = ode45(@SIR_standard_eq, 0:.1:Tfin, [s0,i0,r0], options);
Imax_standard_m(j,i) = max(max(Ys_m(:,2)))*100;
alphai = alpha_standard;



end

end



Dv = zeros(Msteps,Nsteps);

for j = 1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);

Dv(j,i) = Yd(end,1);

end
end







figure(2)
subplot(1,3,2)
grey = .2;
plot(Cv,Imax(4,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1.8);
hold on
line(Cv,Imax(3,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1.4)
hold on
l1=line(Cv,Imax(2,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1);
hold on
line(Cv,Imax(1,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',0.6)
hold on
l2=plot(Cv,Imax_standard(1,:),'r--','LineWidth',1.4);
hold on
l3=plot(Cv,Imax_standard_m(1,:),'b--','Linewidth',1.4);
set(gca,'fontsize',14)
axis([0 Cv(end) 0 100])
xlabel('Contact rate (C)')
ylabel({'Peak infection [% of population]'})
legend([l2,l3,l1],{'\alpha_{max}','\alpha_{min}', '\alpha(H)'},'Fontsize',12,'NumColumns',3,'location','north')





Nsteps = size(Cv,2);
Msteps = size(kv,2);
betav1 = zeros(Msteps,Nsteps);
betav1m = zeros(Msteps,Nsteps);
betav2 = zeros(Msteps,Nsteps);
betav = zeros(Msteps,Nsteps);
alpha_lin = zeros(Msteps,Nsteps);
Dstar = zeros(1,Nsteps);

d0 = 0;
sR= s0;
iR= i0;
rR= r0;


for j=1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);


options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);
 
Dstar(j,i) = Yd(end,1);   


betav1(j,i) =   C*alphai;
betav1m(j,i) = C*alphai_min;

if choice == 1
betav2(j,i) = (alphai_min + (alphai-alphai_min)/2*(1-tanh((Dstar(j,i)-Dmaxi)*epsi)));
betav(j,i) = C*(alphai_min + (alphai-alphai_min)/2*(1-tanh((Dstar(j,i)-Dmaxi)*epsi)));
end

if choice == 2
%linear function T(D)
if Dstar(j,i) <= Dmaxi*2
alpha_lin(j,i) = -(alphai-alphai_min)/(Dmaxi*2)*Dstar(j,i)+alphai;
elseif Dstar(j,i)>Dmaxi*2
alpha_lin(j,i) = alphai_min;
end
betav(j,i) = C*alpha_lin(j,i);
end



end

end






tcy = 0.5; %transparency
p1 = 0.3/gamma*(alpha_standard/0.2);
p2 = .42/gamma*(alpha_standard/0.2);
p3 = 0.6/gamma*(alpha_standard/0.2);
p4 = 1.2/gamma*(alpha_standard/0.2);

x = 0;

ind4 = find(diff(betav(4,:)) < 0);
ind3 = find(diff(betav(3,:)) < 0);
ind2 = find(diff(betav(2,:)) < 0);
ind1 = find(diff(betav(1,:)) < 0);

figure(2)
subplot(1,3,1)
plot([0,Cv(end)],[1 1],'k:','Linewidth',0.5)
hold on
if sum(ind4) > 0
c1=plot([ind4(1)*dc,ind4(end)*dc],[p1-x,p1-x],'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',5);
c1.Color(4) = tcy;
end
hold on
if sum(ind3) > 0
c2=line([ind3(1)*dc,ind3(end)*dc],[p2-x,p2-x],'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',5);
c2.Color(4) = tcy;
end
hold on
if sum(ind2) > 0
c3=line([ind2(1)*dc,ind2(end)*dc],[p3-x,p3-x],'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',5);
c3.Color(4) = tcy;
end
hold on
if sum(ind1) > 0
c4=line([ ind1(1)*dc,ind1(end)*dc],[p4-x,p4-x],'Color',[0.592 0.054 0.325],'Linestyle','-','LineWidth',5);
c4.Color(4) = tcy;
end
h1=plot(Cv,betav1(1,:)/gamma,'Color',[1 0 0],'Linestyle','--','LineWidth',1.4);
hold on
line(Cv,betav(4,:)/gamma,'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',1.8);
hold on 
line(Cv,betav(3,:)/gamma,'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',1.4)
hold on 
h2=line(Cv,betav(2,:)/gamma,'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',1);
line(Cv,betav(1,:)/gamma,'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',0.6);
h3=line(Cv,betav1m(1,:)/gamma,'Color',[0 0 1],'Linestyle','--','LineWidth',1.4);
set(gca,'fontsize',14)
axis([0 Cv(end) 0 3.5/gamma*(alpha_standard/0.2)])
legend([h1,h3,h2],{'\alpha_{max}','\alpha_{min}', '\alpha(H)'},'Fontsize',12,'NumColumns',3,'location','north')
xlabel('Contact rate (C)')
ylabel('Basic reproduction number (R_0)')
set(gca,'fontsize',14)




%%% Disease tolerance



delta = 1;  % rate of decay of microbial diversity
gamma = 1/7;  % recovery rate
alpha_standard = 0.1;
alphai = alpha_standard;  % maximal probability of pathogenic transmission given a contact
epsi = 1; % smoothering factor for function 1
epsm = 1; % smoothering factor for function 2
Dmaxi = 3;  %maximal diversity
Dmaxm = 6;  %minimal diversity
pd_min = 0.02;
pd_max = 0.1;
pi = 1;
pr = 1;
s0 = 1;
i0 = 1e-6;
r0 = 0;
c0 = 0;
k = .5;
Tfin = 300;



Cv = 0:dc:30;
kv= [0.25,0.5,0.75,1] ;
Nsteps = size(Cv,2);
Msteps = size(kv,2);
Imax_sev = zeros(Msteps,Nsteps);
Imax_standard = zeros(Msteps,Nsteps);
Imax_standard_sev_low = zeros(Msteps,Nsteps);
Imax_standard_sev_high = zeros(Msteps,Nsteps);
Cmax = zeros(Msteps,Nsteps);

d0 = 0;
sR= s0;
iR= i0;
rR= r0;



for j = 1:Msteps

    k=kv(j);

for i=1:Nsteps
    C = Cv(i);


options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);
 
Dstar = Yd(end,1);   



if choice ==1
pd = (pd_min + (pd_max-pd_min)/2*(1-tanh((Dstar-Dmaxi)*epsi)));
end


if choice == 2
%linear function alpha(D)
if Dstar<=Dmaxi*2
pd = (-(pd_max-pd_min)/(Dmaxi*2)*Dstar+pd_max);
elseif Dstar>Dmaxi*2
pd = pd_min;
end
end



[T,Y] = ode45(@SIR_square_sev_eq, 0:.1:Tfin, [s0,i0,r0], options);

Imax_sev(j,i) = pd*max(max(Y(:,2)))*100;


[T,Ys] = ode45(@SIR_standard_eq, 0:.1:Tfin, [s0,i0,r0], options);
Imax_standard(j,i) = max(max(Ys(:,2)))*100;
Imax_standard_sev_low(j,i) = pd_min*Imax_standard(j,i);
Imax_standard_sev_high(j,i) = pd_max*Imax_standard(j,i);



end

end



Dv = zeros(Msteps,Nsteps);

for j = 1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Yd] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);

Dv(j,i) = Yd(end,1);

end
end






figure(2)
subplot(1,3,3)
grey = .2;
plot(Cv,Imax_sev(4,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1.8);
hold on
line(Cv,Imax_sev(3,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1.4);
hold on
l1=line(Cv,Imax_sev(2,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',1);
hold on
line(Cv,Imax_sev(1,:),'Color',[grey grey grey],'Linestyle','-','Linewidth',0.6)
hold on
l2=plot(Cv,Imax_standard_sev_high(1,:),'r--','LineWidth',1.4);
hold on
l3=plot(Cv,Imax_standard_sev_low(1,:),'b--','Linewidth',1.4);
set(gca,'fontsize',14)
axis([0 Cv(end) 0 pd_max*101])
xlabel('Contact rate (C)')
ylabel({'Peak hospitalization [% of population]'})
legend([l2,l3,l1],{'p_{max}','p_{min}', 'p(H)'},'Fontsize',12,'NumColumns',3,'location','north')






