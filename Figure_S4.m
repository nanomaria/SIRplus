% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

global alphai C alphai_min delta gamma epsi epsm Dmaxi Dmaxm pi pr k Dstar sR iR rR

delta = 1;  % rate of decay of microbial diversity
gamma = 1/7;  % recovery rate
alphai = 0.1;  % maximal probability of pathogenic transmission given a contact
alphai_min = 0.02; % minimal probability of pathogenic transmission given a contact
epsi = 1; % smoothering factor for function 1
epsm = 1; % smoothering factor for function 2
Dmaxi = 3;  %maximal diversity
Dmaxm = 6;  %minimal diversity
pi = 0.4;
pr = 0.6;
s0 = 1;
i0 = 1e-6;
r0 = 0;
c0 = 0;
Tfin = 300;

dC = 1;  % 0.1 for smooth curves


Cv = 0:dC:30;
kv = [0.25,0.5,0.75,1];%:0.25:1 ;
Nsteps = size(Cv,2);
Msteps = size(kv,2);
Imax = zeros(Msteps,Nsteps);
Imax_standard = zeros(Msteps,Nsteps);
Cmax = zeros(Msteps,Nsteps);

d0 = 0;
sR= s0;
iR= i0;
rR= r0;





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


Nsteps = size(Cv,2);
Msteps = size(kv,2);
betav1 = zeros(Msteps,Nsteps);
betav2 = zeros(Msteps,Nsteps);
betav = zeros(Msteps,Nsteps);

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


betav2(j,i) = (alphai_min + (alphai-alphai_min)/2*(1-tanh((Dstar(j,i)-Dmaxi)*epsi)));



end




end










figure(3);
%subplot(2,1,1)
fig = figure(3);
right_color =  [0 0.462 0.729];
left_color = [0.995 0.681 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);




yyaxis right
plot(Cv,Dv(4,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1.8)
hold on
line(Cv,Dv(3,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1.4)
hold on
line(Cv,Dv(2,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',1)
hold on
line(Cv,Dv(1,:),'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',0.6)
hold on


set(gca,'fontsize',14)
axis([0 Cv(end) 0 7])
xlabel('Contact rate (C)')
ylabel({'Interaction-dependent health benefits',' at equilibrium (H*)'})



yyaxis left
plot(Cv,betav2(4,:),'Color',[0.995 0.681 0],'Linestyle','-.','Linewidth',1.8)
hold on
line(Cv,betav2(3,:),'Color',[0.995 0.681 0],'Linestyle','-.','Linewidth',1.4)
hold on
line(Cv,betav2(2,:),'Color',[0.995 0.681 0],'Linestyle','-.','Linewidth',1)
hold on
line(Cv,betav2(1,:),'Color',[0.995 0.681 0],'Linestyle','-.','Linewidth',0.6)
hold on
axis([0 Cv(end) -.045 0.25])
hold on
set(gca,'fontsize',14)
ylabel({'Probability of disease transmission (\alpha(H))'})
set(gca,'fontsize',14)
tcy = 0.5; %transparency
p1 = -0.0375;
p2 = -.025;
p3 = -0.0125;
p4 = -0.0;
x = 0.005;
thc = 4; %thickness

d1=line([ Dmaxi*delta/kv(1),Cv(end)],[p1,p1],   'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',thc);
d1.Color(4) = tcy;
hold on
d2=line([ Dmaxi*delta/kv(2),Cv(end)],[p2,p2], 'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',thc);
d2.Color(4) = tcy;
hold on
d3=line([ Dmaxi*delta/kv(3),Cv(end)],[p3,p3],   'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',thc);
d3.Color(4) = tcy;
hold on
d4=line([ Dmaxi*delta/kv(4),Cv(end)],[p4,p4], 'Color',[0.995 0.681 0],'Linestyle','-','LineWidth',thc);
d4.Color(4) = tcy;


b1=line([ 0,Dmaxm*delta/kv(1)],[p1+x,p1+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b1.Color(4) = tcy;
hold on
b2=line([ 0,Dmaxm*delta/kv(2)],[p2+x,p2+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b2.Color(4) = tcy;
hold on
b3=line([ 0,Dmaxm*delta/kv(3)],[p3+x,p3+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b3.Color(4) = tcy;
hold on
b4=line([ 0,Dmaxm*delta/kv(4)],[p4+x,p4+x],'Color',[0 0.462 0.729],'Linestyle','-','LineWidth',thc);
b4.Color(4) = tcy;



