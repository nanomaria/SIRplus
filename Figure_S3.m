% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

global C delta epsm Dmaxm k sR iR rR pr pi L

epsm = 1; % smoothering factor for function 2
Dmaxm = 6;  %minimal diversity
L = 0.1;
delta = 1*L;  % rate of decay of microbial diversity
pi = 0.5;
pr = 0.8;
d0 = 0;
sR = 1;
iR = 0;
rR = 0;
rp = 0.1; %percent of people recovered
Tfin = 10/L/2;



% Contact rate versus diversity at equilibrium

Cv = 0.5:0.5:50;
kv = 0.25*L:.25*L:1*L;
Nsteps = size(Cv,2);
Msteps = size(kv,2);

Dv = zeros(Msteps,Nsteps);

for j = 1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);
options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y] = ode45(@onlyD_eq, 0:.1:Tfin, d0, options);


figure(4)
subplot(1,2,1)
if k == 0.25*L && (C==3 || C==5 ||C==8 || C==15 || C==30)
w = 2.2/sqrt(sqrt(C));
figure(4)
plot(T,Y,'k-','Linewidth',w)
hold on
axis([0 Tfin 0 7])
ylabel('Interaction-dependent health benefits (H)')
xlabel('Time [arbitrary unit]')
set(gca,'fontsize',14)
end


figure(4)
subplot(1,2,2)
w = 2.2/sqrt(sqrt(C));
if k == 0.75*L && (C==3 || C==5 ||C==8 || C==15 || C==30)
figure(4)
plot(T,Y,'k-','linewidth',w)
hold on
axis([0 Tfin 0 7])
ylabel('Interaction-dependent health benefits (H)')
xlabel('Time [arbitrary unit]')
set(gca,'fontsize',14)
end

Dv(j,i) = Y(end,1);

end
end


