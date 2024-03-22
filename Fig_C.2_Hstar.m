% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

global C delta Hmax h k 

Hmax = 10;  
delta = .1;  % rate of decay of H
h = .1;
H0 = 0;
Tfin = 30; 
k = 0.25;

%H_star at equilibrium

Cv = [3, 5, 8, 15, 30];
kv = [0.25,0.75];
Nsteps = size(Cv,2);
Msteps = size(kv,2);

Hv = zeros(Msteps,Nsteps);

for j = 1:Msteps
    k = kv(j);
for i=1:Nsteps
    C = Cv(i);
[T,Y] = ode45(@onlyD_eq_updated, 0:.001:Tfin, H0);


figure(4)
subplot(1,2,1)
if k == 0.25 && (C==3 || C==5 ||C==8 || C==15 || C==30)
plot(T,Y,'k-','Linewidth',3/sqrt(i))
hold on
axis([0 Tfin 0 11])
ylabel('Interaction-dependent health benefits (H)')
xlabel('Time [arbitrary unit]')
set(gca,'fontsize',14)
end


figure(4)
subplot(1,2,2)
if k == 0.75 && (C==3 || C==5 ||C==8 || C==15 || C==30)
figure(4)
plot(T,Y,'k-','linewidth',3/sqrt(i))
hold on
axis([0 Tfin 0 11])
ylabel('Interaction-dependent health benefits (H)')
xlabel('Time [arbitrary unit]')
set(gca,'fontsize',14)
end

Hv(j,i) = Y(end,1);

end
end


