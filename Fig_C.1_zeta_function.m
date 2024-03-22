% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------

clear all
close all

Hmax = 10;  
h = .1;
C = 10; 
k = 0.5;


Hv = 0:0.1:15;

zeta = C*k*h*1/2*(1-tanh(3*(Hv-Hmax)));

figure(1);
plot(Hv,zeta,'-k','LineWidth',1.5);
axis([0 Hv(end) 0 0.7])
xlabel('H')
ylabel('\zeta(H)')
set(gca,'fontsize',14)


