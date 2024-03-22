% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = SIR_square_cumul_sev_eq(t,y)

global pd_min pd_max gamma Hmaxi Hstar choice eta betafix zip


S = y(1);
I = y(2);
Sev = y(3);
dy = zeros(3,1);

if choice ==1
pd = (pd_min + (pd_max-pd_min)/2*(1-tanh((Hstar-Hmaxi))));
end
if choice == 2
%linear function pd(D)
if Hstar<=Hmaxi*2
pd = (-(pd_max-pd_min)/(Hmaxi*2)*Hstar+pd_max);
elseif Hstar>Hmaxi*2
pd = pd_min;
end
end
if zip == 3
pd = pd_min;
end
if zip == 4
pd = pd_max;
end

dy(1) = -betafix*(I+Sev)*S;
dy(2) = (1-pd)*betafix*(I+Sev)*S-gamma*I;
dy(3) = pd*betafix*(I+Sev)*S - eta*Sev;

end