% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = Eq_fig3_updated(t,y)

global pd_min pd_max alphai alphai_min gamma eta Hmaxi Hstar beta zip C


S = y(1);
I = y(2);
Sev = y(3);
dy = zeros(3,1);


beta = C*(alphai_min + (alphai-alphai_min)/2*(1-tanh((Hstar-Hmaxi))));


if zip == 1
pd = pd_max;
end
if zip == 2
pd = (pd_min + (pd_max-pd_min)/2*(1-tanh((Hstar-Hmaxi))));
end

dy(1) = -beta*(I+Sev)*S;
dy(2) = (1-pd)*beta*(I+Sev)*S-gamma*I;
dy(3) = pd*beta*(I+Sev)*S - eta*Sev;

end