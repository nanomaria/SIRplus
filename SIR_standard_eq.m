% SIR-square project
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = SIR_standard_eq(t,y)

global beta  C gamma alphai


S = y(1);
I = y(2);
R = y(3);
dy = zeros(3,1);

        
beta = C*(alphai);

dy(1) = -beta*I*S;
dy(2) = beta*I*S-gamma*I;
dy(3) = gamma*I;

end