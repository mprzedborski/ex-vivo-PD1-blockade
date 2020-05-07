function [position,isterminal,direction] = EquilibriumEventFcn(t,x,DE_params_equil,Init,t_scale)
% See https://www.mathworks.com/help/matlab/math/ode-event-location.html
% for more details.
%%%%%%%%%%%%%%%%%%%%
%
% Name: EquilibriumEventFcn.m
% Code authors: Michelle Przedborski, Moriah Pellowe
% Last updated: May 7, 2020
%
% This script is part of the ex vivo simulated treatment protocol on
% the PD-1 network described in the manuscript:
% "Integrating an ex vivo human tumor model with systems biology for the 
% study of PD-1 blockade response dynamics in head and neck squamous cell 
% carcinoma (HNSCC).", authored by:  
% Munisha Smalley, Michelle Przedborski, Saravanan Thiyagaragan, Moriah 
% Pellowe, Amit Verma, Nilesh Brijwani, Debika Datta, Misti Jain, 
% Basavaraja U. Shanthappa, Vidushi Kapoor, Kodaganur S. Gopinath, D.C. 
% Doval, K.S. Sabitha, Gaspar Taroncher-Oldenburg, Biswanath Majumder, 
% Pradip Majumder, Mohammad Kohandel, and Aaron Goldman   
%
% Specifically, this function is used to find equilibrium PD1, PD-L1 and 
% PD1:PD-L1 complex levels in the untreated case.
%
% Note: All equations and parameters have been scaled. In particular, time
%       is measured relative to t_scale (12 hours) and all protein levels
%       are measured relative to initial levels (except the protein 
%       complexes).
%       To convert the parameter values that are here to the same units as
%       stated in the manuscript, use the following conversions:
%       -to get to units of " /day", multiply by 2; 
%       -to get to units of " /min", multiply by (2/(24*60)), i.e. divide
%       by (12*60)
%
%%%%%%%%%%%%%%%%%%%%
isterminal = 1;
direction = [];
threshold = 0.005;

%--------------------------------------------------------------------------
dxdt = integratingfunction_equilibrium(t,x,DE_params_equil,Init);
%--------------------------------------------------------------------------
Init = [Init(1); Init(2); Init(1)];
derivatives_vec=Init.*dxdt/t_scale;
norm_vec=norm(derivatives_vec);

% Logical(cond) evaluates to 1 if cond is true and to 0 if cond is false.
% So we want to make the condition the opposite of what we actually want,
% e.g. check if norm_vec is above the threshold. Then when it is below,
% position will drop to zero and the integration will stop.
if norm_vec > threshold
   position = 1; 
else
   position = 0;
end

end

