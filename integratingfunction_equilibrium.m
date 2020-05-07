function dxdt = integratingfunction_equilibrium(t,x,params,Init)
% See https://www.mathworks.com/help/matlab/math/ode-event-location.html
% for more details.
%%%%%%%%%%%%%%%%%%%%
%
% Name: integratingfunction_equilibrium.m
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
% Specifically, this function is used to simulate equilibration of PD-1 and 
% PD-L1 levels on the PD1-PDL1 pathway (untreated case).
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
%
%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%Proteins:
%PD1=x(1); PD-L1=x(2); PD-1:PD-L1=x(3); 

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
% Read in parameters
%--------------------------------------------------------------------------
% complex association/disassociation

beta_plus_tilde = params(1); beta_minus_tilde = params(2); 
%--------------------------------------------------------------------------
PD1_0 = Init(1); PDL1_0 = Init(2);

%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);
%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------

% Normalized to PD1_0               
%(1) PD-1: 
dxdt(1) = -beta_plus_tilde*PDL1_0*x(1)*x(2) + beta_minus_tilde*x(3);
            
% Normalized to PDL1_0
%(2) PD-L1: 
dxdt(2) = -beta_plus_tilde*PD1_0*x(1)*x(2) + beta_minus_tilde*PD1_0*x(3)/PDL1_0;
 
% Normalized to PD1_0
%(3) PD-1:PD-L1: 
dxdt(3) = beta_plus_tilde*PDL1_0*x(1)*x(2) - beta_minus_tilde*x(3);
            
end
