function dxdt = integratingfunction_flow_cytometry_treatment(t,x,params,Init)
% See https://www.mathworks.com/help/matlab/math/ode-event-location.html
% for more details.
%%%%%%%%%%%%%%%%%%%%
%
% Name: integratingfunction_flow_cytometry_treatment.m
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
% Specifically, this function is used to simulate the flow cytometry 
% experiments.
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

%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%Cell populations:
%TN4=x(1); TN8=x(2); Th1=x(3); Th2=x(4); Tc=x(5); C=x(6);
%--------------------------------------------------------------------------
%Proteins:
%IFN_gamma=x(7); IL-4=x(8); IL-6=x(9); IL-12=x(10); PD1=x(11); 
%PD-L1=x(12); PD-1:PD-L1=x(13); 
%--------------------------------------------------------------------------
%Drugs:
%A=x(14); A:PD-1=x(15); 

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
% Read in all parameters
%--------------------------------------------------------------------------
% net proliferation
n_4_tilde = params(1); n_8_tilde = params(2); n_1_tilde = params(3); 
n_C_tilde = params(4); n_Can_tilde = params(5);
% growth and decay of cells
g_2_tilde = params(6); g_2_4_tilde = params(7); g_C_12_tilde = params(8); 
delta_2_tilde = params(9);
% differentiation
d_1_IFN_tilde = params(10); d_1_12_tilde = params(11); d_2_tilde = params(12); 
d_C_tilde = params(13);
% cell killing
k_C_tilde = params(14);
% production
p_1_IFN_tilde = params(15); p_1_12_tilde = params(16); 
p_2_4_tilde = params(17); p_2_4_6_tilde = params(18); 
p_2_6_tilde = params(19); p_C_IFN_tilde = params(20); 
p_Can_6_tilde = params(21); p_Can_12_tilde = params(22);
% decay
delta_IFN_tilde = params(23); delta_IL4_tilde = params(24); 
delta_IL6_tilde = params(25); delta_IL12_tilde = params(26); 
delta_A_tilde = params(27);
% induction/upregulation
q_1 = params(28); q_IFN_1_tilde = params(29); 
q_IFN_PDL1_tilde = params(30); 
q_gIL4_tilde = params(31); q_dIL4_tilde = params(32); q_IL6_tilde = params(33); 
q_dIL12_tilde = params(34); q_gIL12_tilde = params(35);
% inhibition
r_IFN_tilde = params(36); r_IL4_tilde = params(37); r_IL6_tilde = params(38);
% protein expression
rho_tilde = params(39); lambda_tilde = params(40); 
lambda_Can_IFN_tilde = params(41);
% complex association/disassociation
beta_plus_tilde = params(42); beta_minus_tilde = params(43); 
alpha_plus_tilde = params(44); alpha_minus_tilde = params(45);
% inhibition by PD-1:PD-L1
s_1_tilde = params(46); s_2_tilde = params(47); s_C_tilde = params(48);
%--------------------------------------------------------------------------

PD1_0 = Init(1); PDL1_0 = Init(2); A_0 = Init(3);

%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);
%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------   
%(1) T_N4: 
dxdt(1) = n_4_tilde*x(1)...
        - (d_1_12_tilde*x(1) * (x(10)/(q_dIL12_tilde + x(10)))...
        +  d_1_IFN_tilde*x(1) * (x(7)/(q_IFN_1_tilde + x(7))))*(s_1_tilde/(s_1_tilde+x(13)))...
        - (d_2_tilde*x(1) * (x(8)/(q_dIL4_tilde + x(8))))*(s_2_tilde/(s_2_tilde+x(13)));
            
%(2) T_N8: 
dxdt(2) = n_8_tilde*x(2) - d_C_tilde*x(2) * (x(3)/(q_1+x(3)))*(s_C_tilde/(s_C_tilde + x(13)));
            
%(3) Th1: 
dxdt(3) = n_1_tilde*x(3) + (d_1_12_tilde*x(1) * (x(10)/(q_dIL12_tilde + x(10))) + ...
          d_1_IFN_tilde*x(1) * (x(7)/(q_IFN_1_tilde + x(7))))*(s_1_tilde/(s_1_tilde+x(13)));
            
%(4) Th2: 
dxdt(4) = (g_2_tilde*x(4) + g_2_4_tilde*x(4) * (x(8)/(q_gIL4_tilde + x(8)))) * ...
          (r_IFN_tilde/(r_IFN_tilde + x(7))) + ...
          (d_2_tilde*x(1) * (x(8)/(q_dIL4_tilde + x(8))))*(s_2_tilde/(s_2_tilde+x(13))) - ...
           delta_2_tilde*x(4);
            
%(5) Tc: 
dxdt(5) = n_C_tilde*x(5) + g_C_12_tilde*x(5)*x(10)/(q_gIL12_tilde + x(10))...
        + d_C_tilde*x(2) * (x(3)/(q_1+x(3)))*(s_C_tilde/(s_C_tilde + x(13)));
    
%(6) Cancer cells:
dxdt(6) = n_Can_tilde*x(6) - k_C_tilde*x(6)*x(5);
    
% Normalized to IFNg_0    
%(7) IFN_gamma: 
dxdt(7) = p_1_IFN_tilde*x(3) * (r_IL4_tilde/(r_IL4_tilde + x(8))) * (r_IL6_tilde/(r_IL6_tilde + x(9))) + ...
          p_C_IFN_tilde*x(5) - delta_IFN_tilde*x(7);
            
% Normalized to IL4_0
%(8) IL-4: 
dxdt(8) = p_2_4_tilde*x(4) + p_2_4_6_tilde*x(4) * (x(9)/(q_IL6_tilde + x(9))) - ...
          delta_IL4_tilde*x(8);
            
% Normalized to IL6_0
%(9) IL-6: 
dxdt(9) = p_2_6_tilde*x(4) + p_Can_6_tilde*x(6) - delta_IL6_tilde*x(9);

% Normalized to IL12_0            
%(10) IL-12: 
dxdt(10) = p_1_12_tilde*x(3) + p_Can_12_tilde*x(6) - delta_IL12_tilde*x(10);

% Normalized to PD1_0        
%(11) PD-1: 
dxdt(11) =  rho_tilde*(dxdt(3) + dxdt(4) + dxdt(5))...   
            - beta_plus_tilde*PDL1_0*x(11)*x(12) + beta_minus_tilde*x(13)...
            - alpha_plus_tilde*A_0*x(11)*x(14) + alpha_minus_tilde*A_0*x(15)/PD1_0;
            
% Normalized to PDL1_0   
%(12) PD-L1: 
dxdt(12) =  lambda_tilde*(dxdt(3) + dxdt(4) + dxdt(5) + dxdt(6))...
            + lambda_Can_IFN_tilde*dxdt(6)*(x(7)/(q_IFN_PDL1_tilde + x(7)))...        
            - beta_plus_tilde*PD1_0*x(11)*x(12) + beta_minus_tilde*PD1_0*x(13)/PDL1_0;
            
% Normalized to PD1_0
%(13) PD-1:PD-L1: 
dxdt(13) =  beta_plus_tilde*PDL1_0*x(11)*x(12) - beta_minus_tilde*x(13);
            
% Normalized to A_0
%(14) A (free drug): 
dxdt(14) = - alpha_plus_tilde*PD1_0*x(14)*x(11) + alpha_minus_tilde*x(15)...
           - delta_A_tilde*x(14);
            
% Normalized to A_0
%(15) A:PD-1: 
dxdt(15) = alpha_plus_tilde*PD1_0*x(14)*x(11) - alpha_minus_tilde*x(15);         

end
