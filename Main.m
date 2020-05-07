%%%%%%%%%%%%%%%%%%%%
%
% Name: Main.m
% Code authors: Michelle Przedborski, Moriah Pellowe
% Last updated: May 7, 2020
%
% This script will simulate the ex vivo experimental treatment protocols on
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
% Input: parameter set as described below (parameters with values >= 0)
% Output: (1) cytokine expression levels during PD-1 blockade cytokine 
%         experiments are contained in x_treatment_cyto array
%         (2) T-cell populations during PD-1 blockade flow cytometry
%         experiments are contained in x_treatment_FC array
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

function Main
t_scale = 12*60;  %Time scale

%Experimental parameters:
total_cells_cyto = 1e5;   %Cell numbers
total_cells_flow = 2e6;   %Cell numbers
nivo_0 = 132e-6;     %pg/ml

%--------------------------------------------------------------------------
% Scaled parameter values: (parameters described below)
%--------------------------------------------------------------------------
pars = [0.448, 1.537, 3339.162, 0.114,...
        0.1266, 0.01769, 0.02389, 0.02024,...
        0.003494, 0.01886, 0.01766, 0.01803,...
        0.006100, 0.09426, 0.01790, 0.01047,...
        0.01134, 0.000005470, 0.0006529, 0.00007885,...
        0.009159, 0.0006046, 0.4991, 0.5075,...
        0.4994, 0.3485, 0.02405, 159.05,...
        0.4049, 0.6307, 4.0322, 0.8436,...
        130.04, 0.006307, 0.03447, 0.08895,...
        0.7467, 141.94, 7127.55, 7198.54,...
        0.1266, 0.0008437, 0.7585, 0.001736,...
        48.515, 2.1451, 19.433, 0.8121,...
        0.0003879, 0.0008561];	

%--------------------------------------------------------------------------
% Initial cytokine levels (for cytokine experiments)
%--------------------------------------------------------------------------
IFNg_0 = pars(1); IL12_0 = pars(2); IL6_0 = pars(3); IL4_0 = pars(4);

%--------------------------------------------------------------------------
% Optimization parameters
%--------------------------------------------------------------------------
% net proliferation
n_4_tilde = pars(5); n_8_tilde = pars(6); n_1_tilde = pars(7); 
n_C_tilde = pars(8); n_Can_tilde = pars(9);
% growth and decay of cells
g_2_tilde = pars(10); g_2_4_tilde = pars(11); g_C_12_tilde = pars(12); 
delta_2_tilde = pars(13);
% differentiation
d_1_IFN_tilde = pars(14); d_1_12_tilde = pars(15); d_2_tilde = pars(16); 
d_C_tilde = pars(17);
% cell killing
k_C_tilde = pars(18);
% production
p_1_IFN_tilde = pars(19)/IFNg_0; p_2_4_6_tilde = pars(20)/IL4_0; 
p_Can_6_tilde = pars(21)/IL6_0; p_Can_12_tilde = pars(22)/(24*60*IL12_0);
% decay
delta_IFN_tilde = pars(23); delta_IL4_tilde = pars(24); 
delta_IL6_tilde = pars(25); delta_IL12_tilde = pars(26); 
delta_A_tilde = pars(27);
% induction/upregulation
q_1 = pars(28); q_IFN_1_tilde = pars(29)/IFNg_0; 
q_IFN_PDL1_tilde = pars(30)/IFNg_0; q_gIL4_tilde = pars(31)/IL4_0; 
q_dIL4_tilde = pars(32)/IL4_0; q_IL6_tilde = pars(33)/IL6_0; 
q_dIL12_tilde = pars(34)/IL12_0; q_gIL12_tilde = pars(35)/IL12_0;
% inhibition
r_IFN_tilde = pars(36)/IFNg_0; r_IL4_tilde = pars(37)/IL4_0; 
r_IL6_tilde = pars(38)/IL6_0;
% protein expression
rho = pars(39)/t_scale; lambda = pars(40)/t_scale; lambda_Can_IFN = pars(41)/t_scale;
% complex association/disassociation
beta_plus_tilde = pars(42); beta_minus_tilde = pars(43); 
alpha_plus_tilde = pars(42); %assume that alpha_plus and beta_plus are the 
                             %same and only the dissociation rate changes
alpha_minus_tilde = pars(44);
% inhibition by PD-1:PD-L1
s_1 = pars(45); s_2 = pars(46); s_C = pars(47);
% cell population fractions
C_frac = pars(48); Th1_frac = pars(49); Th2_frac = pars(50);

%--------------------------------------------------------------------------
% Initialize remaning variables
%--------------------------------------------------------------------------
% Experimental T-cell fractions at t=0 (assume relative populations are the
% same as for the flow cytometry experiments)
TN8_frac = 0.6471;
Tc_frac = 0.1048;
CD4_frac = 0.2481;

% 1. Initialize cell populations for cytokine experiments
%--------------------------------------------------------------------------
T_cell_0 = total_cells_cyto;
TN4_0 = (1.0-Th1_frac-Th2_frac)*CD4_frac*T_cell_0;
Th1_0 = Th1_frac*CD4_frac*T_cell_0;
Th2_0 = Th2_frac*CD4_frac*T_cell_0;
TN8_0 = TN8_frac*T_cell_0;
Tc_0 = Tc_frac*T_cell_0;

% 2. Calculate remaining protein parameters:
%--------------------------------------------------------------------------
p_1_12_tilde = delta_IL12_tilde/Th1_0;
p_2_6_tilde = delta_IL6_tilde/Th2_0;
p_2_4_tilde = (delta_IL4_tilde-p_2_4_6_tilde*Th2_0/(q_IL6_tilde+1.0))/Th2_0;
p_C_IFN_tilde = (delta_IFN_tilde-p_1_IFN_tilde*Th1_0...
              *r_IL4_tilde/(r_IL4_tilde + 1.0)*r_IL6_tilde/(r_IL6_tilde + 1.0))/Tc_0;                             
          
% 3. Initialize the total amount of PD-1 and PD-L1
%--------------------------------------------------------------------------
PD1_0 = rho*(Th1_0 + Th2_0 + Tc_0);
rho_tilde = rho/PD1_0;
s_1_tilde = s_1/PD1_0;
s_2_tilde = s_2/PD1_0;
s_C_tilde = s_C/PD1_0;

PDL1_0 = lambda*(Th1_0 + Th2_0 + Tc_0);
lambda_tilde = lambda/PDL1_0;
lambda_Can_IFN_tilde = lambda_Can_IFN/PDL1_0;
     
%==========================================================================
% 1. Simulations to equilibrate PD-1:PD-L1, PD-1, and PD-L1
%==========================================================================
IC_0 = [1.0, 1.0, 0];   
IC_0 = IC_0.';

DE_params_equil = [beta_plus_tilde, beta_minus_tilde];

Init = [PD1_0, PDL1_0];

% An event function is included here to ensure the system reaches an 
% equilibrium state.  
%--------------------------------------------------------------------------
tmin = 0.0;
tmax_equil = 1000000/t_scale;
num_equil_points = 50000;
tspan_equil=linspace(tmin,tmax_equil,num_equil_points);

% Pass the parameters to the event function first:
%--------------------------------------------------------------------------
ParamEventHdl = @(t_equil, x_equil)EquilibriumEventFcn(t_equil,x_equil,DE_params_equil,Init,t_scale); 
opts_equil = odeset('RelTol',1e-10,'AbsTol',1e-12,'Events',ParamEventHdl); % Add the events function here. 

% Run the equilibration simulation: (Call ode45 solver while passing the 
% parameters and ODE options to the integrating function)
%--------------------------------------------------------------------------
[t_equil,x_equil,te,ye,ie] = ode45(@(t_equil,x_equil)...
    integratingfunction_equilibrium(t_equil,x_equil,DE_params_equil,Init),tspan_equil,IC_0,opts_equil);

% Check to make sure that the levels did actually equilibrate.
%--------------------------------------------------------------------------
if max(t_equil)==tmax_equil
   disp('ERROR: PD-1:PD-L1 levels did not equilibrate.');
   return
end 
          
%==========================================================================
% 2. Cytokine production in treated cells
%==========================================================================
% Parameter array to pass to integrating function:
%--------------------------------------------------------------------------
DE_params = [n_4_tilde, n_8_tilde, n_1_tilde, n_C_tilde, n_Can_tilde,...
             g_2_tilde, g_2_4_tilde, g_C_12_tilde,... 
             delta_2_tilde,... 
             d_1_IFN_tilde, d_1_12_tilde, d_2_tilde, d_C_tilde,...
             k_C_tilde,...
             p_1_IFN_tilde, p_1_12_tilde, p_2_4_tilde, p_2_4_6_tilde, p_2_6_tilde,...
             p_C_IFN_tilde, p_Can_6_tilde, p_Can_12_tilde,...
             delta_IFN_tilde, delta_IL4_tilde, delta_IL6_tilde,...
             delta_IL12_tilde, delta_A_tilde,...
             q_1, q_IFN_1_tilde, q_IFN_PDL1_tilde, q_gIL4_tilde,...
             q_dIL4_tilde, q_IL6_tilde, q_dIL12_tilde, q_gIL12_tilde,...
             r_IFN_tilde, r_IL4_tilde, r_IL6_tilde,...
             rho_tilde, lambda_tilde, lambda_Can_IFN_tilde,...
             beta_plus_tilde, beta_minus_tilde,...  
             alpha_plus_tilde, alpha_minus_tilde,...
             s_1_tilde, s_2_tilde, s_C_tilde];

Init = [PD1_0, PDL1_0, nivo_0];

tmin = 0.0;
tmax = 24*60/t_scale;
num_points = 8000;
tspan=linspace(tmin,tmax,num_points);         

%--------------------------------------------------------------------------         
% i) First 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = [TN4_0, TN8_0, Th1_0, Th2_0, Tc_0,...
      1.0, 1.0, 1.0, 1.0,...
      x_equil(end,1), x_equil(end,2), x_equil(end,3),...
      1.0, 0];  
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 
[t_treatment_1d,x_treatment_1d] = ode15s(@(t_treatment_1d,x_treatment_1d)...
       integratingfunction_cytokine_treatment(t_treatment_1d,x_treatment_1d,DE_params,Init),tspan,IC,opts);
   
%--------------------------------------------------------------------------         
% ii) Second 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_1d(end,:);
IC(13) = 1.0;   %Wash out and replenish the drug
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 

[t_treatment_2d,x_treatment_2d] = ode15s(@(t_treatment_2d,x_treatment_2d)...
       integratingfunction_cytokine_treatment(t_treatment_2d,x_treatment_2d,DE_params,Init),tspan,IC,opts);

t_treatment_2d = t_treatment_2d + 24*60/t_scale;

%--------------------------------------------------------------------------         
% iii) Third 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_2d(end,:);
IC(13) = 1.0;    %Wash out and replenish the drug
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 

[t_treatment_3d,x_treatment_3d] = ode15s(@(t_treatment_3d,x_treatment_3d)...
       integratingfunction_cytokine_treatment(t_treatment_3d,x_treatment_3d,DE_params,Init),tspan,IC,opts);
   
t_treatment_3d = t_treatment_3d + 2*24*60/t_scale;

% Append all the treatment arrays
%--------------------------------------------------------------------------
t_treatment_cyto = [t_treatment_1d; t_treatment_2d; t_treatment_3d];
x_treatment_cyto = [x_treatment_1d; x_treatment_2d; x_treatment_3d];

%==========================================================================  
% 3. Flow cytometry experiments
%==========================================================================
% We assume the same production and decay rates for the flow cytometry
% experiments. Since the number of cells will be different, the initial
% protein levels will be different. Thus we re-calculate the initial
% protein levels and new scaled parameters and proceed as before.
C_0 = C_frac*total_cells_flow;

T_cell_0 = (1.0-C_frac)*total_cells_flow;
TN4_0 = (1.0-Th1_frac-Th2_frac)*CD4_frac*T_cell_0;
Th1_0 = Th1_frac*CD4_frac*T_cell_0;
Th2_0 = Th2_frac*CD4_frac*T_cell_0;
TN8_0 = TN8_frac*T_cell_0;
Tc_0 = Tc_frac*T_cell_0;

%--------------------------------------------------------------------------
% Calculate new SS initial protein levels based on cell populations, then
% calculate new scaled parameters:
p_Can_12 = p_Can_12_tilde*IL12_0;
p_1_12 = p_1_12_tilde*IL12_0;

IL12_0_FC = (p_Can_12*C_0+p_1_12*Th1_0)/delta_IL12_tilde;
%---------------------------------
p_2_6 = p_2_6_tilde*IL6_0;
p_Can_6 = p_Can_6_tilde*IL6_0;

IL6_0_FC = (p_2_6*Th2_0 + p_Can_6*C_0)/delta_IL6_tilde;
%---------------------------------
p_2_4 = p_2_4_tilde*IL4_0;
p_2_4_6 = p_2_4_6_tilde*IL4_0;
q_IL6 = q_IL6_tilde*IL6_0;

IL4_0_FC = (p_2_4*Th2_0 + p_2_4_6*Th2_0*IL6_0_FC/(q_IL6+IL6_0_FC))/delta_IL4_tilde;
%---------------------------------
p_1_IFN = p_1_IFN_tilde*IFNg_0;
r_IL4 = r_IL4_tilde*IL4_0;
r_IL6 = r_IL6_tilde*IL6_0;
p_C_IFN = p_C_IFN_tilde*IFNg_0;

IFNg_0_FC = (p_1_IFN*Th1_0*r_IL4/(r_IL4 + IL4_0_FC)*r_IL6/(r_IL6 + IL6_0_FC)...
         + p_C_IFN*Tc_0)/delta_IFN_tilde;
%---------------------------------
PD1_0_FC = rho*(Th1_0 + Th2_0 + Tc_0);
rho_tilde_FC = rho/PD1_0_FC;
s_1_tilde_FC = s_1/PD1_0_FC;
s_2_tilde_FC = s_2/PD1_0_FC;
s_C_tilde_FC = s_C/PD1_0_FC;
%---------------------------------
q_IFN_PDL1 = q_IFN_PDL1_tilde*IFNg_0;

PDL1_0_FC = lambda*(Th1_0 + Th2_0 + Tc_0 + C_0) + ...
            lambda_Can_IFN*C_0*IFNg_0_FC/(q_IFN_PDL1 + IFNg_0_FC);

lambda_tilde_FC = lambda/PDL1_0_FC;
lambda_Can_IFN_tilde_FC = lambda_Can_IFN/PDL1_0_FC;

%--------------------------------------------------------------------------
%Now rescale the parameters appropriately
p_1_IFN_tilde_FC = p_1_IFN_tilde*IFNg_0/IFNg_0_FC;
p_1_12_tilde_FC = p_1_12_tilde*IL12_0/IL12_0_FC;
p_2_4_tilde_FC = p_2_4_tilde*IL4_0/IL4_0_FC;
p_2_6_tilde_FC = p_2_6_tilde*IL6_0/IL6_0_FC;
p_2_4_6_tilde_FC = p_2_4_6_tilde*IL4_0/IL4_0_FC; 
p_C_IFN_tilde_FC = p_C_IFN_tilde*IFNg_0/IFNg_0_FC;
p_Can_6_tilde_FC = p_Can_6_tilde*IL6_0/IL6_0_FC; 
p_Can_12_tilde_FC = p_Can_12_tilde*IL12_0/IL12_0_FC;
q_IFN_1_tilde_FC = q_IFN_1_tilde*IFNg_0/IFNg_0_FC; 
q_IFN_PDL1_tilde_FC = q_IFN_PDL1_tilde*IFNg_0/IFNg_0_FC; 
q_gIL4_tilde_FC = q_gIL4_tilde*IL4_0/IL4_0_FC; 
q_dIL4_tilde_FC = q_dIL4_tilde*IL4_0/IL4_0_FC; 
q_IL6_tilde_FC = q_IL6_tilde*IL6_0/IL6_0_FC; 
q_dIL12_tilde_FC = q_dIL12_tilde*IL12_0/IL12_0_FC; 
q_gIL12_tilde_FC = q_gIL12_tilde*IL12_0/IL12_0_FC;
r_IFN_tilde_FC = r_IFN_tilde*IFNg_0/IFNg_0_FC; 
r_IL4_tilde_FC = r_IL4_tilde*IL4_0/IL4_0_FC; 
r_IL6_tilde_FC = r_IL6_tilde*IL6_0/IL6_0_FC;


%==========================================================================
% 1. Simulations to equilibrate PD-1:PD-L1, PD-1, and PD-L1
%==========================================================================
IC_0 = [1.0, 1.0, 0];   
IC_0 = IC_0.';

DE_params_equil = [beta_plus_tilde, beta_minus_tilde];

Init = [PD1_0_FC, PDL1_0_FC];

% An event function is included here to ensure the system reaches an 
% equilibrium state.  
%--------------------------------------------------------------------------
tmin = 0.0;
tmax_equil = 1000000/t_scale;
num_equil_points = 50000;
tspan_equil=linspace(tmin,tmax_equil,num_equil_points);

% Pass the parameters to the event function first:
%--------------------------------------------------------------------------
ParamEventHdl = @(t_equil, x_equil)EquilibriumEventFcn(t_equil,x_equil,DE_params_equil,Init,t_scale); 
opts_equil = odeset('RelTol',1e-10,'AbsTol',1e-12,'Events',ParamEventHdl); % Add the events function here. 

% Run the equilibration simulation: (Call ode45 solver while passing the 
% parameters and ODE options to the integrating function)
%--------------------------------------------------------------------------
[t_equil,x_equil,te,ye,ie] = ode45(@(t_equil,x_equil)...
    integratingfunction_equilibrium(t_equil,x_equil,DE_params_equil,Init),tspan_equil,IC_0,opts_equil);
          
%==========================================================================
% 2. Flow cytometry in treated cells
%==========================================================================
% Parameter array to pass to integrating function:
%--------------------------------------------------------------------------
DE_params = [n_4_tilde, n_8_tilde, n_1_tilde, n_C_tilde, n_Can_tilde,...
             g_2_tilde, g_2_4_tilde, g_C_12_tilde,... 
             delta_2_tilde,... 
             d_1_IFN_tilde, d_1_12_tilde, d_2_tilde, d_C_tilde,...
             k_C_tilde,...
             p_1_IFN_tilde_FC, p_1_12_tilde_FC, p_2_4_tilde_FC,...
             p_2_4_6_tilde_FC, p_2_6_tilde_FC,...
             p_C_IFN_tilde_FC, p_Can_6_tilde_FC, p_Can_12_tilde_FC,...
             delta_IFN_tilde, delta_IL4_tilde, delta_IL6_tilde,...
             delta_IL12_tilde, delta_A_tilde,...
             q_1, q_IFN_1_tilde_FC, q_IFN_PDL1_tilde_FC,...
             q_gIL4_tilde_FC, q_dIL4_tilde_FC, q_IL6_tilde_FC,...
             q_dIL12_tilde_FC, q_gIL12_tilde_FC,...
             r_IFN_tilde_FC, r_IL4_tilde_FC, r_IL6_tilde_FC,...
             rho_tilde_FC, lambda_tilde_FC, lambda_Can_IFN_tilde_FC,...
             beta_plus_tilde, beta_minus_tilde,...  
             alpha_plus_tilde, alpha_minus_tilde,...
             s_1_tilde_FC, s_2_tilde_FC, s_C_tilde_FC];

Init = [PD1_0_FC, PDL1_0_FC, nivo_0];

tmin = 0.0;
tmax = 24*60/t_scale;
num_points = 8000;
tspan=linspace(tmin,tmax,num_points);         

%--------------------------------------------------------------------------         
% i) First 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = [TN4_0, TN8_0, Th1_0, Th2_0, Tc_0, C_0,...
      1.0, 1.0, 1.0, 1.0,...
      x_equil(end,1), x_equil(end,2), x_equil(end,3),...
      1.0, 0];  
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 

[t_treatment_1d,x_treatment_1d] = ode15s(@(t_treatment_1d,x_treatment_1d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_1d,x_treatment_1d,DE_params,Init),tspan,IC,opts);
   
%--------------------------------------------------------------------------         
% ii) Second 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_1d(end,:);
IC(14) = 1.0;   %Wash out and replenish the drug
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 

[t_treatment_2d,x_treatment_2d] = ode15s(@(t_treatment_2d,x_treatment_2d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_2d,x_treatment_2d,DE_params,Init),tspan,IC,opts);
   
t_treatment_2d = t_treatment_2d + 24*60/t_scale;

%--------------------------------------------------------------------------         
% iii) Third 24 hours:         
%--------------------------------------------------------------------------
% Initial conditions array to pass to integrating function:
%--------------------------------------------------------------------------
IC = x_treatment_2d(end,:);
IC(14) = 1.0;    %Wash out and replenish the drug
IC = IC.';

% Run the treatment simulations
%--------------------------------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12); 

[t_treatment_3d,x_treatment_3d] = ode15s(@(t_treatment_3d,x_treatment_3d)...
       integratingfunction_flow_cytometry_treatment(t_treatment_3d,x_treatment_3d,DE_params,Init),tspan,IC,opts);
   
t_treatment_3d = t_treatment_3d + 2*24*60/t_scale;

% Append all the treatment arrays
%--------------------------------------------------------------------------
t_treatment_FC = [t_treatment_1d; t_treatment_2d; t_treatment_3d];
x_treatment_FC = [x_treatment_1d; x_treatment_2d; x_treatment_3d];

end