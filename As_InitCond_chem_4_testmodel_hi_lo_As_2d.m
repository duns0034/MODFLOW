function [mob_eq_comp, mob_eq_sw_ic, mob_eq_pw_ic, catex_comp, catex_pw_ic, surf_sp_comp, surf_sp_pw_ic] = ...
    As_InitCond_chem_4_testmodel_hi_lo_As_test4(phrq_sim_dir, phrq_exe, tempC, fl_Yexch) 

% _testmodel_hi_lo_As_2D:  Modified to be amenable to 2D MODFLOW scripts
% for GW flow.
%
% ############################################################
% Search for: "CUSTOMIZE TO YOUR COMPUTER!!"
% ############################################################

phrq_infil = fullfile(phrq_sim_dir, 'batch_phrq.dat'); 
phrq_outfil = fullfile(phrq_sim_dir, 'batch_phrq.out'); 
phrq_tblfil = fullfile(phrq_sim_dir, 'batch_phrq.tbl'); 
phrq_sel_outfil = fullfile(phrq_sim_dir, 'batch_phrq.sel'); 

file_databas = fullfile(phrq_sim_dir, 'pht3d_datab.dat');

por = 0.20;
tempC = 10.0;
units = 'mol/kgw';
phrq_sim_dir = sim_dir;

% *** CUSTOMIZE: **********************************************************
% This script uses PHREEQC to prepare (equilibrate and charge-balance)
% different aqueous solutions for use in initial and boundary conditions.  
% Different aq solutions could correspond to different parts of the domain,
% e.g. Cell 2 water, bulk domain groundwater, and recharge water.  These
% different aq solutions will be called 'zones'.
%  
% Use this section to specify the number of zones and to comment on what
% each zone represents.

nzones = 2;
% ** comment here to describe 'zones'
% zone 1: Aquitard water, boundary  and initial concentrations
% zone 2: Aquifer water, boundary and initial concentrations
zone_names = {'Aquitard', 'Aquifer'};
%**************************************************************************

% -- SET MODEL CHEMISTRY
% (Units -- aq: mol/L_w, user-defined immob (e.g. bacteria, napl): mol/L_w, 
% minerals (and gases?): mol/L_v, exchangers and surfaces: mol/L_v)

% - mobile equil components
% (ic from PHREEQC batch run)    
n_mob_eq_max = 50;
mob_eq_comp = cell(n_mob_eq_max,1);

% direct sw and pw obs from spreadsheet (units in spreadsheet)
mob_eq_obs_z = zeros(nzones, n_mob_eq_max);
convert_obs2model = zeros(1, n_mob_eq_max); % convert obs units to mol/L
mob_eq_extra_z = cell(nzones, n_mob_eq_max);
ii = 0; 

%**************************************************************************

ii=ii+1; mob_eq_comp{ii} = 'As(3)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L 
convert_obs2model(ii) = 1/((1e6)*74.992);
ii=ii+1; mob_eq_comp{ii} = 'As(5)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % ug/L 
convert_obs2model(ii) = 1/((1e6)*74.992);
ii=ii+1; mob_eq_comp{ii} = 'C(4)'; % **************************************
zz = 0;
Ka1 = 10^-6.3; Ka2 = 10^-10.3; KH = 10^-14;
pH = 7;
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(421.5, pH, Ka1, Ka2, KH); % meq/L (alk = 207.5 mg/L as CaCO3 at Cell 1, 307.5 mg/L as CaCO3 at Cell 2) 
zz = zz + 1; mob_eq_obs_z(zz,ii) = Alk2DIC(412.8, pH, Ka1, Ka2, KH); % eq/L (alk = 150 mg/L as CaCO3 at MW12)
convert_obs2model(ii) = 1*0.6/((1e3)*100.089);
ii=ii+1; mob_eq_comp{ii} = 'C(-4)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L  
convert_obs2model(ii) = 1/1e3/16.0;
ii=ii+1; mob_eq_comp{ii} = 'Ca'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 139.1; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 112.98; % mg/L
convert_obs2model(ii) = 1/((1e3)*40.078);
ii=ii+1; mob_eq_comp{ii} = 'Cl'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 2.325; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 23.436; % mg/L 
convert_obs2model(ii) = 1/((1e3)*35.453);
ii=ii+1; mob_eq_comp{ii} = 'Fe(2)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 2.855; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/((1e3)*55.845);
ii=ii+1; mob_eq_comp{ii} = 'Fe(3)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/((1e3)*55.845);
ii=ii+1; mob_eq_comp{ii} = 'K'; % *****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 5.885; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 3.276; % mg/L 
convert_obs2model(ii) = 1/((1e3)*39.098);
ii=ii+1; mob_eq_comp{ii} = 'Mg'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 56.2; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 43.18; % mg/L 
convert_obs2model(ii) = 1/((1e3)*24.305);
ii=ii+1; mob_eq_comp{ii} = 'Mn(2)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.08415; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/((1e3)*54.938);
ii=ii+1; mob_eq_comp{ii} = 'Amm'; % ***************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/1e6; % Amm is NH3
ii=ii+1; mob_eq_comp{ii} = 'Na'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 44.9; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 24.494; % mg/L 
convert_obs2model(ii) = 1/((1e3)*22.990);
ii=ii+1; mob_eq_comp{ii} = 'O(0)'; % **************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.01; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.6; % mg/L 
convert_obs2model(ii) = 1/((1e3)*16.0);
ii=ii+1; mob_eq_comp{ii} = 'P'; % *****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.011; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.004; % mg/L 
convert_obs2model(ii) = 1/((1e3)*30.9737);
ii=ii+1; mob_eq_comp{ii} = 'S(6)'; % **************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 239; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 239; % mg/L 
convert_obs2model(ii) = 1/((1e3)*96.066);
ii=ii+1; mob_eq_comp{ii} = 'S(-2)'; % *************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/((1e3)*33.073);
ii=ii+1; mob_eq_comp{ii} = 'Si'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 14.81; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 14.052; % mg/L 
convert_obs2model(ii) = 1/((1e3)*28.086);
ii=ii+1; mob_eq_comp{ii} = 'N(5)'; % **************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.7138; % mg/L 
convert_obs2model(ii) = 1/((1e3)*62.0049);
ii=ii+1; mob_eq_comp{ii} = 'N(3)'; % **************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0.0058; % mg/L 
convert_obs2model(ii) = 1/((1e3)*46.0049);
ii=ii+1; mob_eq_comp{ii} = 'N(0)'; % **************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 0; % mg/L 
convert_obs2model(ii) = 1/((1e3)*28.012);
ii=ii+1; mob_eq_comp{ii} = 'pH'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 7.4; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 7.2; % mg/L 
convert_obs2model(ii) = 1;
ii=ii+1; mob_eq_comp{ii} = 'pe'; % ****************************************
zz = 0;
zz = zz + 1; mob_eq_obs_z(zz,ii) = 12; % mg/L 
zz = zz + 1; mob_eq_obs_z(zz,ii) = 12; % mg/L 
convert_obs2model(ii) = 1;
%**************************************************************************

n_mob_eq = ii;
mob_eq_comp = mob_eq_comp(1:n_mob_eq);
mob_eq_obs_z = mob_eq_obs_z(:,1:n_mob_eq);
mob_eq_extra_z = mob_eq_extra_z(:,1:n_mob_eq);
convert_obs2model = convert_obs2model(1:n_mob_eq);
mob_eq_obs_z = mob_eq_obs_z .* repmat(convert_obs2model, nzones, 1);
%**************************************************************************

% - mineral eq components
% (only allowed for porewater solution)
n_min_eq_max = 10;
min_eq_comp = cell(n_min_eq_max,1);
min_eq_ic_z = zeros(nzones,n_min_eq_max); 
ii = 0;
ii=ii+1; min_eq_comp{ii} = 'Orpiment'; % ***********************************
zz = 0;
zz = zz + 1; min_eq_ic_z(zz,ii) = 1.29e-5;  % based on XANES in Nicholas et al 2017
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % arbitrary conc
ii=ii+1; min_eq_comp{ii} = 'Arsenopyrite'; % ***********************************
zz = 0;
zz = zz + 1; min_eq_ic_z(zz,ii) = 7.79e-5;  % based on XANES in Nicholas et al 2017
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % arbitrary conc
ii=ii+1; min_eq_comp{ii} = 'Fe(OH)3(a)'; % ***********************************
zz = 0;
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % based on XANES in Nicholas et al 2017
zz = zz + 1; min_eq_ic_z(zz,ii) = 0;  % arbitrary conc
% ii=ii+1; min_eq_comp{ii} = 'Pyrite'; % ***********************************
% min_eq_pw_ic(:,:,:,ii) = 0; % (mol/Lv)
 %ii=ii+1; min_eq_comp{ii} = 'Goethite'; % ***********************************
 %min_eq_ic(:,:,:,ii) = 50e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Siderite'; % ***********************************
% min_eq_ic(:,:,:,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Rhodochrosite'; % ***********************************
% min_eq_ic(:,:,:,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'FeS(ppt)'; % ***********************************
% min_eq_ic(:,:,1:50,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Mackinawite'; % ***********************************
% min_eq_ic(:,:,1:50,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Gypsum'; % ***********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Hydroxyapatite'; % ***********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Hematite'; % ***********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*1864; % (mol/Lv)
% ii=ii+1; min_eq_comp{ii} = 'Vivianite'; % ***********************************
% min_eq_ic(:,:,51:100,ii) = 1e-6*1864; % (mol/Lv)

n_min_eq = ii;
min_eq_comp = min_eq_comp(1:n_min_eq);
min_eq_ic_z = min_eq_ic_z(:,1:n_min_eq);  

% - catex components (catex_exch is exchanger site)
n_catex_exch_max = 10;
catex_exch_comp = cell(n_catex_exch_max,1);
catex_exch_pw_ic = zeros(1,n_catex_exch_max); 
ii = 0;
if fl_Yexch
    ii=ii+1; catex_exch_comp{ii} = 'Y'; % ***********************************
    catex_exch_pw_ic(ii) = 0.020725*por;  % Doug's xls "Sorption Notes" mol/L_w, email 4/12/13
end
n_catex_exch = ii;
catex_exch_comp = catex_exch_comp(1:n_catex_exch);
catex_exch_pw_ic = catex_exch_pw_ic(1:n_catex_exch); 

% - catex components (catex is exchanger sorbed compound)
% (ic from PHREEQC batch run)
n_catex_max = 10;
catex_comp = cell(n_catex_max,1);
catex_pw_ic = zeros(1,n_catex_max); 
ii = 0;
if fl_Yexch
    ii=ii+1; catex_comp{ii} = 'HY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'NaY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'KY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'AmmHY'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'CaY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'MgY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'FeY2'; % ***********************************
    ii=ii+1; catex_comp{ii} = 'MnY2'; % ***********************************
%     ii=ii+1; catex_comp{ii} = 'SrY2'; % ***********************************
%     ii=ii+1; catex_comp{ii} = 'BaY2'; % ***********************************
end
n_catex = ii;
catex_comp = catex_comp(1:n_catex);
catex_pw_ic = catex_pw_ic(1:n_catex);



% - surface complexation sites
% **** WARNING: not set up for surfaces w/
%               multiple sites (would need to enter only 1 surf_calc_type
%               per surface in _ph file)
n_surf_max = 10;
surf_comp = cell(n_surf_max,1);
% surf_ic: number of sites [mol/Lv] if not coupled to phase or reactant,
%          (number of sites per mol x porosity) [mol/mol] if coupled
surf_ic = zeros(1,n_surf_max); 
% surf_par(1,ii): SurfArea ([m2/g] if not coupled to phase or reactant,
%              [m2/mol] if coupled to phase or reactant)
% surf_par(2,ii): mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_par = zeros(2,n_surf_max); 
% surf_cpl{1,ii}: Empty if not coupled to phase or reactant
% surf_cpl{1,ii}: Name of pure phase or kin reactant coupled to
% surf_cpl{2,ii}: 'equilibrium_phase' to couple to pure phase, or
%                    'kinetic_reactant' to couple to kinetic reactant
surf_cpl = cell(2,n_surf_max); 
surf_calc_type = ''; % 3 options: '', '-no_edl', 'diffuse_layer'
ii = 0;
ii=ii+1; surf_comp{ii} = 'Hfo_wOH'; % ***********************************
surf_ic(ii) = 0.2*por;  % number of sites, [mol/mol]*por if coupled
surf_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_cpl{1,ii} = 'Fe(OH)3(a)';
surf_cpl{2,ii} = 'equilibrium_phase';
ii=ii+1; surf_comp{ii} = 'Hfo_sOH'; % ***********************************
surf_ic(ii) = 0.005*por;  % number of sites, [mol/mol]*por if coupled
surf_par(1,ii) = 1243.5/0.00012;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_cpl{1,ii} = 'Fe(OH)3(a)';
surf_cpl{2,ii} = 'equilibrium_phase';  
n_surf = ii;
surf_comp = surf_comp(1:n_surf);
surf_ic = surf_ic(:,1:n_surf); 
surf_par = surf_par(:,1:n_surf); 

%( Take out this section, unless you want to see equilibrated surf)
% - surface complexation specie components (surf_comp is for sites only)
% **** WARNING: not set up for surfaces w/
%               multiple sites (would need to enter only 1 surf_calc_type
%               per surface in _ph file)
n_surf_sp_max = 10;
surf_sp_comp = cell(n_surf_sp_max,1);
% surf_sp_ic: number of sites [mol/Lv] if not coupled to phase or reactant,
%          (number of sites per mol x porosity) [mol/mol] if coupled
surf_sp_pw_ic = zeros(1,n_surf_sp_max); 
% surf_sp_par(1,ii): SurfArea ([m2/g] if not coupled to phase or reactant,
%              [m2/mol] if coupled to phase or reactant)
% surf_sp_par(2,ii): mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_sp_par = zeros(2,n_surf_sp_max); 
% surf_sp_cpl{1,ii}: Empty if not coupled to phase or reactant
% surf_sp_cpl{1,ii}: Name of pure phase or kin reactant coupled to
% surf_sp_cpl{2,ii}: 'equilibrium_phase' to couple to pure phase, or
%                    'kinetic_reactant' to couple to kinetic reactant
surf_sp_cpl = cell(2,n_surf_sp_max); 
% surf_calc_type = ''; % 3 options: '', '-no_edl', 'diffuse_layer'
ii = 0;
ii=ii+1; surf_sp_comp{ii} = 'Hfo_wOH'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';
ii=ii+1; surf_sp_comp{ii} = 'Hfo_wH2AsO4'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';
ii=ii+1; surf_sp_comp{ii} = 'Hfo_wHAsO4-'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';
ii=ii+1; surf_sp_comp{ii} = 'Hfo_wOHAsO4-3'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
%              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';
% ii=ii+1; surf_sp_comp{ii} = 'Hfo_wOH2+'; % ***********************************
% surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
% surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
% %              include but ignored if coupled to phase or reactant)
% surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
% surf_sp_cpl{2,ii} = 'equilibrium_phase';
% ii=ii+1; surf_sp_comp{ii} = 'Hfo_wOCa+'; % ***********************************
% surf_sp_par(1,ii) = 1243.5/0.004665;  % SurfArea ([m2/mol] if coupled to phase or reactant)
% surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
% %              include but ignored if coupled to phase or reactant)
% surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
% surf_sp_cpl{2,ii} = 'equilibrium_phase'; 
ii=ii+1; surf_sp_comp{ii} = 'Hfo_sOH'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.00012;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
% %              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';  
ii=ii+1; surf_sp_comp{ii} = 'Hfo_sOFe+'; % ***********************************
surf_sp_par(1,ii) = 1243.5/0.00012;  % SurfArea ([m2/mol] if coupled to phase or reactant)
surf_sp_par(2,ii) = 4145;  % mass ([g/L] if not coupled to phase or reactant,
% %              include but ignored if coupled to phase or reactant)
surf_sp_cpl{1,ii} = 'Fe(OH)3(a)';
surf_sp_cpl{2,ii} = 'equilibrium_phase';

n_surf_sp = ii;
surf_sp_comp = surf_sp_comp(1:n_surf_sp);
surf_sp_pw_ic = surf_sp_pw_ic(:,1:n_surf_sp); 
surf_sp_par = surf_sp_par(:,1:n_surf_sp); 


% -- (done setting model chemical components)

% - select output list
phrq_sel_outlist.total = [mob_eq_comp(~strncmp('p',mob_eq_comp,1))];
phrq_sel_outlist.mol = [catex_comp; surf_sp_comp];
phrq_sel_outlist.equilphase = min_eq_comp;
phrq_sel_outlist.si = []; phrq_sel_outlist.gas = [];    

sel_outlist.total = [mob_eq_comp(~strncmp('p',mob_eq_comp,1))];
sel_outlist.mol = catex_comp;  % better to put desired surf species in add_sel_outlist.mol
% sel_outlist.mol = [catex_comp; surf_comp];
sel_outlist.equilphase = min_eq_comp;
% sel_outlist.si = {'Siderite'}; 
sel_outlist.si = {'Fe(OH)3(a)'}; 
% sel_outlist.si = []; 
sel_outlist.gas = []; 
% sel_outlist.si = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 
% sel_outlist.gas = {'CH4(g)', 'CO2(g)', 'Ntwo(g)', 'O2(g)'}; 

% iunclude any additional components
% sel_outlist.total = [sel_outlist.total; addl_sel_outlist.total(:)];
% sel_outlist.mol = [sel_outlist.mol; addl_sel_outlist.mol(:)];
% sel_outlist.equilphase = [sel_outlist.equilphase; addl_sel_outlist.equilphase(:)];
% sel_outlist.si = [sel_outlist.si; addl_sel_outlist.si(:)]; 
% sel_outlist.gas = [sel_outlist.gas; addl_sel_outlist.gas(:)]; 



% -- Get initial conditions w/ PHREEQC

%     % (re-)run PHREEQC to equilibrate ic conc for each zone
%     fprintf('Have not finished code of fl_init_phrq=1!  Need to do for each zone!  Exiting...\n');
%     return

mob_eq.name = mob_eq_comp; 


catex = []; % sw: no cation exchange
surf = []; % sw: no surface complexation
fl_force_redox = 1;  % force redox equilibration, not just charge balance

% catex_ic_equil = nan(1,n_catex);
mob_eq_ic_equil = zeros(2, n_mob_eq);  % sw, pw
for ii = 1:2  % 1: sw, 2: pw
    if ii == 1 % sw
        fprintf('\n--- Surface water ---');
        mob_eq.ic = mob_eq_sw_ic;
        phase.name = min_eq_comp; phase.si = zeros(n_min_eq,1); 
        phase.ic = zeros(1,n_min_eq);
        phase.ic(:,strcmp(phase.name, 'Arsenopyrite')) = 6.23e-5; phase.ic(:,strcmp(phase.name, 'Orpiment')) = 1.03e-5; 
%         tempC = tempC_sw;
        mob_eq.extra = mob_eq_extra_z(1,:);
        if n_catex > 0 % cation exchange components
            % assume exchanger is same for all zones
            catex.name = catex_exch_comp; 
            catex.cec = catex_exch_pw_ic / por;
        end
        if n_surf_sp > 0 && n_surf > 0 % surface complexation components, see equilibrated species results
            surf.name = surf_comp; 
            surf.num_sites = surf_ic / por;
            surf.spec_area = surf_par(1,:);
            surf.mass = surf_par(2,:);
            surf.cpl_phase = surf_cpl(1,:);
            surf.cpl_opt = surf_cpl(2,:);
            surf.type = surf_calc_type;
        end
    elseif ii == 2 % pw
        fprintf('\n--- Pore water ---');
        mob_eq.ic = mob_eq_pw_ic;
        phase.name = min_eq_comp; phase.si = zeros(n_min_eq,1); 
        phase.ic = zeros(1,n_min_eq);
        phase.ic(:,strcmp(phase.name, 'Fe(OH)3(a)')) = 0; %phase.ic(:,strcmp(phase.name, 'Pyrite')) = 0;
%         phase.name = {}; phase.si = [];
%         phase.ic = [];
%         tempC = tempC_pw;
        mob_eq.extra = mob_eq_extra_z(2,:);
        if n_catex > 0 % cation exchange components
            % assume exchanger is same for all zones
            catex.name = catex_exch_comp; 
            catex.cec = catex_exch_pw_ic / por;
        end
        if n_surf_sp > 0 && n_surf > 0 % surface complexation components, see equilibrated species results
            surf.name = surf_comp; 
            surf.num_sites = surf_ic / por;
            surf.spec_area = surf_par(1,:);
            surf.mass = surf_par(2,:);
            surf.cpl_phase = surf_cpl(1,:);
            surf.cpl_opt = surf_cpl(2,:);
            surf.type = surf_calc_type;
        end
    end    
    
    [mob_eq_ic_equil(ii,:), b, c] = generate_ic_PHREEQC_f_160606a(...
        phrq_exe, phrq_infil, phrq_outfil, file_databas, phrq_sel_outfil, phrq_sel_outlist, ...
        por, tempC, units, ...
        mob_eq, phase, catex, surf, ...
        mob_eq_comp, catex_comp, surf_sp_comp, fl_force_redox, ...
        phrq_tblfil);
    if ~isempty(b)
%         catex_sw_ic = b;
        catex_pw_ic = b;
    end    
    if ~isempty(c)
        surf_sp_pw_ic = c;
    end    
  
    movefile(phrq_infil, [phrq_infil, num2str(ii)]);
    movefile(phrq_outfil, [phrq_outfil, num2str(ii)]);
    movefile(phrq_sel_outfil, [phrq_sel_outfil, num2str(ii)]);
    movefile(phrq_tblfil, [phrq_tblfil, num2str(ii)]);
    
end

mob_eq_sw_ic = mob_eq_ic_equil(1,:);
mob_eq_pw_ic = mob_eq_ic_equil(2,:);


% rmpath(SC_dir);
