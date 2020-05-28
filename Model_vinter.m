% 'g' = den forede gravitationsledning, 'tryk' = trykledningen, 'PVC' = PVC-gravitationsledningen
%% IMPORT AF OPHOLDSTIDER OG VANDF?RING
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "opholdstider_samlet";
opts.DataRange = "A2:D8641";

% Specify column names and types
opts.VariableNames = ["tid", "LilleQ", "MellemQ", "StoreQ", "sum_lille", "sum_mellem", "sum_stor", "Relativ_flow", "Relativ_flow_tryk"];
opts.SelectedVariableNames = ["tid", "LilleQ", "MellemQ", "StoreQ", "sum_lille", "sum_mellem", "sum_stor", "Relativ_flow", "Relativ_flow_tryk"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
OpholdstiderS3 = readtable("Opholdstider.xlsx", opts, "UseExcel", false);
% Convert to output type
OpholdstiderS3 = table2array(OpholdstiderS3);

% Clear temporary variables
clear opts

% Tidspunkter hvor pumpen, der leder vand til trykledningen, k?rer
pumpeaktiv_lille = OpholdstiderS3(:,5);
pumpeaktiv_mellem = OpholdstiderS3(:,6);
pumpeaktiv_stor = OpholdstiderS3(:,7);

% Opholdstider til n?r pumpen, der leder vand til trykledningen, k?rer
opholdstid_lille = OpholdstiderS3(:,2)*(60*60*24).*pumpeaktiv_lille; % [s]
opholdstid_mellem = OpholdstiderS3(:,3)*(60*60*24).*pumpeaktiv_mellem; % [s]
opholdstid_stor = OpholdstiderS3(:,4)*(60*60*24).*pumpeaktiv_stor; % [s]

% Vandf?ring i l?bet af dagen i den forede ledning
variation_byspildevand_g = OpholdstiderS3(:,8);
flow_g_mellem = (6281.5.*variation_byspildevand_g)/10; % vinter [m3/s]
%flow_g_mellem_best = (14165.5.*variation_byspildevand_g)/10; % vinter [m3/s]
%flow_g_mellem_worst = (3526.5.*variation_byspildevand_g)/10; % vinter [m3/s]

% Hvor meget vand der udledes n?r pumpen k?rer
% Pumpekapaciteten er 12 L/s
variation_tryk = OpholdstiderS3(:,9);
flow_tryk_lille =  12*0.001.*pumpeaktiv_lille; % [m3/s] 
flow_tryk_mellem =  12*0.001.*pumpeaktiv_mellem; %[m3/s]
flow_tryk_stor = 12*0.001.*pumpeaktiv_stor; %[m3/s] 
variation_flow_tryk = (85.*variation_tryk)/10; %[m3/s]
%% TRYKLEDNINGEN
% Geometriske beregninger for trykledningen
Areal_tryk = pi*(49.94+74.16+1497.73+1591.16+506.07)*0.151 + pi*41.44*0.148; %[m2]
Volume_tryk = ((0.151/2)^2)*pi*(49.94+74.16+1497.73+1591.16+506.07) + ((0.148/2)^2)*pi*41.44; %[m3]
ArealVolume_tryk = Areal_tryk/Volume_tryk; %[1/m]

% V?rdier for sulfiddannelsesraten i trykledninger [g S/m2*h]
%Rate_industri = 0.1/(60*60); % industrispildevand [gS/m2*s]

%Udregning af raten for sulfiddannelse. 'a' for industri med meget org.
%stof er mellem 0,007-0,010 g S gO2 / m*h(s. 234 i SP)
%COD antages 700 g O2/m3
temperatur = 8; % grader
Rate_a_industri = (0.008/(60*60))*((950-50)^0.5)*1.03^(temperatur - 20);

% Koncentrationen af sulfid der dannes i trykledningen beregnes ud fra den
% ovenfor n?vnte rate, areal/volumenforholdet og opholdstider [gS/m3]

C_tryk=zeros(length(opholdstid_mellem),1);
for i=1:length(opholdstid_mellem)

     C_tryk(i) = Rate_a_industri*ArealVolume_tryk.*opholdstid_mellem(i); 

     if C_tryk(i) > 28.3
     C_tryk(i) = 28.3;
     else
     C_tryk(i) = Rate_a_industri*ArealVolume_tryk.*opholdstid_mellem(i);      
     end
end
%% PVC-GRAVITATIONSLEDNINGEN - VOLUMEN AF VAND- OG GASFASE 
% Fysiske parametre for PVC-gravitationsledning
Laengde_PVC = 1523.6; % [m]
Haeldning_PVC = 0.00119; % [m/m]
Diameter_PVC = 0.291; % [m] 0,315 som ydre diameter
Ruhed_PVC = 0.00003; %[m] Denne er antaget fra "normal" PVC r??r = 0,01-0,05 mm 
Manningtal_PVC = 25.4/(Ruhed_PVC)^(1/6); %[-]

% Vandf?ring n?r PVC-gravitationsledningen er fuldtl?bende
Hydraulisk_radius_fuld_PVC = Diameter_PVC/4; %[m]
v_fuld_PVC = Manningtal_PVC*(Hydraulisk_radius_fuld_PVC)^(2/3)*Haeldning_PVC^0.5; %[m/s] 
Areal_fuld_PVC = pi*((Diameter_PVC/2)^2); %[m2] radius for fuldl?bende r?r = diameter/4
Q_fuld_PVC = v_fuld_PVC*Areal_fuld_PVC; %[m3/s]

% Vanddybde afh?ngig af vandf?ring (Q), n?r PVC-gravitationsledningen er
% delfyldt
vanddybde_PVC=zeros(length(flow_tryk_mellem),1);
 for i=1:length(flow_tryk_mellem)
     
     vanddybde_PVC(i) = (0.3183098862*(3.141592654-1.*acos((0.125*(-25.*Q_fuld_PVC+sqrt(800.*Q_fuld_PVC*flow_tryk_mellem(i)+289.*Q_fuld_PVC^2)))/Q_fuld_PVC)))*Diameter_PVC;
     
     if flow_tryk_mellem(i) >0 
     else
     vanddybde_PVC(i) = 0;
     end
 end 
% y = 0.3183098862*(3.141592654 - 1.*arccos(0.1250000000*(-25.*Q__f + sqrt(800.*Q__f*Q + 289.*Q__f^2))/Q__f))*d

% Nedenst?ende er beregninger der er n?dvendige for at kunne beregne
% geniltningskoefficienten (KLa) for PVC-gravitationsledningen

% Teta_PVC er vinklen mellem centrum af r?ret til hvor vandspejlet rammer
% PVC-gravitationsledningens overflade
teta_PVC = acos(1-((2.*vanddybde_PVC)./Diameter_PVC));% [rad]
% Tv?rsnitsarealet af vandfasen
Areal_vand_PVC = ((Diameter_PVC^2)/8)*(2.*teta_PVC-sin(2.*teta_PVC));% [m2]

% Hastigheden af vandet i delfyldt r?r
v_delfuld_PVC = flow_tryk_mellem./Areal_vand_PVC; % [m/s]
v_delfuld_PVC(isnan(v_delfuld_PVC))=0; 

% Middelvanddybden i delfyldt r?r
bredde_vandspejl_PVC = sqrt(1-(1-((2*vanddybde_PVC)./Diameter_PVC)).^2).*Diameter_PVC; % [m]
middel_vanddybde_PVC = Areal_vand_PVC./bredde_vandspejl_PVC; % [m]
middel_vanddybde_PVC(isnan(middel_vanddybde_PVC))=0;

% Froudes tal bestemmer om str?mningen er strygende eller str?mmende, 
% under 1: str?mmende
% over 1: strygende
Froude_PVC = v_delfuld_PVC.*(9.816.*middel_vanddybde_PVC).^(-0.5); % [-]

% Geniltningskoefficient g?ldende for PVC-gravitationsledningen 
KLa_ilt_PVC = (0.86.*(1+0.2.*(Froude_PVC).^2).*(Haeldning_PVC.*v_delfuld_PVC).^(3/8).*(1./middel_vanddybde_PVC)*1.024^(temperatur-20))/3600; %[1/s] 
KLa_ilt_PVC(isnan(KLa_ilt_PVC))=0;

% Overfladeareal af betonr?rets overflade i kontakt med gasfasen
areal_gasfase_delfuld_PVC = (Diameter_PVC*pi - teta_PVC.*Diameter_PVC).*5.8; %[m2]
% Volumen af vandfasen i delfyldt r?r
volume_vandfase_delfuld_PVC = ((Diameter_PVC^2)/8.*(2.*teta_PVC - sin(2.*teta_PVC))).*5.8; %[m3]
% Volumen af vandfasen i fyldt r?r
volume_total_fuld_PVC = ((pi*(Diameter_PVC/2)^2).*5.8); %[m3]
% Volumen af gasfasen i delfyldt r?r
volume_gasfase_delfuld_PVC = volume_total_fuld_PVC - volume_vandfase_delfuld_PVC; %[m3]

% Der laves matricer
areal_gasfase_delfuld_PVC_matrix = areal_gasfase_delfuld_PVC.*ones(8640,263);
volume_gasfase_delfuld_PVC_matrix = volume_gasfase_delfuld_PVC.*ones(8640,263);
%% PVC-GRAVITATIONSLEDNING - SULFID OXIDATION

% Maksimale m?ngde ilt, som kan optages i vandet i
% PVC-gravitationsledningen
ilt_max_optag_PVC = 14.652 - 0.41022*temperatur+0.00799*temperatur^2-0.0000777*temperatur^3; %mg/L eller g/m3; og T i grader

% Geniltningsrate i PVC-gravitationsledningen
geniltningsrate = KLa_ilt_PVC*(ilt_max_optag_PVC-0); %mg/L*s eller g/m3*s - 0 er iltkoncentration til start
total_ilt_optag_i_PVC = (Laengde_PVC./v_delfuld_PVC).*geniltningsrate;% g/m3
total_ilt_optag_i_PVC(isnan(total_ilt_optag_i_PVC))=0;

% Hvor meget sulfid der oxideres i PVC-gravitationsledningen:
% Der optages 5,347 g ilt pr. 1 m3. 
% Reaktion for oxidation af sulfid til sulfat: 2HS- + 4O2 -> 2SO4^2- + 2H+
% Der forbruges 4 iltmolekyler til 2 sulfidmolekyler. Forhold: 4:2 = 2:1.
% Det antages at 50% af den optagede ilt i vandet forbruges til sulfidoxidation
ilt_forbrug_til_sulfid_oxidation = 0.5*total_ilt_optag_i_PVC; %gO2/m3 
Mw_oxygen = 15.999; % [g/mol]
Mw_svovl = 32.065; % [g/mol]
forbrug_sulfid_pr_ilt = Mw_svovl/(2*(Mw_oxygen*2)); % g sulfid som forbruges pr. 1 g ilt 
forbrug_sulfid_pr_m3 = forbrug_sulfid_pr_ilt * ilt_forbrug_til_sulfid_oxidation; %g S/m3

% M?ngde sulfid der kommer ud af PVC-gravitationsledningen: 
% m?ngde sulfid der oxideres, tr?kkes fra den totale m?ngde sulfid, 
% der blev dannet i trykledningen
C_PVC=zeros(length(C_tryk),1);
for i=1:length(C_tryk)
    if C_PVC(i) < 0
    C_PVC(i) = 0;
    else 
    C_PVC(i) = C_tryk(i) - forbrug_sulfid_pr_m3(i);
    end
end
%% PVC-GRAVITATIONSLEDNING - EMISSION FRA VAND- TIL GASFASE
% V?rdien for pH er antaget
pH = 6.9; % [-]

% Svovlbrinteudvekslingskoefficient beregnes p? baggrund af
% iltudvekslingskoefficienten og en omregningsfaktor
KLa_h2s_PVC = max(KLa_ilt_PVC)*0.84; % [1/s]
KLa_h2s_matrix_PVC= KLa_h2s_PVC.*ones(8640,263); %[1/s]

% Fraktion af H2S (aq) i vandfasen; pKa = 7.1
H2S_aq_PVC = C_PVC./(exp(-log(10)*7.1 + log(10)*pH)+1); %[g H2S/m3]
H2S_vandfase_PVC = [H2S_aq_PVC.*ones(8640,2) zeros(8640,261)]; 

% Konversionsfaktor af ppm til g/m3 og omvendt
volume_idealgas_temperatur = (1*8.314472*(273.15+temperatur))/101325;
konversionsfaktor = (1/volume_idealgas_temperatur)*32.065*10^-6;

% Henrys konstant korrigeret efter temperatur
henrys_konstant = 1./ (0.1.*exp( 2100*(1./(temperatur+273.15) - 1./298.15) ) ) .*ones(8640,263); %[atm]

% Frigivelse af H2S fra vand- til gasfase; der ganges med tidsskridtet p?

H2S_frigivelse_PVC = [((KLa_h2s_PVC.*(H2S_aq_PVC - 0).*1.024^(temperatur-20)).*10).*ones(8640,2) zeros(8640,261)]; %Antages at for nedre rand er gasfasen nul ppm
H2S_eq_PVC = zeros(8640,263);
H2S_gasfase_PVC = zeros(8640,263);
H2S_gasfase_ppm_PVC = zeros(8640,263);

H2S_oxidation_PVC_gm2s = zeros(8640,263);
H2S_oxidation_PVC = zeros(8640,263);

for i=2:8640
    for j=2:263 % F?rste kolonne er en fiktiv nedre rand
        H2S_vandfase_PVC(i,j) = H2S_vandfase_PVC(i-1,j-1) - H2S_frigivelse_PVC(i-1,j-1); %[g S/m3]
                
        H2S_frigivelse_PVC(i,j) = KLa_h2s_matrix_PVC(i,j).*(H2S_vandfase_PVC(i,j)- H2S_eq_PVC(i,j) ).*1.024^(temperatur-20).*10; %[g S/m3]
            
        H2S_gasfase_PVC(i,j) = H2S_frigivelse_PVC(i-1,j) + H2S_gasfase_PVC(i-1,j-1) - H2S_oxidation_PVC(i,j); %[g S/m3]

        H2S_gasfase_ppm_PVC(i,j) = H2S_gasfase_PVC(i-1,j)./konversionsfaktor; %[ppm]
        
        H2S_eq_PVC(i,j) = 32.065.*10.^-3.*(H2S_gasfase_ppm_PVC(i,j)./henrys_konstant(i,j)) ; %[g S/m3] 
        
        H2S_oxidation_PVC_gm2s(i,j) = ( 0.4.*H2S_gasfase_ppm_PVC(i,j).^0.62) .*0.001*(1/3600); %[g S/m2*s] 0.4 = K_F og n=0.62
        
        H2S_oxidation_PVC(i,j) = H2S_oxidation_PVC_gm2s(i,j).*(areal_gasfase_delfuld_PVC_matrix(i,j)./volume_gasfase_delfuld_PVC_matrix(i,j).*10); %[g S/m3]
 
        
    end
end
%% DEN FOREDE LEDNING - VOLUMEN AF VAND- OG GASFASE
% Fysiske parametre
Laendge_g = 61.77; % [m]
Haeldning_g = 0.00258; % [m/m]
Diameter_g = 0.90-(2*0.01); % [m] Medregnet tykkelse af m?rtel
Ruhed_g = 0.003; %[m] Denne er antaget fra "normal" betonr?r
Manningtal_g = 25.4/(Ruhed_g)^(1/6); %[-]

% Vandf?ring i den forede ledning, n?r den er fuldtl?bende
Hydraulisk_radius_fuld_g = Diameter_g/4; %[m]
v_fuld_g = Manningtal_g*(Hydraulisk_radius_fuld_g)^(2/3)*Haeldning_g^0.5; %[m/s] 
Areal_fuld_g = pi*(Diameter_g/2)^2; %[m2]
Q_full_g = v_fuld_g*Areal_fuld_g; %[m3/s]

% Vanddybde afh?ngig af vandf?ring (Q)
vanddybde_g = (.3183098862*(3.141592654-1.*acos((.1250000000*(-25.*Q_full_g+sqrt(800.*Q_full_g*flow_g_mellem+289.*Q_full_g^2)))/Q_full_g)))*Diameter_g;

% Nedenst?ende beregninger leder hen til beregningerne af volumen af vand-
% og gasfase

% Teta_g er vinklen mellem centrum af r?ret til hvor vandspejlet rammer
% betonr?rets overflade
teta_g = acos(1-((2.*vanddybde_g)/Diameter_g));% [rad]

% Tv?rsnitsarealet af vandfasen
areal_vand_g = Diameter_g^2/8*(2*teta_g-sin(2*teta_g));%[m2]

% Hastighed af vandet i delfyldt r?r
v_delfuld_g = flow_g_mellem./areal_vand_g; %[m/s]

% Middelvanddybde i delfyldt r?r
bredde_vandspejl_g = sqrt(1-(1-2*vanddybde_g/Diameter_g).^2)*Diameter_g; % [m]
middel_vanddybde_g = areal_vand_g./bredde_vandspejl_g; %[m]

% Froudes tal bestemmer om str?mningen er strygende eller str?mmende, 
% under 1: str?mmende
% over 1: strygende
Froude_g = v_delfuld_g.*(9.182.*middel_vanddybde_g).^-0.5; %[-]

% Reduceret areal af betonr?rets overflade i kontakt med gasfasen;
% s?ledes der tages h?jde for at der kun sker oxidation af svovlbrinte til
% svovlsyre i de fugtige omr?der - IKKE MEDTAGET ENDNU
% biofilm_av = (teta.*Diameter_gravitationsledning)./crossarea_vand; %[1/m]

% Overfladeareal af betonr?rets overflade i kontakt med gasfasen
areal_gasfase_delfuld_g = (Diameter_g*pi - teta_g.*Diameter_g).*5; %[m2]
% Volumen af vandfasen i delfyldt r?r
volume_vandfase_delfuld_g = ((Diameter_g^2)/8.*(2.*teta_g - sin(2.*teta_g))).*5; %[m3]
% Volumen af vandfasen i fyldt r?r
volume_total_fuld_g = ((pi*(Diameter_g/2)^2).*5); %[m3]
% Volumen af gasfasen i delfyldt r?r
volume_gasfase_delfuld_g = volume_total_fuld_g - volume_vandfase_delfuld_g; %[m3]

% Der laves matricer
areal_gasfase_delfuld_g_matrix = areal_gasfase_delfuld_g.*ones(8640,11);
volume_gasfase_delfuld_g_matrix = volume_gasfase_delfuld_g.*ones(8640,11);
%% DEN FOREDE LEDNING - KONCENTRATION AF SULFID
% Koncentration af sulfid i husholdningsspildevand; antaget!
C_byspildevand = 0.1; %[gS/m3]

tid_forskydning_matrix = zeros(round(1522/max(v_delfuld_PVC)/10) ,1);
C_PVC_tid_forskydning_bund =  C_PVC((length(C_PVC) - length(tid_forskydning_matrix)):end);
C_PVC_tid_forskydning_top = C_PVC(1: (length(C_PVC) - length(tid_forskydning_matrix)-1 ));

flow_PVC_tid_forskydning_bund = flow_tryk_mellem((length(flow_tryk_mellem) - length(tid_forskydning_matrix)):end);
flow_PVC_tid_forskydning_top = flow_tryk_mellem(1: (length(flow_tryk_mellem) - length(tid_forskydning_matrix)-1 ));

C_PVC_tid_korregeret = [C_PVC_tid_forskydning_bund; C_PVC_tid_forskydning_top];
flow_PVC_tid_korregeret = [flow_PVC_tid_forskydning_bund; flow_PVC_tid_forskydning_top];

%Samlet sulfid koncentration efter emission af svovlbrinte til gasfasen
C_PVC_ud = H2S_vandfase_PVC(:,263).*10^(pH-7.1)+H2S_vandfase_PVC(:,263);

% Koncentration af sulfid i spildevandet i den forede ledning
C_total = ((C_byspildevand.*(flow_g_mellem-flow_PVC_tid_korregeret))+(C_PVC_ud.*flow_PVC_tid_korregeret))./(flow_g_mellem);
%% DEN FOREDE LEDNING - KONCENTRATIONER I BR?NDEN 
volume_broend = pi*0.4^(2)*3.01; %0,4 er br?ndens radius
gas_i_broend = volume_broend - flow_g_mellem.*(0.8./v_delfuld_g); %0,8 er diamteren af br?nden

masse_H2S_PVC = H2S_gasfase_PVC(:,263).*volume_gasfase_delfuld_PVC;
koncentration_H2S_broend = masse_H2S_PVC./gas_i_broend;
%% DEN FOREDE LEDNING - EMISSION FRA VAND- TIL GASFASE
% V?rdien for pH er antaget
pH = 7.5; % [-]

% Svovlbrinteudvekslingskoefficient beregnes p? baggrund af
% iltudvekslingskoefficienten og en omregningsfaktor
KLa_ilt_g = 0.86.*(1+0.2.*(Froude_g).^2).*(Haeldning_g.*v_delfuld_g).^(3/8).*(1./middel_vanddybde_g)*(1/3600); %[1/s]
KLa_h2s_g = KLa_ilt_g*0.84; % [1/s]
KLa_h2s_matrix_g = [KLa_h2s_g.*2.*ones(8640,2) KLa_h2s_g.*ones(8640,9)]; %[1/s]

% Fraktion af H2S (aq) i vandfasen; pKa = 7.1
H2S_aq = C_total./(exp(-log(10)*7.1 + log(10)*pH)+1); %[g H2S/m3]
H2S_vandfase = [H2S_aq.*ones(8640,2) zeros(8640,9)]; 

% Konversionsfaktor af ppm til g/m3 og omvendt
volume_idealgas_temperatur = (1*8.314472*(273.15+temperatur))/101325;
konversionsfaktor = (1/volume_idealgas_temperatur)*32.065*10^-6;

% Henrys konstant korrigeret efter temperatur
henrys_konstant = 1./ (0.1.*exp( 2100*(1./(temperatur+273.15) - 1./298.15) ) ) .*ones(8640,11); %[atm]

% Kinetiske parametre
K_F = 15.52*ones(8640,11);
n = 0.55*ones(8640,11);

% Frigivelse af H2S fra vand- til gasfase; der ganges med tidsskridtet p?
% 10 sekunder
% Gennemsnitshastigheder der afg?r antal stedsskridt (kolonner), j: 
% lille vandf?ring forede ledning er 0,5 m/s = 5m = 13 stedsskridt
% middel vandf?ring forede ledning er 0,62 m/s = 6,2m = 10 stedsskridt
% stor vandf?ring forede ledning 0,8 m/s = 8,2m = 8 stedsskridt


H2S_gasfase = [ 0.75.*koncentration_H2S_broend.*ones(8640,2) zeros(8640,9) ] ;
H2S_gasfase_ppm = [ 0.75.*(koncentration_H2S_broend./konversionsfaktor).*ones(8640,2) zeros(8640,9) ];
H2S_eq = 32.065.*10.^-3.*(H2S_gasfase_ppm./henrys_konstant);
H2S_frigivelse = ((KLa_h2s_matrix_g.*(H2S_vandfase - H2S_eq).*1.024^(temperatur-20)).*10); %Antages at for nedre rand er gasfasen nul ppm

H2S_oxidation_gm2s = [( 15.52.*(0.75.*(koncentration_H2S_broend./konversionsfaktor)).^0.55).*0.001*(1/3600).*ones(8640,2) zeros(8640,9)]; %[g S/m2*s]
H2S_oxidation = zeros(8640,11);

% Simulering af frigivelse af H2S og adsorption samt oxidation p? betonr?rets overflade ned gennem den forede ledning
for i=2:8640
    for j=2:11 % F?rste kolonne er en fiktiv nedre rand
        
        H2S_vandfase(i,j) = H2S_vandfase(i-1,j-1) - H2S_frigivelse(i,j-1); %[g S/m3]
                
        H2S_frigivelse(i,j) = KLa_h2s_matrix_g(i,j).*(H2S_vandfase(i,j)- H2S_eq(i,j) ).*1.024^(temperatur-20).*10; %[g S/m3]

        if (H2S_frigivelse(i,j) + 0.8381*H2S_gasfase(i-1,j) + 0.0809.*H2S_gasfase(i-1,j-1) ) > H2S_oxidation(i-1,j) %+ 0.1.*H2S_gasfase(i-1,j-1)
            
            H2S_gasfase(i,j) = H2S_frigivelse(i,j)  + 0.8381.*H2S_gasfase(i-1,j) + 0.0809.*H2S_gasfase(i-1,j-1) -  H2S_oxidation(i-1,j) ; %[g S/m3] + H2S_gasfase(i,j=2)
            
        else
           H2S_gasfase(i,j) = 0;
        end
        
        H2S_gasfase_ppm(i,j) = H2S_gasfase(i,j)./konversionsfaktor; %[ppm]
        
        H2S_oxidation_gm2s(i,j) = ( K_F(i,j).*H2S_gasfase_ppm(i,j).^n(i,j)) .*0.001*(1/3600); %[g S/m2*s]
        
        H2S_oxidation(i,j) = H2S_oxidation_gm2s(i,j).*(areal_gasfase_delfuld_g_matrix(i,j)./volume_gasfase_delfuld_g_matrix(i,j).*10); %[g S/m3]
        
        H2S_eq(i,j) = 32.065.*10.^-3.*(H2S_gasfase_ppm(i,j)./henrys_konstant(i,j)) ; %[g S/m3] 
        
    end
end

% Her laves matricer hvor nedre rand ikke er inkluderet
H2S_roerprofil_oxidation_vinter = H2S_oxidation(:, 2:end);
H2S_roerprofil_vandfase_vinter = H2S_vandfase(:, 2:end);
H2S_roerprofil_gas_vinter = H2S_gasfase(:, 2:end);
H2S_roerprofil_gas_ppm_vinter = H2S_gasfase_ppm(:, 2:end);
%% DEN FOREDE LEDNING - KORROSIONSRATE
% Faktor for hvor meget af syren der reagerer med betonen; antagelse: mellem 0.3-1
k_korrosion = 0.7; % [-] 
% N?r F stiger s? falder k ned til 0.3-0.4 - K_F??

% Alkalinitet af den p?f?rte m?rtel; antagelse: omkring 0.6-0.8, men normalt er den 0,2
A = 0.8; % [g CaCO3/g beton]

% Korrosionsrate i den forede ledning, angivet for hvert kontrolvolumen
korrosionsrate_vinter = 11.4.*k_korrosion.*((H2S_oxidation_gm2s.*3600)./A); %[mm/?r]
korrosionsrate_vinter_mean = mean(korrosionsrate_vinter(:,2:end))
%korrosionsrate_vinter_mean_best = mean(korrosionsrate_vinter(:,2:end));
%korrosionsrate_vinter_mean_worst = mean(korrosionsrate_vinter(:,2:end));
%% VARIABLER TIL WORST/BEST-CASE SCENARIER

nedstroems_ppm_vinter = H2S_roerprofil_gas_ppm_vinter(:,10);
opstroems_ppm_vinter = H2S_roerprofil_gas_ppm_vinter(:,1);

%nedstroems_ppm_vinter_best = H2S_roerprofil_gas_ppm_vinter(:,10);
%opstroems_ppm_vinter_best = H2S_roerprofil_gas_ppm_vinter(:,1);

%nedstroems_ppm_vinter_worst = H2S_roerprofil_gas_ppm_vinter(:,10);
%opstroems_ppm_vinter_worst = H2S_roerprofil_gas_ppm_vinter(:,1);
%% PLOTS
%PLOT AF DE TRE OPHOLDSTIDER
% figure(1)
% plot(opholdstid_lille/(60*60));
% hold on
% plot(opholdstid_mellem/(60*60));
% hold on
% plot(opholdstid_stor/(60*60));
% legend('15 m$^3$/d\o gn', '85 m$^3$/d\o gn', '160 m$^3$/d\o gn', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Opholdstid [time]', 'interpreter', 'latex')
% title('Opholdstid af det udpumpede vand fra trykledningen', 'fontSize', 14, 'interpreter', 'latex')
% hold off
% 
% %PLOT AF OPHOLDSTID OG SULFID KONCENTRATION
% figure(2)
% subplot(2,1,1)
% plot(opholdstid_mellem/(60*60), 'color', '#3FA663')
% title('Opholdstid af det udpumpede vand fra trykledningen', 'fontSize', 14, 'interpreter', 'latex')
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Opholdstid [time]', 'interpreter', 'latex')
% 
% subplot(2,1,2)
% plot(C_tryk, 'color', '#D43049')
% title('Sulfidkoncentration som funktion af opholdstid','fontSize', 14, 'interpreter', 'latex')
% xticks([0 1 2 3 4 5])
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Koncentration [g S/m$^3$]', 'interpreter', 'latex')
% %legend('Opholdstid', 'Sulfidkoncentration','location', 'southoutside', 'Orientation', 'horizontal')
% 
%PLOT AF SULDFIDKONCENTRATIONER SOM FUNKTION AF LEDNINGERNE
% figure(3)
% plot(C_tryk,'color', '#3FA663')
% hold on
% plot(C_PVC,'color', '#D43049')
% hold on
% plot(C_total,'color', '#0000e6')
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Koncentration [g S/m$^3$]', 'interpreter', 'latex')
% title('Sulfidkoncentrationer igennem r\o rforl\o bet', 'fontSize', 14, 'interpreter', 'latex')
% legend('Trykledningen', 'PVC-gravitationsledningen', 'Den forede ledning', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
% hold off
% 
% %PLOT AF FLOW I FORET LEDNING OG TRYKLEDNING
% figure(4)
% plot(flow_g_mellem*1000, 'color', '#3FA663')
% hold on
% plot(flow_tryk_mellem*1000, 'color', '#D43049')
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Vandf\o ring [L/s]', 'interpreter', 'latex')
% legend('Vandf\o ring i den forede ledning', 'Vandf\o ring i trykledningen', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
% title('Vandf\o ring i ledninger', 'fontSize', 14, 'interpreter', 'latex')
% hold off

%3D plot over koncentration H2S i v?skefase, tid og r?rprofil(kontrolvolumener)
% figure(5)
% surf((1:10), (1:8640), H2S_roerprofil_gas_ppm_vinter, 'FaceAlpha','0.9', 'EdgeColor','none')
% xlabel('R\o rprofil [m]', 'interpreter', 'latex')
% ylabel('Klokkesl\ae t', 'interpreter', 'latex')
% zlabel('Koncentration H$_{2}$S [ppm]', 'interpreter', 'latex')
% xticks([0 2 4 6 8 10 12])
% xticklabels({'0', '10', '20', '30', '40', '50', '60'})
% yticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% yticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% colormap(parula(8))
% colorbar
% 
% %PLOT OVER NEDSTR?MS OG OPSTR?MS SVOVLBRINTEKONC
% figure(6)
% nedstroems_ppm = H2S_roerprofil_gas_ppm_vinter(:,10);
% opstroems_ppm = H2S_roerprofil_gas_ppm_vinter(:,1);
% plot(nedstroems_ppm,'color', '#3FA663', 'LineWidth', 1.5)
% hold on
% plot(opstroems_ppm,  'color', '#D43049', 'LineWidth', 1.5)
% legend('Nedstr\o ms','Opstr\o ms', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
% xticks([0 720 1440 2160  2880  3600  4320  5040 5760 6480 7200 7920 8640])
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Koncentration H$_{2}$S [ppm]', 'interpreter', 'latex')
% title('Forl\o b af opstr\o ms og nedstr\o ms svovlbrintekocentrationer', 'interpreter', 'latex')
% hold off
% % 
% % KUMULERET PLOT OVER SVOVLBRINTEKONC
% 
% figure(7)
% sortet_nedstroems = sortrows(nedstroems_ppm);
% sortet_opstroems = sortrows(opstroems_ppm);
% fraktil_rank =( (1:length(nedstroems_ppm))' - 0.5) /length(nedstroems_ppm);
% plot(sortet_nedstroems,fraktil_rank, 'color',  '#D43049', 'LineWidth', 1.5)
% hold on
% plot(sortet_opstroems, fraktil_rank,  'color', '#3FA663', 'LineWidth', 1.5)
% legend('Nedstr\o ms','Opstr\o ms', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
% ylabel('Fraktion af tiden', 'interpreter', 'latex')
% xlabel('Koncentration H$_{2}$S [ppm]', 'interpreter', 'latex')
% axis([0 150 0 1])
% xticks([0 25 50 75 100 125 150])
% yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
% grid on
% hold off

% Plot over variation i trykledning
% figure(8)
% plot(variation_flow_tryk*1000, 'color', '#3FA663', 'Linewidth' , 2)
% xticks(0:720:8640)
% xticklabels({'00', '02', '04', '06', '08', '10', '12', '14', '16', '18', '20', '22', '00'})
% xlabel('Klokkesl\ae t', 'interpreter', 'latex')
% ylabel('Vandf\o ring [L/s]', 'interpreter', 'latex')
% axis([0 8640 0 2])
% title('Vandf\o ring til trykledningen', 'fontSize', 14, 'interpreter', 'latex')
