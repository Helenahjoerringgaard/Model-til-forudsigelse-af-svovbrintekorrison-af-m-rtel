%% Samlet korrosionsrate

korrosionrate_samlet_opstroems = (korrosionsrate_vinter_mean(:,1) + korrosionsrate_sommer_mean(:,1)+ korrosionsrate_foraar_mean(:,1) + korrosionsrate_efteraar_mean(:,1))/4;
korrosionrate_samlet_nedstroems = (korrosionsrate_vinter_mean(:,10) + korrosionsrate_sommer_mean(:,13)+ korrosionsrate_foraar_mean(:,12) + korrosionsrate_efteraar_mean(:,10))/4;
