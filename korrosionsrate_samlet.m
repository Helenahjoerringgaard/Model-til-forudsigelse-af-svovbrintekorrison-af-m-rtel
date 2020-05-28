%% Samlet korrosionsrate

korrosionrate_samlet_opstroems = (korrosionsrate_vinter_mean(:,1) + korrosionsrate_sommer_mean(:,1)+ korrosionsrate_foraar_mean(:,1) + korrosionsrate_efteraar_mean(:,1))/4
korrosionrate_samlet_nedstroems = (korrosionsrate_vinter_mean(:,10) + korrosionsrate_sommer_mean(:,13)+ korrosionsrate_foraar_mean(:,12) + korrosionsrate_efteraar_mean(:,10))/4

korrosionrate_samlet_opstroems_best = (korrosionsrate_vinter_mean_best(:,1) + korrosionsrate_sommer_mean_best(:,1)+ korrosionsrate_foraar_mean_best(:,1) + korrosionsrate_efteraar_mean_best(:,1))/4
korrosionrate_samlet_nedstroems_best = (korrosionsrate_vinter_mean_best(:,8) + korrosionsrate_sommer_mean_best(:,10)+ korrosionsrate_foraar_mean_best(:,9) + korrosionsrate_efteraar_mean_best(:,9))/4

korrosionrate_samlet_opstroems_worst = (korrosionsrate_vinter_mean_worst(:,1) + korrosionsrate_sommer_mean_worst(:,1)+ korrosionsrate_foraar_mean_worst(:,1) + korrosionsrate_efteraar_mean_worst(:,1))/4
korrosionrate_samlet_nedstroems_worst = (korrosionsrate_vinter_mean_worst(:,12) + korrosionsrate_sommer_mean_worst(:,13)+ korrosionsrate_foraar_mean_worst(:,12) + korrosionsrate_efteraar_mean_worst(:,12))/4

samlet_opstroems_worst = sortrows((opstroems_ppm_vinter_worst + opstroems_ppm_foraar_worst + opstroems_ppm_sommer_worst + opstroems_ppm_efteraar_worst)/4);
samlet_nedstroems_worst = sortrows((nedstroems_ppm_vinter_worst + nedstroems_ppm_foraar_worst + nedstroems_ppm_sommer_worst + nedstroems_ppm_efteraar_worst)/4);

samlet_opstroems_best = sortrows((opstroems_ppm_vinter_best + opstroems_ppm_foraar_best + opstroems_ppm_sommer_best + opstroems_ppm_efteraar_best)/4);
samlet_nedstroems_best = sortrows((nedstroems_ppm_vinter_best + nedstroems_ppm_foraar_best + nedstroems_ppm_sommer_best + nedstroems_ppm_efteraar_best)/4);

samlet_opstroems = sortrows((opstroems_ppm_vinter + opstroems_ppm_foraar + opstroems_ppm_sommer + opstroems_ppm_efteraar)/4);
samlet_nedstroems = sortrows((nedstroems_ppm_vinter + nedstroems_ppm_foraar + nedstroems_ppm_sommer + nedstroems_ppm_efteraar)/4);

fraktil_rank =((1:8640)' - 0.5) /8640;

% %PLOT OVER FORSKELLIGE SCENARIER
figure(9)
set(figure(9),'defaultAxesTickLabelInterpreter','latex')
plot(samlet_opstroems,fraktil_rank, 'color',  '#3FA663', 'LineWidth', 1.5)
hold on
plot(samlet_nedstroems, fraktil_rank,  'color', '#D43049', 'LineWidth', 1.5)
hold on
plot(samlet_opstroems_worst,fraktil_rank, '--', 'color',  '#3FA663', 'LineWidth', 1.5)
hold on
plot(samlet_nedstroems_worst, fraktil_rank,'--',  'color', '#D43049', 'LineWidth', 1.5)
hold on
plot(samlet_opstroems_best,fraktil_rank,':', 'color',  '#3FA663', 'LineWidth', 1.5)
hold on
plot(samlet_nedstroems_best, fraktil_rank,':',  'color', '#D43049', 'LineWidth', 1.5)
legend('Reference (Opstr\o ms)','Reference (Nedstr\o ms)','Worst-case (Opstr\o ms)', 'Worst-case (Nedstr\o ms)','Best-case (Opstr\o ms)', 'Best-case (Nedstr\o ms)', 'location', 'southeast', 'interpreter', 'latex', 'fontsize', 10)
ylabel('Procentdel af tiden [\%]', 'interpreter', 'latex', 'fontsize', 12)
xlabel('Svovlbrintekoncentration [ppm]', 'interpreter', 'latex', 'fontsize', 12)
axis([0 400 0.7 1])
xticks(0:50:400)
yticks(0.7:0.05:1)
yticklabels({ '70', '75', '80','85', '90', '95', '100'} )
title('Kumuleret fordeling af svovlbrintekoncentrationer afh\ae ngig af scenarie','interpreter', 'latex', 'fontsize', 14)
grid on
hold off

nedstroems_ppm_sommer = H2S_roerprofil_gas_ppm_sommer(:,13);
opstroems_ppm_sommer = H2S_roerprofil_gas_ppm_sommer(:,1);

nedstroems_ppm_efteraar = H2S_roerprofil_gas_ppm_efteraar(:,10);
opstroems_ppm_efteraar = H2S_roerprofil_gas_ppm_efteraar(:,1);

nedstroems_ppm_foraar = H2S_roerprofil_gas_ppm_foraar(:,12);
opstroems_ppm_foraar = H2S_roerprofil_gas_ppm_foraar(:,1);

nedstroems_ppm_vinter = H2S_roerprofil_gas_ppm_vinter(:,10);
opstroems_ppm_vinter = H2S_roerprofil_gas_ppm_vinter(:,1);

nedstroems_samlet = ( nedstroems_ppm_sommer + nedstroems_ppm_vinter + nedstroems_ppm_foraar + nedstroems_ppm_efteraar )/4;
opstroems_samlet = ( opstroems_ppm_vinter + opstroems_ppm_sommer + opstroems_ppm_foraar + opstroems_ppm_efteraar )/4;

figure(7)
set(figure(7),'defaultAxesTickLabelInterpreter','latex'); 
sortet_nedstroems = sortrows(nedstroems_samlet);
sortet_opstroems = sortrows(opstroems_samlet);
fraktil_rank =( (1:length(nedstroems_samlet))' - 0.5) /length(nedstroems_samlet);
plot(sortet_nedstroems,fraktil_rank, 'color',  '#D43049', 'LineWidth', 1.5)
hold on
plot(sortet_opstroems, fraktil_rank,  'color', '#3FA663', 'LineWidth', 1.5)
legend('Nedstr\o ms','Opstr\o ms', 'location', 'southoutside', 'Orientation', 'horizontal', 'interpreter', 'latex')
ylabel('Procentdel af tiden [\%]', 'interpreter', 'latex')
xlabel('Svovlbrintekoncentration [ppm]', 'interpreter', 'latex', 'fontsize', 12)
axis([0 350 0.8 1])
xticks([0:25:350])
yticks([0.8:0.05:1])
yticklabels({'80','85', '90', '95', '100'})
title('Kumuleret fordeling af svovlbrintekoncentrationer','interpreter', 'latex', 'fontsize', 14)
grid on
hold off