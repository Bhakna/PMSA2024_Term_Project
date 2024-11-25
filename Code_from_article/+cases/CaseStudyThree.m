function [mu, sigma, Cs, T, init_transport_velocity, lung_dose_ug, plotFunction] = CaseStudyThree(i)

% Case study three : Dosage
init_transport_velocity = 3.6 / 100 * 60; % 3.6 mm/min - > dm/h
mu    = 3/2*1e-5; % [dm] (radius)
sigma = 0.6/2*1e-5; % [dm] (radius)
Cs = 100;
T = 36;


switch i
    case 1
        lung_dose_ug = 100; % [ug] per man (i.e. not per kg)

    case 2
        lung_dose_ug = 400; % [ug] per man (i.e. not per kg)
   
    case 3
        lung_dose_ug = 1600; % [ug] per man (i.e. not per kg)

    otherwise 
        error('i out of bounds')
        
end

plotFunction = @ PlotFunction;



function []  = PlotFunction(conc, mass, vmass, data, drug_data, lung_dose_ug,  x_lower, x_upper, t, i)

% Lung plots


c_a   = conc.art/drug_data.b_p*drug_data.fu;
c_v   = conc.vein/drug_data.b_p*drug_data.fu;

c_plot{1}   = [conc.f_u./c_a, conc.f_l./c_v]*drug_data.fu_f;
c_plot{2}   = [conc.ep_u./c_a, conc.ep_l./c_v]/drug_data.kp_u;
c_plot{3}   = [conc.sub_u./c_a, conc.sub_l./c_v]/drug_data.kp_u;

[~,ind] = max(t*60 > 5); % Remove singular start.
for k = 1:3
    c_plot{k} = c_plot{k}(ind:end,:);
end
t_sub = t(ind:end);



x_plot = [x_upper; x_lower];
titles = {'Fluid targeting', 'Epithelial targeting', 'Sub-epithelial targeting'};
subtitles = {'LDD = 100 microg', 'LDD = 400 microg', 'LDD = 1600 microg'};
z_label = 'Lung-targeting';

T_lim = 6;

figure(1)

for k = 1:3

    subplot(3,3,i + (k-1)*3)
    hold off
    surf(x_plot,t_sub,c_plot{k})
    hold on
    xlim([x_plot(1),x_plot(end)])
    ylim([0,T_lim])
    shading flat
    view([138.9 62]);
    
    if k == 1
        title(subtitles{i}, 'FontSize',8)
    end
    
    if i == 1
        zlabel(z_label, 'FontSize',8);
    end
    
    
    if (k == 1)
        m = [0.9, 0.57, 0.27];

            annotation('textbox', [0 m(i) 1.0 0.1], ...
                'String', titles{i}, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center')

    end
    
end


% Plasma/Lung
col = {'b','r','k'};
figure(2)

% semilogy(t,conc.vein,[col{i}, '--'])
% hold on
% C_lung = (mass.p_u + mass.p_l + mass.f_u + mass.f_l + ...
%             + mass.ep_u + mass.ep_l + mass.sub_u + mass.sub_l)/data.V.lung; 
% semilogy(t,C_lung,col{i});

semilogy(t,conc.vein/lung_dose_ug,[col{i}, '--'])
hold on

xlim([t(1),t(end)]);
ylim([1e-5, 1e-2]);
if i == 3
    for i = 1:3
        p2(i) = semilogy(0,1,col{i});
    end
   legend(p2, subtitles{:},'Location','southwest' );
   ylabel('Dose-normalized concentration', 'FontSize',8);
end


% Mass

labels = fields(mass);

figure(3)
subplot(1,3,i)
hold off
for k = 1:numel(labels)
    plot(t, mass.(labels{k}))
    hold on
end
title(subtitles{i}, 'FontSize',8)
xlim([t(1),t(end)])
ylim([0, mass.total(1)*1.1])
if i == 3
    legend(fields(mass))
end

% Sync z axis
if i == 3
    figure(1)
    S = [subplot(331), subplot(332), subplot(333);
        subplot(331+3), subplot(332+3), subplot(333+3);
        subplot(331+6), subplot(332+6), subplot(333+6);
        ];
    for j = 1:3
        zmax = max([S(j,1).CLim(2), S(j,2).CLim(2), S(j,3).CLim(2)]);
        for k = 1:3
           S(j,k).ZLim = [0,zmax];
           S(j,k).CLim = [0,zmax];
        end
    end
    
end