% regime_map_plots.m
%
% Script that can be used to gain a basic plot for the raw data underlying
% the various figures in our manuscript.
%
% Source:  https://github.com/OxfordFluidsLab/MovingPoolImpact
% Licence: GPL-3.0 (see the Git repo).
%
% T.C. Sykes (t.c.sykes@outlook.com)
% University of Warwick (2025)

% Read in the data from a csv file
t_out = readtable('postProcessing\paper_outcomesm.csv');

% Separate data by outcome
t_v = t_out(strcmp(t_out.outcome,'v'),:);
t_o = t_out(strcmp(t_out.outcome,'o'),:);


%% Figure 1a

% Setup figure
figure(1)
clf

% Plot
scatter(t_out.Re,t_out.We,10,'filled', 'Marker','v', 'MarkerFaceColor','r')

% Figure formatting
box on
xlabel('Reynolds number')   % Label x axis
ylabel('Weber number')      % Label y axis
title('Figure 1a')


%% Figure 4a

% Setup figure
figure(2)
clf

% Plot
scatter(t_o.quantity,t_o.Ca,10, ...
    'filled', 'Marker','o', 'MarkerFaceColor','r')
hold on
scatter(t_v.quantity,t_v.Ca,10, ...
    'filled', 'Marker','v', 'MarkerFaceColor','g')

% Figure formatting
box on
xlabel('Square root of pool-droplet velocity ratio')                                      % Label x axis
ylabel('Capillary number')    % Label y axis
title('Figure 4a')


