%% clearing
clear; close all; clc;

%% create blue2red colormap
% Define the anchor colors [R, G, B]
colors = [
    0, 0, 1;   % Blue
    1, 1, 1;   % White
    1, 0, 0    % Red
];

% Interpolate to create a smooth 256-step gradient
n = 256;
blue2red_map = interp1(linspace(0, 1, size(colors, 1)), colors, linspace(0, 1, n));
% Apply it to future figures
colormap(blue2red_map);
close

%% Load data
load DNIDEP.mat
data = DNIDEP.Obs;

%% print DNI for months of interest (for reference) and create variables for 2022 and 2023
% 1/15/2022	1.6186
% 2/15/2022	1.4094
% 3/15/2022	-0.7017
% 4/15/2022	-0.404
% 5/15/2022	0.0725
dni_2022 = [1.6186, 1.4094, -0.7017, -0.404, 0.0725];

% 1/15/2023	1.5475
% 2/15/2023	0.3561
% 3/15/2023	1.5711
% 4/15/2023	1.8446
% 5/15/2023	1.9782
dni_2023 = [1.5475, 0.3561, 1.5711, 1.8446, 1.9782];

dni_vals = [dni_2022; dni_2023];
csvwrite('dni_vals.csv', dni_vals);

%% make histogram and color it
figure
nbins = 40;
[counts, edges] = histcounts(data, nbins);
centers = edges(1:end-1) + diff(edges)/2;
colormap(blue2red_map)
b = bar(centers, counts, 'BarWidth', 1);
b.FaceColor = 'flat';
b.CData = centers'; % Note the transpose to match the number of bars

ax = gca;
ax.YLabel.FontSize = 12;
set(ax,'YLim',[-5 60]);
%set(ax,'YTickLabel', {'','0','10','20','30','40','50','60'})
%ylabel('Frequency')

ax.XLabel.FontSize = 12;
%strX = 'Dakar Ni$$\tilde{n}$$o/Ni$$\tilde{n}$$a Index $(^\circ C)$';
strX = 'Dakar Ni$$\tilde{n}$$o/Ni$$\tilde{n}$$a Index ($^\circ$C)';
xlabel(strX,'Interpreter','latex');

hold on
y_height = 0;   % How high the line should go
y_lowerLimit = -1.8; % How low the line should go
symbol_height = -2.35;
triangle_adjustment = -0.1;
circle_adjustment = 0.1;
line_adjustment = 0.3;

% Insert markers for 2023
x_pos = dni_2023(3); % March
line([x_pos, x_pos], [y_lowerLimit, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height, 'ks', 'MarkerFaceColor','k')

x_pos = dni_2023(4); % April
line([x_pos, x_pos], [y_lowerLimit, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height+triangle_adjustment, 'k^', 'MarkerFaceColor','k')

x_pos = dni_2023(5); % May
line([x_pos, x_pos], [y_lowerLimit, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height+circle_adjustment, 'ko', 'MarkerFaceColor','k')

print 'DNIhist_2023.png' -dpng -r300
print 'DNIhist_2023.tiff' -dtiff -r300

% Insert markers for 2022
x_pos = dni_2022(3); % March
line([x_pos, x_pos], [y_lowerLimit, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height, 'ks')

x_pos = dni_2022(4); % April
line([x_pos, x_pos], [y_lowerLimit, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height+triangle_adjustment, 'k^')

x_pos = dni_2022(5); % May
line([x_pos, x_pos], [y_lowerLimit + line_adjustment, y_height], ...
    'Color', 'black', ...
    'LineWidth', 1, ...
    'LineStyle', '-');
plot(x_pos, symbol_height+circle_adjustment, 'ko')

print 'DNIhist_2022_2023.png' -dpng -r300
print 'DNIhist_2022_2023.tiff' -dtiff -r300


%% Calculate "percentage below" using all data (all months) for each of the threshold values 

data = DNIDEP.Obs;
nbins = 40;

[m n] = size(dni_vals);
percentile_simple = zeros(m, n)*nan;
percentile_interp = zeros(m, n)*nan;

ctr=0
f = figure;
set(f,'OuterPosition',[50 50 1800 1000]);
yr = [2022, 2023]
for i = 1:m %loop about year, n=2

    for j = 1:n %loop about month, m=5
    ctr = ctr+1;

    [counts, edges] = histcounts(data, nbins);
    threshold = dni_vals(i,j)
    
    % Method A: Summing full bins (Simple Approach)
    % Find bins where the right (upper) edge is less than or equal to the threshold
    bins_fully_below = edges(2:end) <= threshold;
    counts_below = sum(counts(bins_fully_below));
    total_counts = sum(counts);
    percentage_simple = (counts_below / total_counts) * 100;
    
    % Method B: Interpolation (More Accurate)
    % This accounts for the partial bin if the threshold falls between edges
    bin_widths = diff(edges);
    cumulative_counts = [0, cumsum(counts)];
    % Interpolate the cumulative distribution at the threshold
    count_interp = interp1(edges, cumulative_counts, threshold, 'linear', 'extrap');
    percentage_interp = (count_interp / total_counts) * 100;
    
    % Visualization
    subplot(2,5,ctr);
    histogram(data, nbins);
    xline(threshold, 'r', 'LineWidth', 2, 'Label', 'Threshold');
    title([num2str(yr(i)) '-AllMonths %blw ' num2str(threshold) ': ' num2str(percentage_interp, '%.2f') '%']);
    
    % Output results to Command Window
    fprintf('Total points: %d\n', total_counts);
    fprintf('Percentage below %.2f (Simple): %.2f%%\n', threshold, percentage_simple);
    fprintf('Percentage below %.2f (Interpolated): %.2f%%\n', threshold, percentage_interp);
    
    percentile_simple(i,j) = percentage_simple;
    percentile_interp(i,j) = percentage_interp;

    end
end
    
eval(['print ', '''threshpct_AllMonths.png''', ' -r300 -dpng'])


%csvwrite('percentile_simple.csv', percentile_simple)
csvwrite('percentile_interp.csv', percentile_interp)

%% Calculate "percentage below" using all data (within months) for each of the threshold values 

% create sub-datasets for each month
Ijan = find(month(DNIDEP.time)==1);
Ifeb = find(month(DNIDEP.time)==2);
Imar = find(month(DNIDEP.time)==3);
Iapr = find(month(DNIDEP.time)==4);
Imay = find(month(DNIDEP.time)==5);

DNI_jan = DNIDEP(Ijan,:);
DNI_feb = DNIDEP(Ifeb,:);
DNI_mar = DNIDEP(Imar,:);
DNI_apr = DNIDEP(Iapr,:);
DNI_may = DNIDEP(Imay,:);


[m n] = size(dni_vals);
percentile_simple_month = zeros(m, n)*nan;
percentile_interp_month = zeros(m, n)*nan;

ctr=0
f = figure;
set(f,'OuterPosition',[50 50 1800 1000]);
yr = [2022, 2023]
for i = 1:m %loop about year, n=2

    for j = 1:n %loop about month, m=5

        ctr = ctr+1

        %select the correct month, 1-5, Jan - May
        if i==1
            data = DNI_jan.Obs;
        elseif j==2
            data = DNI_feb.Obs;
        elseif j==3
            data = DNI_mar.Obs;
        elseif j==4
            data = DNI_apr.Obs;
        elseif j==5
            data = DNI_may.Obs;
        end


        [counts, edges] = histcounts(data, nbins);
        threshold = dni_vals(i,j)
        
        % Method A: Summing full bins (Simple Approach)
        % Find bins where the right (upper) edge is less than or equal to the threshold
        bins_fully_below = edges(2:end) <= threshold;
        counts_below = sum(counts(bins_fully_below));
        total_counts = sum(counts);
        percentage_simple = (counts_below / total_counts) * 100;
    
        % Method B: Interpolation (More Accurate)
        % This accounts for the partial bin if the threshold falls between edges
        bin_widths = diff(edges);
        cumulative_counts = [0, cumsum(counts)];
        % Interpolate the cumulative distribution at the threshold
        count_interp = interp1(edges, cumulative_counts, threshold, 'linear', 'extrap');
        percentage_interp = (count_interp / total_counts) * 100;
    
        % Visualization
        subplot(2,5,ctr);
        histogram(data, nbins);
        xline(threshold, 'r', 'LineWidth', 2, 'Label', 'Threshold');
        title([num2str(yr(i)) '-' num2str(j) ' %blw ' num2str(threshold) ': ' num2str(percentage_interp, '%.2f') '%']);
    
        % Output results to Command Window
        fprintf('Total points: %d\n', total_counts);
        fprintf('Percentage below %.2f (Simple): %.2f%%\n', threshold, percentage_simple);
        fprintf('Percentage below %.2f (Interpolated): %.2f%%\n', threshold, percentage_interp);
        
        percentile_simple_month(i,j) = percentage_simple;
        percentile_interp_month(i,j) = percentage_interp;

    end
end

eval(['print ', '''threshpct_months.png''', ' -r300 -dpng'])

%csvwrite('percentile_simple.csv', percentile_simple)
csvwrite('percentile_interp_month.csv', percentile_interp_month)


