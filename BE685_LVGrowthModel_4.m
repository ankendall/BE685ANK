clear
clc

% plot of original data for comparison

trim_o = [0, 2, 3, 4];
h_o = [4.6, 4.8, 5.1, 4.5];

% parameters (biomarkers)

trim = [0, 2, 3, 4];
strain = [19.3, 19.5, 19.7, 19.0];
EDV = [108, 123, 126, 115];
ESV = [46.4, 48.7, 54.5, 48.6];
h0 = 4.6;
V0 = 46.4;

amaxg = 0.10; % max growth rate
amaxp = 0.15; % max regression rate
k1 = 1; % rate of growth
k2 = 60; % rate of regression
tm = 2.0; % midpoint of growth
tp = 4.0; % midpoint of regression
strain0 = 19.3; % non-pregnant strain control
h0 = 4.6; % non-pregnant ventricular wall thickness control
hmax = 5.1; % maximum measured ventricular wall thickness
EDV0 = EDV(1); % non-pregnant end diastolic volume control
ESV0 = ESV(1); % non-pregnant end systolic volume control

% calculating gr(t)
gr = zeros(size(trim));
gr(1) = 1;

for i = 2:length(trim)
    growth = 1/(1 + exp(-k1*(trim(i)-tm)));
    regress = 1/(1 + exp(-k2*(trim(i)-tp)));
    alphaeqn = strain(i) - strain0;
    growthrate = (alphaeqn * ((growth*amaxg) - (regress*amaxp)));
    gr(i) = gr(i-1) + growthrate * (trim(i) - trim(i-1));
end

% wall thickness
ht = h0 + (hmax - h0) * (gr - min(gr)) / (max(gr) - min(gr));

% calculating gc(t)

rEDV = (EDV/EDV0).^(1/3);
rESV = (ESV/ESV0).^(1/3);

rEDV0 = rEDV(1);
rESV0 = rESV(1);

gcEDV = rEDV / rEDV0;
gcESV = rESV / rESV0;

gcavg = (gcEDV + gcESV) / 2;


% plot

figure;
subplot(5,1,1)
plot(trim_o, h_o, 'LineWidth', 2, 'Color', [1,0,0]);
xlabel('Pregnancy Duration (Trimester)');
ylabel('LV Wall Thickness (mm)');
title('LV Wall Thickness Growth and Regression');

subplot(5,1,2)
plot(trim, gr, 'LineWidth', 2, 'Color', [1,0,1]);
xlabel('Pregnancy Duration (trimester)');
ylabel('Normalized Growth');
title('Growth Rate Over Time')

subplot(5,1,3)
plot(trim, ht, 'LineWidth', 2, 'Color', [0,1,0]);
xlabel('Pregnancy Duration (trimester)');
ylabel('Wall Thickness (mm)');
title('Model of LV Wall Thickness Growth and Regression');

subplot(5,1,4)
plot(trim, rEDV, 'LineWidth', 2);
hold on;
plot(trim, rESV, 'LineWidth', 2);
xlabel('Pregnancy Duration (trimester)');
ylabel('Radius (mm)')
legend('Radius from EDV', 'Radius from ESV');
title('Calculated Ventricular Radius from Volumes');

subplot(5,1,5)
plot(trim, (gcavg * rEDV0), 'LineWidth', 2, 'Color', [1,1,1]);
xlabel('Pregnancy Duration (trimester)');
ylabel('Estimated Cavity Radius (from ESV)');
title('Average Ventricular Radius Growth');