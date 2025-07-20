clear
clc

trim = [0, 2, 3, 4];
strain = [19.3, 19.5, 19.7, 19.0];
EDV = [108, 123, 126, 115];
ESV = [46.4, 48.7, 54.5, 48.6];
h0 = 4.6;
V0 = 46.4;

% parameters (biomarkers)

amaxg = 0.10869565; % max growth rate
k1 = 0.16667; % rate of growth
tm = 3; % midpoint of growth
strain0 = 19.3;
h0 = 4.6;
EDV0 = EDV(1);
ESV0 = ESV(1);

% calculating gr(t)
dt = trim(4) - trim(1);
gr = zeros(size(trim));
gr(1) = 1;

for i = 2:length(trim)
    growth = 1/(1 + exp(-k1*(trim(i)-tm)));
    alphaeqn = max(0, strain(i) - strain0);
    growthrate = alphaeqn * growth * amaxg;
    gr(i) = gr(i-1) + growthrate * dt;
end

% wall thickness
ht = h0 * gr;

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
subplot(3,1,1)
plot(trim, ht, 'LineWidth', 2);
xlabel('Pregnancy Duration (trimester)');
ylabel('Wall Thickness (mm)');
title('Model of LV Wall Thickness Growth');

subplot(3,1,2)
plot(trim, rEDV, 'LineWidth', 2);
hold on;
plot(trim, rESV, 'LineWidth', 2);
xlabel('Pregnancy Duration (trimester)');
ylabel('Radius (mm)')
legend('Radius from EDV', 'Radius from ESV');
title('Calculated Ventricular Radius from Volumes');

subplot(3,1,3)
plot(trim, (gcavg * rEDV0), 'LineWidth', 2);
xlabel('Pregnancy Duration (trimester)');
ylabel('Estimated Cavity Radius (from ESV)');
title('Average Ventricular Radius Growth');