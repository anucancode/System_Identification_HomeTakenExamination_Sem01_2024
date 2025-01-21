%%
clear
clear all
clc

%% reading the excel file
data = xlsread('DataSet_2132290.xls');

%%
x = data(:, 1);
y = data(:, 2);

%% plotting the data
figure('Name','Temperature vs Boiling point Elevation');
plot(x,y,'o');
xlabel('Temperature (C)'); ylabel('Boiling point (K)')
grid on 

%% Create and Plot a Selection of Polynomials and exponential eqn
[boilpointdata2,gof2] = fit(x,y,'poly2');
[boilpointdata3,gof3] = fit(x,y,'poly3','Normalize','on');
[boilpointdata4,gof4] = fit(x,y,'poly4','Normalize','on');
[boilpointdata5,gof5] = fit(x,y,'poly5','Normalize','on');
[boilpointdata6,gof6] = fit(x,y,'poly6','Normalize','on');
[boilpointdataExp,gofExp] = fit(x,y,'exp1');

figure;
hold on
plot(boilpointdata2,x,y);
plot(boilpointdata3,'b');
plot(boilpointdata4,'g');
plot(boilpointdata5,'m');
plot(boilpointdata6,'b--');
plot(boilpointdataExp,'r--');
hold off
legend('temperature v boiling point','poly2','poly3','poly4','poly5','poly6','exp1', ...
    'Location','NorthWest');
title('Temperature vs Boiling Point Elevation');
xlabel('Temperature (C)');
ylabel('Boiling Point (K)');
grid on

%% plot residuals

% Plot residuals for boilpointdata2
figure; % Create a new figure
plot(boilpointdata2, x, y, 'residuals');
title('Residuals for Poly2');
xlabel('Temperature');
ylabel('Residuals');

% Plot residuals for boilpointdata3
figure; % Create a new figure
plot(boilpointdata3, x, y, 'residuals');
title('Residuals for Poly3');
xlabel('Temperature');
ylabel('Residuals');

% Plot residuals for boilpointdata4
figure; % Create a new figure
plot(boilpointdata4, x, y, 'residuals');
title('Residuals for Poly4');
xlabel('Temperature');
ylabel('Residuals');

% Plot residuals for boilpointdata5
figure; % Create a new figure
plot(boilpointdata5, x, y, 'residuals');
title('Residuals for Poly5');
xlabel('Temperature');
ylabel('Residuals');

% Plot residuals for boilpointdata6
figure; % Create a new figure
plot(boilpointdata6, x, y, 'residuals');
title('Residuals for Poly6');
xlabel('Temperature');
ylabel('Residuals');

% Plot residuals for boilpointdataExp
figure; % Create a new figure
plot(boilpointdataExp, x, y, 'residuals');
title('Residuals for Exponential Fit');
xlabel('Temperature');
ylabel('Residuals');

%% Examine Fits Beyond the Data Range

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata2); 
hold off;
title('Fit for Poly2 Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata3); 
hold off;
title('Fit for Poly3 Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata4); 
hold off;
title('Fit for Poly4 Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata5); 
hold off;
title('Fit for Poly5 Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata6);
hold off;
title('Fit for Poly6 Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdataExp); 
hold off;
title('Fit for Exponential Model Beyond Data Range');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

%% Plot Prediction Intervals

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata2, 'predobs'); 
hold off;
title('Prediction Intervals for Poly2');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata3, 'predobs'); 
hold off;
title('Prediction Intervals for Poly3');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata4, 'predobs');
hold off;
title('Prediction Intervals for Poly4');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata5, 'predobs'); 
hold off;
title('Prediction Intervals for Poly5');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdata6, 'predobs'); 
hold off;
title('Prediction Intervals for Poly6');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

figure;
plot(x, y, 'o'); 
xlim([100, 300]); 
hold on;
plot(boilpointdataExp, 'predobs'); 
hold off;
title('Prediction Intervals for Exponential Model');
xlabel('X-Axis');
ylabel('Y-Axis');
grid on

%% goodness of fit characteristics

gof2
gof3
gof4
gof5
gof6
gofExp

%%

boilpointdata2
boilpointdata3
boilpointdata4
boilpointdata5
boilpointdata6

%%
ci2 = confint(boilpointdata2)
ci3 = confint(boilpointdata3)
ci4 = confint(boilpointdata4)
ci5= confint(boilpointdata5)
ci6 = confint(boilpointdata6)

%%
xFuture = (180:10:250).';
yFuture = boilpointdata2(xFuture)
ci2 = predint(boilpointdata2,xFuture,0.95,'observation')

%%
plot(x,y,'o');
xlim([100,300])
hold on
plot(boilpointdata2)
h = errorbar(xFuture,yFuture,yFuture-ci2(:,1),ci2(:,2)-yFuture,'.');
hold off
legend('temp v bp','poly2','prediction', ...
    'Location','NorthWest')
grid on

%% xdesign matrix

n = 2; % for a quadratic model
x_design = ones(length(x), n + 1); % Initialize the design matrix with ones
for i = 1:n
    x_design(:, i + 1) = x.^i;
end

%%
XTX = x_design' * x_design; % X_transpose times X
XTX_inv = inv(XTX); % Inverse of X_transpose times X
XTY = x_design' * y; % X_transpose times Y

%% Least squares estimator
p_hat = XTX_inv * XTY

y_est = p_hat(1).*x.^2+p_hat(2).*x+p_hat(3);

% Plotting OLS estimation fit
figure('Name','OLS estimation fit');
scatter(x,y,'o','filled');
hold on
plot(x,y_est,'r');
legend('Data','Parameter estimation');
title('OLS estimation plot for fitted curve');
xlabel('Temperature [deg C]');
ylabel('Boiling point [K]');
grid on;
%%
residuals = y - x_design * p_hat
RMS = sum(residuals.^2) / (length(y) - length(p_hat));

%% Variance covariance matrix
var_cov_matrix = RMS * inv(x_design' * x_design)
%%
% Assuming residuals are already computed
RSS = sum(residuals.^2)
% Number of observations
n = length(y);
% Number of parameters estimated (including the intercept)
m = length(p_hat);
% Estimating sigma_e^2
sigma_e_squared = RSS / (n - m)

%%
variances = diag(var_cov_matrix)
 
std_devs = sqrt(variances)
%% SVD
[U, S, V] = svd(x_design);
%%
[U, S, V] = svd(x_design, 'econ');
sigma_inv = diag(1./diag(S));
p_hat_svd = V * sigma_inv * U' * y

y_est_svd = p_hat_svd(1).*x.^2+p_hat_svd(2).*x+p_hat_svd(3);

% Plotting SVD estimation fit
figure('Name','SVD estimation fit');
scatter(x,y,'o','filled');
hold on
plot(x,y_est,'r');
legend('Data','Parameter estimation');
title('SVD estimation plot for fitted curve');
xlabel('Temperature [deg C]');
ylabel('Boiling point [K]');
grid on;