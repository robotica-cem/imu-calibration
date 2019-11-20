%% Collects some data, calibrates imu using Frosio's method and checks the result.

% Kjartan Halvorsen 2016-01-26

% Collect calibration data
% Instructions: Let the phone rest on the table for 5 seconds. Pick it up
% and rotate it gently around all three axes
%[tCalib, accCalibData, gyrCalibData] = collect_imu_data(1200,1);
load calibdata

% Remove missing data
tCalib(find(isnan(tCalib))) = [];
[i,j] = find(isnan(accCalibData));
accCalibData(i,:) = [];
[i,j] = find(isnan(gyrCalibData));
gyrCalibData(i,:) = [];

%save calibdata tCalib accCalibData gyrCalibData

[accCalibFro6, accDataCalibratedFro6] = Frosio_calib(tCalib, accCalibData, 6);
[accCalibFro9, accDataCalibratedFro9] = Frosio_calib(tCalib, accCalibData, 9);
[accCalibFro12, accDataCalibratedFro12] = Frosio_calib(tCalib, accCalibData, 12);

% Take the first three seconds of data to be IMU at rest. Use this to
% calibrate for gyro bias.
before3 = find(tCalib < (tCalib(1)+3));
gyrCalib.bhat = mean(gyrCalibData(before3,:))';


% Collect some validation data. 
% Instructions: Let the phone rest on the table. Slide it on the table in a
% pure translation a known distance and back again.
%[tValid, accValidData, gyrValidData] = collect_imu_data(200,1);
load validdata
%save validdata tValid accValidData gyrValidData

acc.t = tValid;
acc.dta = accValidData';
acc.dt = diff(tValid);

gyr.t = tValid;
gyr.dta = gyrValidData';
gyr.dt = acc.dt;

% Apply the calibration and plot
[accC, gyrC] = apply_calib(acc, gyr, accCalibFro9, gyrCalib);
h = plot_acc_gyro(acc,gyr);
h = plot_acc_gyro(accC,gyrC, h);

% Find orientation
% Initial orientation quaternian q0 and g
% Estimate roll and pitch
before3 = find(tValid < (tValid(1)+3));
r0meas = atan2(accC.dta(2,before3),sqrt(accC.dta(3,before3).^2+accC.dta(1,before3).^2));
p0meas = atan2(accC.dta(1, before3),sqrt(accC.dta(3,before3).^2+accC.dta(2,before3).^2));
q0 = angle2quat(mean(r0meas),mean(p0meas),0,'XYZ')';
g = -mean(accC.dta(:,before3), 2);
Ra = eye(3); Rg=eye(3);
for i=1:3
    Ra(i,i) = var(accC.dta(i,before3));
    Rg(i,i) = var(gyrC.dta(i, before3));
end
q = ekfOrientation(q0,accC,gyrC,Ra, Rg, norm(g));

d0 = zeros(3,1); v0=zeros(3,1);
[d, v] = integrateAcc(accC, d0, v0, q, g);

figure(1)
clf
plot(accC.t, d);





[t, accTest, gyrTest] = collect_imu_data(1500,1, accCalibFro12, gyrCalib);



