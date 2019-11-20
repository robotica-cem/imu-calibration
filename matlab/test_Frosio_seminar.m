%% Collects some data, calibrates imu using Frosio's method and checks the result.

% Kjartan Halvorsen 2019-11-20


% Collect calibration data
% Instructions: Let the phone rest on the table for 5 seconds. Pick it up
% and rotate it gently around all three axes for about 30 seconds
[tCalib, accCalibData, gyrCalibData] = collect_imu_data(35*50,1);
%load calibdata

% Remove missing data
[ia,j] = find(isnan(accCalibData));
[ig,j] = find(isnan(gyrCalibData));
i = union(ia,ig);
accCalibData(i,:) = [];
gyrCalibData(i,:) = [];
tCalib(i) = [];

%save calibdata tCalib accCalibData gyrCalibData

[accCalibFro6, accDataCalibratedFro6] = Frosio_calib(tCalib, accCalibData, 6);
[accCalibFro9, accDataCalibratedFro9] = Frosio_calib(tCalib, accCalibData, 9);
[accCalibFro12, accDataCalibratedFro12] = Frosio_calib(tCalib, accCalibData, 12);

% Take the first three seconds of data to be IMU at rest. Use this to
% calibrate for gyro bias.
before3 = find(tCalib < (tCalib(1)+3));
gyrCalib.bhat = mean(gyrCalibData(before3,:))';


% Check the calibration 
[t, accTest, gyrTest] = collect_imu_data(1500,1, accCalibFro9, gyrCalib);



