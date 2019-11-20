%% Collects some data, calibrates imu using Fredrik's method and checks the result.

% Kjartan Halvorsen 2016-02-09

% Collect calibration data
% Instructions: Let the phone rest on the table for 5 seconds. Pick it up
% and rotate it gently around all three axes
[tCalib, accCalibData, gyrCalibData] = collect_imu_data(2000,1);
%load calibdata

% Remove missing data
tCalib(find(isnan(tCalib))) = [];
[i,j] = find(isnan(accCalibData));
accCalibData(i,:) = [];
[i,j] = find(isnan(gyrCalibData));
gyrCalibData(i,:) = [];

%save calibdata tCalib accCalibData gyrCalibData
[accCalibFO, gyrCalibFO, accDataCalibrated, gyrDataCalibrated] = FO_calib(tCalib, accCalibData, gyrCalibData);
[accCalibFro9, accDataCalibratedFro9] = Frosio_calib(tCalib, accCalibData, 9);

% Take the first three seconds of data to be IMU at rest. Use this to
% calibrate for gyro bias.
before3 = find(tCalib < (tCalib(1)+3));
gyrCalib.bhat = mean(gyrCalibData(before3,:))';


% Check the calibration 
%[t, accTest, gyrTest] = collect_imu_data(2000,1, accCalibFro9, gyrCalib);
[t, accTest, gyrTest] = collect_imu_data(2000,1, accCalibFO, gyrCalibFO);



