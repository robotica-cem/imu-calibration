%% Collects some data and calibrates imu using Fredrik's algorithm.

% Kjartan Halvorsen 2016-01-21

[t, accCalibData, gyrCalibData] = collect_imu_data(800,1);

save calibdata t accCalibData gyrCalibData
%load calibdata

[accCalibFro, accDataCalibratedFro] = Frosio_calib(t, accCalibData, 9);

[accCalib, gyrCalib, accDataCalibrated, gyrDataCalibrated] = FO_calib(t, accCalibData, gyrCalibData);

acc.t = t;
acc.dta = accCalibData';

gyr.t = t;
gyr.dta = gyrCalibData';

[accC, gyrC] = FO_apply_calib(acc, gyr, accCalib, gyrCalib);

h = plot_acc_gyro(acc,gyr);
h = plot_acc_gyro(accC,gyrC, h);


[t, accTest, gyrTest] = collect_imu_data(800,1, accCalibFro);

[t, accTest, gyrTest] = collect_imu_data(800,1, accCalib, gyrCalib);


