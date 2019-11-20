function acc_calibration(fname, modelorder)
%% Will load a file with imu data. Pick out the acceleration data and run the calibration method by
% Frosio et al ieee J Sensors
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6123172

if nargin == 0
    run_test;
    return
end

if nargin < 2
    modelorder = 2;
end

switch lower(fname[end-2:end])
    case 'log'
        [t, acc] = read_imu_logdata(fname);
    case 'mat'
        load(fname); % Should load t and acc
end
sfname = fname
sfname(end-2:end) = '.mat';
save(sfname, t, acc)
   
%figure(1)
%clf
%plot(t, acc)

[cframes] = choose_data([], acc)

calibdata = []
for i=1:2:length(cframes)-1
    calibdata = [calibdata; acc(cframes(i):cframes(i+1))];
end

model = get_model(modelorder);




keyboard

end %function

function run_test
    fname = 'sensorLog_20160118T172322.txt';

     acc_calibration(fname)
     
end %function

function model=get_model(modelorder)
    switch modelorder
        case 2
            model.residualfn = 'twoparam_residuals';
            model.calibratefn = 'twoparam_calibrate';
    end
end % function

function twoparam_residuals

end %function

function [tAcc, acc, tGyro, gyro, tMag, mag] = read_imu_logdata(fname)
    
    fid = fopen(fname, 'r')
    
    acc = [];
    gyro = [];
    mag = [];
    tAcc = [];
    tGyro = [];
    tMag = [];
    
    tline = fgets(fid);
    while ischar(tline)
        c = strsplit(tline,'\t');
        t = str2double(c{1})/1000; % The timestamp in s
        type = c{2};
        sample = str2double(c(3:5));
        switch upper(type)
            case 'ACC'
                acc = [acc;sample];
                tAcc = [tAcc;t];
            case 'GYR'
                gyro = [gyro;sample];
                tGyro = [tGyro;t];
            case 'MAG'
                mag = [mag;sample];
                tMag = [tMag;t];
        end % switch
        tline = fgets(fid);
    end %while
    
    fclose(fid);
end %function
    