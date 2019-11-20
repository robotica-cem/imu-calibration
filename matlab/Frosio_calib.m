function [accCalib, accDataCalibrated] = Frosio_calib(t, accData, modelorder)
% function [accCalib, accDataCalibrated] = Frosio_calib(t, accData, modelorder)
% Implementation of the calibration method by Frosio et al
% For the time being, only modelorder 6 and 9 is implemented, corresponding to
% the model
%  a_m = Sa + o, 
% where a_m is the measurement, a is the true acceleration, o is an offset
% (bias) and S is the sensitivity matrix. 
% In the model with 6 parameters S is diagonal.

% Kjartan Halvorsen
% 2016-01-21

%g = 9.818; % In Stockholm
g = 9.779; % In Mexico City

% Use matlab's nonlinear least squares solver lsqnonlin from the
% optimization toolbox

switch modelorder
    case 6
        fn = 'residuals6';
        % x = [S11, S22, S33, o1, o2, o3]
        x0 = [1,1,1,0,0,0]';
        x = lsqnonlin(@(x) residual6(x,accData',g), x0);
        accCalib.Dhat = diag(x(1:3));
        accCalib.bhat = x(4:6);
    case 9
        fn = 'residuals9';
        % x = [S11, S12, S13, S22, S23, S33, o1, o2, o3]
        x0 = [1,0,0, 1,0,1, 0,0,0]';
        x = lsqnonlin(@(x) residual9(x,accData',g), x0);
        accCalib.Dhat = [x(1), x(2), x(3); x(2), x(4), x(5); x(3), x(5), x(6)];
        accCalib.bhat = x(7:9);
    case 12
        fn = 'residuals12';
        % x = [S11, S12, S13, S21, S22, S23, s31, s32, S33, o1, o2, o3]
        x0 = [1,0,0, 0,1,0, 0,0,1, 0,0,0]';
        x = lsqnonlin(@(x) residual12(x,accData',g), x0);
        accCalib.Dhat = reshape(x(1:9), 3, 3)';
        accCalib.bhat = x(10:12);
    otherwise
        warning('Modelorder not supported')
            
end

% Apply the calibration parameters to the data
accDataCalibrated.t = t;
accDataCalibrated.dta = accCalib.Dhat \ (accData' - repmat(accCalib.bhat(:), 1, length(t)));
end

function r = residual6(x, accData, g)
    % Model
    %  ameasured = S*(atrue + o), 
    % with
    % S = diag(x(1:3));
    % o = x(4:6);
    anormsquared = ((1/x(1))*(accData(1,:)-x(4))).^2 + ((1/x(2))*(accData(2,:)-x(5))).^2 + ((1/x(3))*(accData(3,:) - x(6))).^2;
    r = anormsquared-g^2;
end

function r = residual9(x, accData, g)
    % Model
    %  ameasured = S*(atrue + o), 
    % with
    % S = [x(1) x(2) x(3) 
    %      x(2) x(4) x(5)
    %      x(3) x(5) x(6)];
    % o = x(7:9);
    S = [x(1) x(2) x(3);x(2) x(4) x(5); x(3) x(5) x(6)];
    o = x(7:9);
    a = S \ (accData - repmat(o, 1, size(accData,2)) );
    anormsquared = sum(a.^2);
    r = anormsquared-g^2;
end

function r = residual12(x, accData, g)
    % Model
    %  ameasured = S*(atrue + o), 
    % with
    % S = [x(1) x(2) x(3) 
    %      x(4) x(5) x(6)
    %      x(7) x(8)  x(9)];
    % o = x(10:12);
    S = reshape(x(1:9), 3,3)';
    o = x(10:12);
    a = S \ (accData - repmat(o, 1, size(accData,2)) );
    anormsquared = sum(a.^2);
    r = anormsquared-g^2;
end
