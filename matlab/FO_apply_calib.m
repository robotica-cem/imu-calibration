function [accData, gyrData] = FO_apply_calib(accRawData, gyrRawData, accCalib, gyrCalib)
%% Applies the calibration parameters in accCalib and gyrCalib to the data in accRawData and gyrRawData.
% 


    % Apply calibration
    
    accData = accRawData;
    gyrData = gyrRawData;

    N = length(accRawData.t);
    accData.dta = accCalib.Dhat\(accRawData.dta - repmat(accCalib.bhat, 1, N));
    gyrData.dta = gyrRawData.dta - repmat(gyrCalib.bhat, 1, N);
    
end