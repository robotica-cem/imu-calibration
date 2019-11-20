function [acc, gyro, t] = imu_receive_acc_gyro()

    Y = judp('RECEIVE',5555,200);
    Y(Y==44) = ' '; % Espacio como separador.

    temp  = cell2mat(textscan(char(Y),'%f'));
    
    % Check length of UDP packet
    if length(temp) < 9
        warning('Too few data in UDP packet. Returning empty matrices')
        acc = [];
        gyro = [];
        t = [];
        return
    end
    
    t = temp(1); % Timestamp in seconds
    
    % Check sensor id
    switch temp(2)
        case 3
            acc = temp(3:5);
        case 4
            gyro = temp(3:5);
    end
    
    switch temp(6)
        case 3
            acc = temp(7:9);
        case 4
            gyro = temp(7:9);
    end
    
    %acc
end
       