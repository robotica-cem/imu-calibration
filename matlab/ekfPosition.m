function [d] = ekfPosition(d0, P0, acc, q, Ra, Rg, g)
% 
% kf for estimating displacement of IMU given already estimated orientation q. 

% Kjartan Halvorsen
% 2016-01-26

%% Initialize
N = length(acc.dta);
d = zeros(3,N); % position
P = zeros(3,3,N); % Covariance matrices

%% Run the EKF
d(:,1) = d0;
P(:,:,1) = P0;
for k = 1:N-1
    % Measurement update
    [d(:,k), P(:,:,k)] = acc_meas_update(q(:,k), P(:,:,k), acc.dta(:,k), Ra, g);
    q(:,k) = q(:,k)/norm(q(:,k)); % Normalize
    
    % Time update
    [q(:,k+1), P(:,:,k+1)] = time_update(q(:,k), P(:,:,k), gyr.dta(:,k), gyr.dt(k), Rg);
end


function [x, P] = time_update(x, P, omega, dt, Rw)

w1 = omega(1); w2 = omega(2); w3 = omega(3);
Sw = [0  -w1  -w2  -w3;
    w1   0   w3  -w2;
    w2 -w3    0   w1;
    w3  w2  -w1    0];

q0 = x(1); q1 = x(2); q2 = x(3); q3 = x(4);
Sq = [-q1 -q2 -q3;
    q0 -q3  q2;
    q3  q0 -q1;
    -q2  q1  q0];

F = eye(4)+(dt/2)*Sw;
Q = (dt^2/4)*(Sq*Rw*Sq');

x = F*x;
P = F*P*F'+Q;

end

function [x, P] = acc_meas_update(x, P, acc, Ra, g)

q0 = x(1); q1 = x(2); q2 = x(3); q3 = x(4);
Q = [2*(q0^2+q1^2) - 1  2*(q1*q2-q0*q3)    2*(q1*q3+q0*q2);
    2*(q1*q2+q0*q3)    2*(q0^2+q2^2) - 1  2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2)    2*(q2*q3+q0*q1)    2*(q0^2+q3^2) - 1];

dQ = [-2*q2 2*q3 -2*q0 2*q1;
    2*q1 2*q0  2*q3 2*q2;
    4*q0    0     0 4*q3];

H = dQ*g;
yres = acc - Q'*[0;0;g];
S = H*P*H' + Ra;
K = (P*H')/S;
x = x + K*yres;
P = (eye(size(K*H)) - K*H)*P;
end



