function [q] = ekfOrientation(q0,acc,gyr,Ra, Rg, g)
% 
% Quaternion-based ekf for estimating orientation of IMU. 
% Based on Fredrik Olsson's function ekfFunc, but simplified.
% The returned sequence of quaternions q is (4 x N). Each quaternion when
% operating on a spatial vector vs using quatrotate(q, vs) gives the vector
% in the body frame. 

% Kjartan Halvorsen
% 2016-01-26

%% Initialize
N = length(acc.dta);
q = zeros(4,N); % Quaternion states
P = zeros(4,4,N); % Quaternion state covariance matrices

%% Run the EKF
q(1,:) = 1;
q(:,1) = q0;
P(:,:,1) = eye(4);
for k = 1:N-1
    k0 = (k-1)*3+1;
    % Measurement update
    [q(:,k), P(:,:,k)] = acc_meas_update(q(:,k), P(:,:,k), acc.dta(:,k), Ra, g);
    q(:,k) = q(:,k)/norm(q(:,k)); % Normalize
    
    % Time update
    [q(:,k+1), P(:,:,k+1)] = time_update(q(:,k), P(:,:,k), gyr.dta(:,k), gyr.dt(k), Rg);
end
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

[Q,dQ] = quat2rotmat(x); % Q when operating on body vector gives spatial vector
H = dQ*g;
yres = acc - Q'*[0;0;g];
S = H*P*H' + Ra;
K = (P*H')/S;
x = x + K*yres;
P = (eye(size(K*H)) - K*H)*P;
end



