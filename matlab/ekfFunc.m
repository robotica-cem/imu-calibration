function [Jac,res,Stot,q,fval] = ekfFunc(z,q0,acc,gyr,h)
%% Initialize
npar = length(z);
N = length(acc.y);
[D,ba,g,bg,Ra,Rg] = structureParams(z);

q = zeros(4,N); % Quaternion states
P = zeros(4,4,N); % Quaternion state covariance matrices
res = zeros(3*(N-1),npar+1); % Residuals used to calculate numerical gradients
Jac = zeros(3*(N-1),npar); % Jacobian matrix
Stot = sparse(3*(N-1),3*(N-1)); % Residual covariance
fval = zeros(1,npar+1); % Value of the cost function

%% Run the EKF using current estimated parameters
q(1,:) = 1;
q(:,1) = q0;
P(:,:,1) = eye(4);
fval(1) = 0;
accCal = D\(acc.y - repmat(ba, 1, N)); % Calibrate accelerometer using current parameters
for k = 1:N-1
    k0 = (k-1)*3+1;
    % Measurement update
    accCal_k = accCal(:,k);
    [q(:,k), P(:,:,k), res(k0:k0+2,1), Stot(k0:k0+2,k0:k0+2), fval(1)] = acc_meas_update(q(:,k), P(:,:,k), accCal_k, Ra, g, fval(1));
    q(:,k) = q(:,k)/norm(q(:,k)); % Normalize
    
    % Time update
    gyrCal = gyr.y(:,k) - bg;
    [q(:,k+1), P(:,:,k+1)] = time_update(q(:,k), P(:,:,k), gyrCal, gyr.dt(k), Rg);
end

%% Find the numerical gradients
% Input h = 0 to skip Jacobian evaluation
if h > 0
    for i = 1:npar
        dq = zeros(4,N); % Quaternion states
        dq(1,:) = 1;
        dq(:,1) = q0;
        dP = zeros(4,4,N); % Quaternion state covariance matrices
        dP(:,:,1) = eye(4);
        [dD,dba,dg,dbg,dRa,dRg] = diffParams(D,ba,g,bg,Ra,Rg,h,i);
        accCal = dD\(acc.y - repmat(dba,1,N)); % Calibrate accelerometer using current parameters
        for k = 1:N-1
            k0 = (k-1)*3+1;
            % Measurement update
            [dq(:,k), dP(:,:,k), res(k0:k0+2,i+1), ~, fval(i+1)] = acc_meas_update(dq(:,k), dP(:,:,k), accCal(:,k), dRa, dg, fval(i+1));
            dq(:,k) = dq(:,k)/norm(dq(:,k)); % Normalize
            
            % Time update
            gyrCal = gyr.y(:,k) - dbg;
            [dq(:,k+1), dP(:,:,k+1)] = time_update(dq(:,k), dP(:,:,k), gyrCal, gyr.dt(k), dRg);
        end
        % Calculate numerical jacobian
        Jac(:,i) = (res(:,i+1) - res(:,1))/h;
    end
end

fval = fval(1); %Return the current value of the cost function
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

function [x, P, yres, S, fval] = acc_meas_update(x, P, acc, Ra, g, fval)

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

mval = (1/2)*(yres'/S*yres)^2 + log(det(S));
fval = fval + mval; % Update the cost function
% disp(['likelihood: ',num2str(exp(-mval))]);

end



