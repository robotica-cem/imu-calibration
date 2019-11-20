function [accCalib, gyrCalib, accDataCalibrated, gyrDataCalibrated] = FO_calib(t, accData, gyrData)
% [accCalib, gyrCalib] = FO_calib(t, accData, gyroData))
% Runs Fredrik's calibration algorithm on provided data. Make sure the data
% starts with at least 3 seconds of the imu lying still. THen move the imu
% slowly to different orientations and keep it still for a few seconds in
% each different orientation.
%
% Input
%   t         -  vector (Nx1) of timestamps (seconds)
%   accData   -  matrix (Nx3) of acceleration data (m/s^2)
%   gyroData  -  matrix (Nx3) of gyro data (rad/s)
%
% Output
%   accCalib  -  Struct with the fields
%                       Dhat - calibration matrix
%                       bhat - sensor bias
%                       ghat - gravitational dir
%                       Rhat - estimated rotation matrix
%   gyrCalib  -  Struct with the fields
%                       bhat - sensor bias
%                       Rhat - estimated rotation matrix
%   accDataCalibrated  -  Struct with the fields
%                       t  - time vector
%                       dta  - data (3 x N)
%   gyrDataCalibrated  -  Struct with the fields
%                       t  - time vector
%                       dta  - data (3 x N)


%% Settings
g0 = 9.779; % Vertical (z-axis) gravitational acceleration (radius of the sphere that the data is calibrated to fit)
%g0 = 9.818; % Vertical (z-axis) gravitational acceleration (radius of the sphere that the data is calibrated to fit)

h = 0.01; % Differentiation step
alpha = 0.01; % Alpha parameter used in the line search algorithm
beta = 0.5; % Beta parameter used in the line search algorithm
max_iter = 100; % Maximum number of iterations allowed
eps1 = 0.00001; % Newton decrement minimum (Convergence threshold 1)
eps2 = 1e-8; % Minimum step length (Convergence treshold 2)
npar = 22; % Number of parameters
mc_iterations = 1; % Number of Monte Carlo iterations (set to 1 when using real data)

% Choose which parameters that are optimized using Gauss-Newton
% Calibration matrix, D - 1:9
% Accelerometer bias, ba - 10:12
% Gyroscope bias, bg - 13:15
% Accelerometer noise covariance matrix diagonal, Ra - 16:18
% Gyroscope noise covariance matrix diagonal, Rg - 19:21
% (Gravitational acceleration, g0 - 22 (not working))
unknown = 1:12; % Sets which parameters that are optimized
known = unknown(end)+1:npar; % Sets which parameters that are only initialized

% Create data structs. Let the first 3 seconds be init data used to
% estimate covariance matrices for the measurement noise
before3 = find(t < (t(1)+3));
after3 = find(t >= (t(1)+3));

dt = diff(t);
dt = [dt;dt(end)];

acc.t_init = t(before3);
acc.y_init = accData(before3,:)';

acc.t = t(after3);
acc.y = accData(after3,:)';
acc.dt = dt(after3);

gyr.t_init = t(before3);
gyr.y_init = gyrData(before3,:)';

gyr.t = t(before3);
gyr.y = gyrData(after3,:)';
gyr.dt = acc.dt;

N = length(acc.t);


%% Initialize result variables
z_results = zeros(npar,1);
q0_results = zeros(4,1);
q_results = zeros(4,N,1);
fval_results = zeros(1,1);
niter_results = zeros(1,1);
Hes_results = zeros(length(unknown),length(unknown),1);

%% Run calibration algorithm
% Orientation quaternion states
q = zeros(4,N); % Quaternions
q(1,:) = 1;
P = zeros(4,4,N); % Covariances

% Parameters
D = eye(3); % Accelerometer calibration matrix
ba = zeros(3,1); % Accelerometer bias vector
bg = zeros(3,1); % Gyroscope bias vector
Ra = eye(3); % Accelerometer covariance matrix
Rg = eye(3); % Gyroscope covariance matrix

% Initialization
% Initialize bg, Ra and Rg
% Estimate covariances and gyroscope bias
for i = 1:3
    bg(i) = mean(gyr.y_init(i,:));
    Ra(i,i) = var(acc.y_init(i,:));
    Rg(i,i) = var(gyr.y_init(i,:));
end

% Initialize q0
% True value
% q0 = gyr.quat(:,1);
% No knowledge (random initialization)
%     q0 = angle2quat(2*pi*randn(),2*pi*randn(),2*pi*randn(),'XYZ')'; %Random orientation initialization
% Estimate roll and pitch
r0meas = atan2(acc.y_init(2,:),sqrt(acc.y_init(3,:).^2+acc.y_init(1,:).^2));
p0meas = atan2(acc.y_init(1,:),sqrt(acc.y_init(3,:).^2+acc.y_init(2,:).^2));
q0 = angle2quat(mean(r0meas),mean(p0meas),0,'XYZ')';

%% Calibration loop
z0 = vectorizeParams(D,ba,g0,bg,Ra,Rg); % Vectorize the parameters
z = z0;
lam = 2*eps1 + 1;
niter = 0;
while lam/2 > eps1 && niter < max_iter
    % Approximate the Hessian
    [Jac,res,Stot,~,fval] = ekfFunc(z,q0,acc,gyr,h); % Run the EKFs and find the numeric Jacobian matrix
    Jac = Jac(:,unknown);
    z2 = z(known);
    z = z(unknown);
    
    gra = Jac'/Stot*res(:,1); % Gradient vector
    Hes = Jac'/Stot*Jac; % Hessian matrix
    
    % Calculate step direction
    sdir = -Hes\gra; %Search direction
    lamOld = lam;
    lam = gra'/Hes*gra;
    
    if lam/2 > eps1
        % Line search
        t = 1; %Step length
        zline = [z+t*sdir;z2];
        [~,~,~,~,linefval] = ekfFunc(zline,q0,acc,gyr,0);
        while linefval > fval + alpha*t*gra'*sdir %Wolfe condition
            t = beta*t;
            if t < eps2
                t = 0;
                lam = eps1;
                break;
            end
            zline = [z+t*sdir;z2];
            [~,~,~,~,linefval] = ekfFunc(zline,q0,acc,gyr,0);
        end
        
        % Obtain updated estimates
        z = z + t*sdir;
    end
    
    % Iterate until convergence
    z = [z;z2];
    niter = niter + 1;
    
    % Display interesting stuff
    disp([num2str(niter),' out of maximum ',num2str(max_iter),' iterations performed. fval = ',num2str(fval),'. lam/2 (Newton decrement) = ',num2str(lam/2),'. t (step length) = ',num2str(t),'.']);
    
end

% Return the calibration values
[accCalib.Dhat,accCalib.bhat,accCalib.ghat,gyrCalib.bhat,accCalib.Rhat,gyrCalib.Rhat] = structureParams(z);

acc.dta = acc.y;
gyr.dta = gyr.y;
%[accDataCalibrated, gyrDataCalibrated] = FO_apply_calib(acc,gyr, accCalib, gyrCalib);
accDataCalibrated = [];
gyrDataCalibrated = [];
end

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




function [D,ba,g,bg,Ra,Rg] = structureParams(z)

if numel(z) == 22
    D = reshape(z(1:9),3,3);
    ba = z(10:12);
    %     g = z(13);
    %     bg = z(14:16);
    %     Ra = diag(z(17:19));
    %     Rg = diag(z(20:22));
    bg = z(13:15);
    Ra = diag(z(16:18));
    Rg = diag(z(19:21));
    g = z(22);
else
    D = [];
    ba = [];
    g = [];
    bg = [];
    Ra = [];
    Rg = [];
end
end

function z = vectorizeParams(D,ba,g,bg,Ra,Rg)

if numel(D) == 9 && numel(ba) == 3 && numel(g) == 1 && numel(bg) == 3 && numel(Ra) == 9 && numel(Rg) == 9
    z = zeros(22,1);
    z(1:9) = D(:);
    z(10:12) = ba;
    %     z(13) = g;
    %     z(14:16) = bg;
    %     z(17:19) = diag(Ra);
    %     z(20:22) = diag(Rg);
    z(13:15) = bg;
    z(16:18) = diag(Ra);
    z(19:21) = diag(Rg);
    z(22) = g;
else
    z = [];
end
end


function [dD,dba,dg,dbg,dRa,dRg] = diffParams(D,ba,g,bg,Ra,Rg,h,i)

z = vectorizeParams(D,ba,g,bg,Ra,Rg);
if i <= length(z)
    z(i) = z(i) + h;
end
[dD,dba,dg,dbg,dRa,dRg] = structureParams(z);
end
