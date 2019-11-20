function [d, v] = integrateAcc(acc, d0, v0, q, g)
% 
% Integrates acceleration to obtain displacement of IMU given already estimated orientation q. 

% Kjartan Halvorsen
% 2016-01-26

if nargin == 0
    unit_test();
    return
end

%% Initialize
N = length(acc.dt);

%% Apply rotation to rotate acceleration measurements to static frame

accS = quatrotate(quatinv(q'), acc.dta') + repmat(g(:)', N, 1); % N x 3


v = cumsum(accS).*repmat(acc.dt, 1, 3) + repmat(v0', N, 1);
d = cumsum(v).*repmat(acc.dt, 1, 3) + repmat(d0', N, 1);

end

function unit_test()

d0 = zeros(3,1);
v0 = zeros(3,1);
g = [0;0;9.82];

N = 20;
acc.dta = repmat(-g + [1;0;0], 1, N);
acc.dt = ones(N,1);

q = zeros(4,N);
q(1,:) = 1;

[d,v] = integrateAcc(acc, d0, v0, q, g);

figure
clf

plot(d(:,1))
hold on
plot(v(:,1), 'r')
end
