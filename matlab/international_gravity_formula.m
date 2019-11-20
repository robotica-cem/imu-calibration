function g = international_gravity_formula(latitude, altitude)
% g = international_gravity_formula(latitude, altitude)
% Returns the value of the gravitational acceleration, according to the
% international gravity formula
% 
% Input
%    latitude  - local latitude in degrees
%    altitude  - local altitude in meters above sea level

% Kjartan Halvorsen
% 2016-02-03

if nargin == 0
    disp('g at CEM')
    international_gravity_formula(19+25/60, 2360)
    disp('g in Uppsala')
    international_gravity_formula(59+52/60, 10)
    return
end

g = 9.7803327 * (1 + 5.3024e-3 * sind(latitude)^2  - 5.8e-6 * sind(2*latitude)^2 ) - 3.086e-6 * altitude;


