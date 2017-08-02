%%%%%%%%%%%%%%%%%%%%%%%%%%%% createCDMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates and returns steering vector(d_omega_theta)        %
% CDMA geometry and mics' location angles(fi) and delay times(tau)        %
% according to its parameters; number of Mics, radius of the CDMA and     %
% incident angle of raw source signal                                     %
% Author: Burak Parim                                                     %
% Date of start: 28.04.2017                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SENIC GmbH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_omega_theta,tau,phi] = createCDMA(numOfMics,radiusOfArray,incidentAngle)
    c = 343.6; % speed of sound at 21 deg. Celcius [m]
    M = numOfMics;
    r = radiusOfArray; % radius of the Array [m]
    theta = linspace(-pi,pi,200); % incident angle swipe for FR
    theta_s = incidentAngle; % angle between horizontal axis and normal of incident wavefront
    omega = linspace(-pi, pi, 200); % random angular angle CHANGE THIS! for frequency response
    phi = 2 * pi * ((1:M) - 1) / M; % angular position of the array elements [rad]
    tau = r * cos(theta_s - phi(:))/c; % time delays between mics and the origin of the array [sec.]
    inter_dist = 2 * r * sin(pi/M); %interelement distance of the array
        for k = 1:M
        steer_mic(k) = cos(theta_s - tau(k));
        end
    % create steering vector d(omega,teta) %% WARNING: This loop is not
    % able to adapt different number of microphones yet. --> FIX IT!!!
    for i = 1:length(omega)        
        for k = 1:M
            d_omega_theta(i,k) = exp(1j * omega(i) * r * (1/c) * steer_mic(k)).'; % steering vector describes the array geometry           
        end
    end
    d_omega_theta = d_omega_theta.'; %Transpose

end