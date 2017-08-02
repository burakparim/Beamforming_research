%%%%%%%%%%%%%%%%%%%%%%% SENIC Lamp/Hub CDMA Design %%%%%%%%%%%%%%%%%%%%%%%%
% On this file simulations and evaluation of CDMA and ACDMA are runned    %
% Date of Start: 27.04.2017                                               %
% Author: Burak Parim                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SENIC GmbH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
M = 6; % Number of Mics on the array
r = 0.01025; % radius of circular Mic. Array [m]
incidentAngle = pi/2; % DOA given parameter (DULL FOR NOW!!!)
[d_omega_theta,tau,phi] = createCDMA(M,r,incidentAngle);
MicLocation = r * (cos(phi(:))+ 1j * sin(phi(:)));
plot(real(MicLocation), imag(MicLocation),'o')
xlim([(-r - 0.2*r),(r + 0.2*r)])
ylim([(-r - 0.2*r),(r + 0.2*r)])
axis equal
grid on

P1 = eye(M);
P2 = vertcat(zeros(1,M),eye(M)); P2 = horzcat(P2, vertcat(1,zeros(M,1)));
P3 = P2^2;
P4 = P2^3;
P5 = P2^4;
P6 = P2^5;
%define mic_m input signal
%define frequency weights H_m(omega)
    %define and calculate the parameters needed for H_m(omega)
%define beampattern B_N(theta - theta_s)
    %define coefficients a for beampattern equation
% After DOA use only 2 mics for 1st order LDMA or 3 for 2nd order LDMA
% For adaptivity, use LNMS algorithm and check DSPLab_project and Elko
    % paper for examples