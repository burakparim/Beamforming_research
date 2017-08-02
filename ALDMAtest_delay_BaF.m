% Plot frequency response (=beamformer response evaluated for one specific angle over all possible frequencies),
% directional response (=beamformer response evaluated for one specific frequency
% over all possible directions), and beampattern (=beamformer response evaluated
% for all combinations of angles and frequencies)

clear all;
close all;
% clc;
%% ----------------------------------------------------------------------------------------
% Directional response of forward and backward cardioids with additional delay element T
%----------------------------------------------------------------------------------------
% Microphone spacing
d = 0.0215; % in m

% speed of sound
c = 343; % in m/s

% Define look directions 'theta' in ° or rad (depending on which trigonometric functions you use) for which the frequency response is evaluated
theta = linspace(-pi,pi,1000);

% Define fixed frequency value 'f' in Hz for which the frequency response is evaluated
% f = [50 200 500 1000 3000 4000]; %[Hz]
f = 1000;
% Angular frequency omega (depends on f)
omega = 2*pi*f;

% Angular frequency over speed of sound
k = omega/c;

% Delay element
T= d/c;

figure;
%Calculate directional response by evaluating (4) or (5) in G. W. Elko, A Simple Adaptive First-Order Differential Microphone
car_fw = 2*1j*1*exp(-1j*omega*(T/2))*sin(k*d*(1 + cos(theta))/2);polar(theta,abs(car_fw),'--r');hold on
car_bw = 2*1j*1*exp(-1j*omega*(T/2))*sin(k*d*(1 - cos(theta))/2);polar(theta,abs(car_bw),':');

%Plot directional response over angles. 

title('Directivity pattern of forward and backward cardioid')
ylabel('response in Magnitude')
legend('forward cardioid', 'backword cartdioid')

%% ----------------------------------------------------------------------------------------
% Frequency response of forward and backward cardioids without additional delay element T
%----------------------------------------------------------------------------------------
% Microphone spacing
d = 0.0215; % in m

% speed of sound
c = 343; % in m/s

% Define fixed look direction 'theta' in ° or rad (depending on which trigonometric functions you use) for which the frequency response is evaluated
theta = [pi (3*pi/2) pi/2 pi/4 0];
% theta = pi/2;

% Define frequency values 'f' in Hz for which the frequency response is evaluated
f = linspace(0,16000,1000);

% Angular frequency omega (depends on f)
omega = 2*pi*f;

% Delay element
T= d/c;

% Angular frequency over speed of sound
k = omega/c;

figure;
%Calculate frequency response by evaluating (4) or (5) in G. W. Elko, A Simple Adaptive First-Order Differential Microphone
%Plot frequency response over frequency. 
car_fw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k*d*(1 + cos(theta(1)))/2);plot(f,abs(car_fw),'LineWidth',2.5);hold on
car_bw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k*d*(1 - cos(theta(1)))/2);plot(f,abs(car_bw),'--r','LineWidth',2.5);
car_fw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k*d*(1 + cos(theta(2)))/2);plot(f,abs(car_fw),'m','LineWidth',2.5);
car_bw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k(1:end)*d*(1 - cos(theta(2)))/2);plot(f,abs(car_bw),'--c','LineWidth',2.5);
car_fw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k(1:end)*d*(1 + cos(theta(3)))/2);plot(f,abs(car_fw),'y','LineWidth',2.5);
car_bw = 2*1j*1*exp(-1j*omega*(T/2)).*sin(k(1:end)*d*(1 - cos(theta(3)))/2);plot(f,abs(car_bw),'--k','LineWidth',2.5);

grid on
title('Frequency response of the forward and backward cardioid')
xlabel('frequency axis in Hz')
ylabel('|Y(omega,theta)|')
legend('forward cardioid(theta = pi/2)', 'backward cardioid(theta = pi/2)','forward cardioid(theta = pi/4)', 'backward cardioid(theta = pi/4)','forward cardioid(theta = 0)', 'backward cardioid(theta = 0)')
%% ----------------------------------------------------------------------------------------
% Beampattern of forward and backward cardioids with additional delay element T
%----------------------------------------------------------------------------------------
% Microphone spacing
d = 0.0215; % in m

% speed of sound
c = 343; % in m/s

% Define look directions 'theta' in ° or rad (depending on which trigonometric functions you use) for which the frequency response is evaluated
theta = linspace(-pi,pi,1000);

% Define frequency values 'f' in Hz for which the frequency response is evaluated
f =  linspace(0,16000,1000);

% Angular frequency omega (depends on f)
omega = 2*pi*f;

% Delay element
T= d/c;

% Angular frequency over speed of sound
k = omega/c;

%Evaluate (4) or (5) in G. W. Elko, A Simple Adaptive First-Order
%Differential Microphone for all combinations of specified angles theta and
%all specified frequencis in f
Z_fw = zeros(length(omega),length(theta)); %init forward cardioid
Z_bw = zeros(length(omega),length(theta)); %init backward cardioid

for i = 1:length(omega)
    for ii = 1:length(theta)
Z_fw(i,ii) = 2*1j*exp(-1j*omega(i)*(T/2))*sin(k(i)*d*(1 + cos(theta(ii)))/2);
Z_bw(i,ii) = 2*1j*exp(-1j*omega(i)*(T/2))*sin(k(i)*d*(1 - cos(theta(ii)))/2);
    end
end
%Plot the beampattern in dB -> 20*log10(Z) over angles and frequency
%
figure;
imagesc(theta,f,20*log10(abs(Z_fw)));
title('Beam pattern of forward cardioid')
xlabel('Angle-axis')
ylabel('Frequency-axis')
colorbar; %switch colorbar on
caxis([-50 10]); %set limit of c-axis
figure;
imagesc(theta,f,20*log10(abs(Z_bw)));
title('Beam pattern of backward cardioid')
colorbar; %switch colorbar on
colormap('default')
caxis([-50 10]); %set limit of c-axis
xlabel('Angle-axis')
ylabel('Frequency-axis')
