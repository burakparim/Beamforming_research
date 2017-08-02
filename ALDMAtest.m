%% ----------------------------------------------------------------------------------------
% Directional response of combined forward-backward cardioid without additional delay element T
%----------------------------------------------------------------------------------------
% Microphone spacing
d = 0.0215; % in m

% speed of sound
c = 343; % in m/s

% Define look directions 'theta' in ° or rad (depending on which trigonometric functions you use) for which the frequency response is evaluated
theta = linspace(-pi,pi,1000);

% Define fixed frequency value 'f' in Hz for which the frequency response is evaluated
f = [50 200 500 1000 3000 4000]; %[Hz]
% f = 1000;
% Angular frequency omega (depends on f)
omega = 2*pi*f;

% Angular frequency over speed of sound
k = omega/c;
% Delay element
T= d/c;

% Weight parameter
beta = [.01 .1 .3 .5 .7 1];

%Calculate directional response by evaluating (6) in G. W. Elko, A Simple Adaptive First-Order Differential Microphone
%Plot directional response over angles. 

% Plots are showing the directivity pattern with different beta values on
% same figure
figure;
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb));
hold on
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(2)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'--r');
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(3)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),':c');
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(4)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'-.m');
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(5)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'r');
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(6)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'--k');
title('Directivity pattern')
ylabel('response in Magnitude')
legend('beta = 0.01','beta = 0.1','beta = 0.3','beta = 0.5','beta = 0.7','beta = 1.0')

figure;
car_comb = 2*abs(sin(k(6)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(6)*d*(1 - cos(theta))/2)));polar(theta,(car_comb));
hold on
car_comb = 2*abs(sin(k(5)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(5)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'--r');
car_comb = 2*abs(sin(k(4)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(4)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),':c');
car_comb = 2*abs(sin(k(3)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(3)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'-.m');
car_comb = 2*abs(sin(k(2)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(2)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'r');
car_comb = 2*abs(sin(k(1)*d*(1 + cos(theta))/2) - (beta(1)*sin(k(1)*d*(1 - cos(theta))/2)));polar(theta,(car_comb),'--k');
title('Directivity pattern')
ylabel('response in Magnitude')
legend('Freq. = 4000Hz','Freq. = 3000Hz','Freq. = 1000Hz','Freq. = 500Hz','Freq. = 200Hz','Freq. = 50Hz')
%% ----------------------------------------------------------------------------------------
% Frequency response of combined forward-backward cardioid without additional delay element T
%----------------------------------------------------------------------------------------
% Microphone spacing
d = 0.0215; % in m
dd = [0.0154 .0215 .0300 .0350 .0428 .0504];
% speed of sound
c = 343; % in m/s

% Define fixed look direction 'theta' in ° or rad (depending on which trigonometric functions you use) for which the frequency response is evaluated
theta = [pi/2, pi/4, 0];

% Define frequency values 'f' in Hz for which the frequency response is evaluated
f = linspace(0,16000,1000);

% Angular frequency omega (depends on f)
omega = 2*pi*f;

% Angular frequency over speed of sound
k = omega/c;

% Weight parameter
beta = [.01 .1 .3 .5 .7 1.0];

%Calculation of frequency response by evaluating (6) in G. W. Elko, A Simple Adaptive First-Order Differential Microphone


figure;
subplot(311)
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(1)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),'LineWidth',2.5);
hold on
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(2)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),'--r','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(3)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),':g','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(4)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),'-.m','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(5)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),'c','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(1)))/2) - (beta(6)*sin(k*d*(1 - cos(theta(1)))/2)));plot(f,abs(car_comb),'--k','LineWidth',2.5);
title('Frequency response(theta=pi/2)')
xlabel('frequency axis in Hz')
ylabel('|Y(omega,theta)|')
legend('beta = 0.01','beta = 0.1','beta = 0.3','beta = 0.5','beta = 0.7','beta = 1.0')

subplot(312)
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(1)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'LineWidth',2.5);
hold on
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(2)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'--r','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(3)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),':g','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(4)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'-.m','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(5)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'c','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(6)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'--k','LineWidth',2.5);
title('Frequency response(theta=pi/4)')
xlabel('frequency axis in Hz')
ylabel('|Y(omega,theta)|')
legend('beta = 0.01','beta = 0.1','beta = 0.3','beta = 0.5','beta = 0.7','beta = 1.0')

subplot(313)
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(1)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'LineWidth',2.5);
hold on
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(2)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'--r','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(3)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),':g','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(4)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'-.m','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(5)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'c','LineWidth',2.5);
car_comb = 2*abs(sin(k*d*(1 + cos(theta(2)))/2) - (beta(6)*sin(k*d*(1 - cos(theta(2)))/2)));plot(f,abs(car_comb),'--k','LineWidth',2.5);
title('Frequency response(theta=0)')
xlabel('frequency axis in Hz')
ylabel('|Y(omega,theta)|')
legend('beta = 0.01','beta = 0.1','beta = 0.3','beta = 0.5','beta = 0.7','beta = 1.0')
figure;
car_comb = 2*abs(sin(k*dd(1)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(1)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),'LineWidth',2.5);
hold on
grid on
car_comb = 2*abs(sin(k*dd(2)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(2)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),'--r','LineWidth',2.5);
car_comb = 2*abs(sin(k*dd(3)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(3)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),':g','LineWidth',2.5);
car_comb = 2*abs(sin(k*dd(4)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(4)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),'-.m','LineWidth',2.5);
car_comb = 2*abs(sin(k*dd(5)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(5)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),'c','LineWidth',2.5);
car_comb = 2*abs(sin(k*dd(6)*(1 + cos(theta(3)))/2) - (beta(1)*sin(k*dd(6)*(1 - cos(theta(3)))/2)));plot(f,20*log10((car_comb)),'--k','LineWidth',2.5);
title('Frequency response(theta=0)')
xlabel('frequency axis in Hz')
ylabel('20*log(Y(omega,theta)) [dB]')
legend('Int.Dist. = 0.0154m','Int.Dist. = 0.0215m','Int.Dist. = 0.03m','Int.Dist. = 0.035m','Int.Dist. = 0.0428m','Int.Dist. = 0.0504m')
%% ----------------------------------------------------------------------------------------
% Beampattern of combined forward-backward cardioid without additional delay element T
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
% T= d/c;
T = 0;
% Angular frequency over speed of sound
k = omega/c;

% Weight parameter
beta =0.1;

%Evaluate (6) in G. W. Elko, A Simple Adaptive First-Order
%Differential Microphone for all combinations of specified angles theta and
%all specified frequencis in f
Z = zeros(length(omega),length(theta));
for i = 1:length(omega)
    for ii = 1:length(theta)
Z(i,ii) = 2*abs(sin(k(i)*d*(1 + cos(theta(ii)))/2) - (beta*sin(k(i)*d*(1 - cos(theta(ii)))/2)));
    end
end
%Plot the beampattern in dB -> 20*log10(Z) over angles and frequency

figure;
imagesc(theta,f,20*log10(abs(Z)));
colorbar %switch colorbar on
colormap('default')
caxis([-50 10]); %set limit of c-axis
title('Beampattern of Combination of Forward & Backward Cardioid')
xlabel('Angle-axis[1/pi]')
ylabel('Frequency-axis[Hz]')