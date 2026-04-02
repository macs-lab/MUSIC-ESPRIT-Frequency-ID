% Multiple narrow band frequency identification
% Xu Chen
% 2010-05-07

%% Data load
clear all
close all

% ======load the data here======
load('y1'); 
Ts = 1;
y1          = y;
datalength  = length(y1);

% if you are using this sample y, input p = 3 in the following.
while(1)
    p = input('Input the number of narrow band signals: ');
    if p >= 1
        break;
    else
        disp('Wrong selection, reselect please.');
    end
end

% The real part and the imaginary part (should both be sum of sinusoidals)
rey1        = real(y1); imy1 = imag(y1);
% Magnitude and phases of the complex exponentials
magy1       = abs(y1); phy1 = angle(y1);
%% Plot the data
% y1 can be either complex or real. 
% Real parts and imaginary parts
% y1
figure;
subplot(211);
plot(rey1);ylabel('real part');title('y1');
hold on;plot(mean(rey1)*ones(datalength,1),'r--')
subplot(212);plot(imy1);ylabel('imaginary part');
hold on;plot(mean(imy1)*ones(datalength,1),'r--')
xlabel('sample')

% Amplitudes and phases
% y1
figure;
subplot(211);
plot(magy1);ylabel('magnitude');title('y1');
hold on;plot(mean(magy1)*ones(100,1),'r--')
subplot(212);plot(phy1);ylabel('normalized phase');
hold on;plot(mean(phy1)*ones(100,1),'r--')
xlabel('sample')
%%
specCal(y1,1/Ts)
%%
disp('============MUSIC=============')
wregion     = 0:0.001:2*pi;
if p==1
    M = 6;
else
    M = p^2;
end
Pmusic1     = m_music(y1,p,M,wregion);
[Pmusic1_sorted,ind1]       = sort(Pmusic1);
figure;
plot(wregion,Pmusic1),xlabel('Frequency (rad)')
% plot(wregion/pi,Pmusic1),xlabel('Frequency (\pi)')
grid,title('MUSIC Spectrum'),ylabel('Gain (dB)')
hold on
[ymax,imax,ymin,imin]       = extrema(Pmusic1);
plot(wregion(imax(1:p)),ymax(1:p),'g.')
%% ESPRIT
disp('============ESPRIT=============')
w1          = m_esprit(y1,p,M);
%% AF with Cadzow's denoising method
disp('============AF with denoising=============')
[w1_AF,sigma_AF]            = AF_with_denoising(y1,p,M);
figure;semilogy(sigma_AF,'*-');
xlabel('Index on the diagonal');ylabel('Singular values \sigma')
%% Result Comparing
disp('Peak frequenceis from MUSIC (in rad):')
sort(wregion(imax(1:p)))'
disp('Calculated frequencies from ESPRIT (in rad):')
sort(w1)
disp('Calculated frequencies from annihilating filter method with Catzow denoising:')
sort(w1_AF)