% 利用MATLAB实现FMCW雷达的距离多普勒估计
% https://mp.weixin.qq.com/s/-ES3mstm2O_aLv0ljMOnJA
clear all;
clc;
close all;
%% Radar Specifications 
c = 3e8;
fc= 77e9;             %carrier freq                                                        
Nd=128;                
Nr=256; 
lamuda = c/fc;
%% User Defined Range and Velocity of target
range = 100;
vel = 1;
Tchirp = 40e-6;
range_res = 1;
max_vel = lamuda/(4*Tchirp);% m/s

%% FMCW Waveform Generation
B = c / (2*range_res);
vel_res = lamuda/(2*Tchirp*Nd);
slope = B/Tchirp;
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
r_t=zeros(1,length(t));

td = 2 * (range + vel*t)/ c;
Tx  = exp(1j*2*pi*(fc*t + (slope*t.*t)/2 ) );
Rx  = exp(1j*2*pi*(fc*(t -td) + (slope * (t-td).*(t-td))/2 ) );
Mix = Tx.*conj(Rx);%混频

Mix_fft1d =Mix+awgn(Mix,20);%增加高斯白噪声
Mix_fft1d = reshape(Mix_fft1d, [Nr, Nd]);

%%原始信号
figure(1);
subplot(2,1,1);
plot(abs((Mix_fft1d(:,1))));
xlabel('采样点数(N)')
ylabel('幅度(A)')
subplot(2,1,2);
plot(abs(fft(Mix_fft1d(:,1))));
xlabel('采样点数(N)')
ylabel('幅度(A)')

signal_fft = fft(Mix_fft1d, Nr);
figure(2)
mesh(db(abs(signal_fft./max(signal_fft))));
title('Range from First FFT');
ylabel('Range [m]');
xlabel('chirps[N]')
zlabel('Amplitude (dB)');
title('距离维FFT')
for i=1:Nr
    signal_fft2(i,:) = fftshift(fft(signal_fft(i,:)));
end
doppler_axis = vel_res*((-Nd/2:(Nd-1)/2)-1);
range_axis = range_res*((1:Nr)-1);

figure(3);
mesh(doppler_axis,range_axis,db(abs(signal_fft2)));
xlabel('Speed[m/s]');
ylabel('Range[m]');
zlabel('Amplitude[dB]');
title('Amplitude and Range From FFT2');

%% CFAR implementation
RDM =(abs(signal_fft2));
n_train_cells = 8;
n_train_bands = 8;
%test (CUT) for accurate estimation
n_guard_cells = 4;
n_guard_bands = 4;
offset = 1.4;
noise_level = zeros(1,1);
RDM = RDM / max(RDM(:));
for row0 = n_train_cells + n_guard_cells + 1 : Nr - (n_train_cells + n_guard_cells)
  for col0 = n_train_bands + n_guard_bands + 1 : (Nd) - (n_train_bands + n_guard_bands)
    noise_level = zeros(1, 1);
    for row1 = row0 - (n_train_cells + n_guard_cells) : row0 + (n_train_cells + n_guard_cells)
      for col1 = col0 - (n_train_bands + n_guard_bands) : col0 + (n_train_bands + n_guard_bands)
        if (abs(row0 - row1) > n_guard_cells || abs(col0 - col1) > n_guard_bands)
          noise_level = noise_level + db2pow(RDM(row1, col1));
        end
      end
    end
    % Calculate threshold from noise average then add the offset
    thresh = pow2db(noise_level / (2 * (n_train_bands + n_guard_bands + 1) * 2 * (n_train_cells + n_guard_cells + 1) - (n_guard_cells * n_guard_bands) - 1));
    thresh = thresh + offset;

    CUT = RDM(row1-(n_train_cells + n_guard_cells) ,col1-(n_train_bands + n_guard_bands));
    if (CUT < thresh)
      RDM(row0, col0) = 0;
else
      RDM(row0, col0) = 1;
    end

  end
end

RDM(RDM~=0 & RDM~=1) = 0;
figure('Name', 'CA-CFAR Filtered RDM')
mesh(doppler_axis,range_axis,RDM);
title( 'CA-CFAR Filtered RDM surface plot');
xlabel('Speed[m/s]');
ylabel('Range[m]');
zlabel('Amplitude[dB]');
%%$ ending