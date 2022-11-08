
%% ���ߣ���Ƥ������
%% ���ںţ���Ƥ��������
%% ʱ��2022��04��
%https://mp.weixin.qq.com/s/8xOw8_YsFrv7BKnb1m-MZg
%%
clc;
close all;
clear all;

%% �״����
Tx_Number = 2;               %��������
Rx_Number = 4;               %��������
Range_Number = 128;          %���������ÿ��chirp 128���㣩
Doppler_Number = 128;         %������ͨ����
global Params;
Params.NChirp = Doppler_Number;               %1֡���ݵ�chirp����
Params.NChan =  Rx_Number;                    %RxAn��,ADCͨ����
Params.NSample = Range_Number;                %ÿ��chirp ADC������
Params.Fs = 2.5e6;                          %����Ƶ��
Params.c = 3.0e8;                     %����
Params.startFreq = 77e9;              %��ʼƵ�� 
Params.freqSlope = 60e12;             %chirp��б��
Params.bandwidth = 3.072e9;           %��ʵ����
Params.lambda=Params.c/Params.startFreq;    %�״��źŲ���
Params.Tc = 144e-6;                         %chirp����
global FFT2_mag;

%% �������
[X,Y] = meshgrid(Params.c*(0:Params.NSample-1)*Params.Fs/2/Params.freqSlope/Params.NSample, ...
    (-Params.NChirp/2:Params.NChirp/2 - 1)*Params.lambda/Params.Tc/Params.NChirp/2);   
       
adc_data =load('angle_15.mat');
Data_dec=(adc_data.prompt_1);  %��16����ת��Ϊ10����

%% ���ݶ�ȡ����֡����
Data_zuhe=zeros(1,Tx_Number*Rx_Number*Doppler_Number*Range_Number*2); %��������洢���ݵĿվ���
for i=1:1:Tx_Number*Rx_Number*Doppler_Number*Range_Number*2
    
    Data_zuhe(i) = Data_dec((i-1)*2+1)+Data_dec((i-1)*2+2)*256;%�����ֽ����һ�������ڶ����ֽڳ���256�൱������8λ��
    if(Data_zuhe(i)>32767)
        Data_zuhe(i) = Data_zuhe(i) - 65536;  %���Ʒ���
    end
end

%% �ַ�����
ADC_Data=zeros(Tx_Number,Doppler_Number,Rx_Number,Range_Number*2); %��������洢���ݵĿվ���
for t=1:1:Tx_Number
    for i=1:1:Doppler_Number
        for j=1:1:Rx_Number
            for k=1:1:Range_Number*2 %ʵ���鲿
                ADC_Data(t,i,j,k) = Data_zuhe(1,(((t-1)*Doppler_Number+(i-1))*Rx_Number+(j-1))*Range_Number*2+k);%ʱ����������˳��Ϊ TX1 TX2
            end
        end
    end
end

%% ��ӡȫ����ʵ������
Re_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %��������洢���ݵĿվ���
Im_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %��������洢���ݵĿվ���

% �鲿ʵ���ֽ�
for i=1:1:Tx_Number*Rx_Number*Doppler_Number*Range_Number
    Im_Data_All(1,i) = Data_zuhe(1,(i-1)*2+1);
    Re_Data_All(1,i) = Data_zuhe(1,(i-1)*2+2);
end

% ԭʼ�ź�ʵ�����鲿ͼ�λ��� 
% figure()
% subplot(2,1,1);
% plot(Im_Data_All(1,1:3000));title('ʵ������');
% xlabel('��������');
% ylabel('����');
% subplot(2,1,2);
% plot(Re_Data_All(1,1:3000),'r');title('�鲿����');
% xlabel('��������');
% ylabel('����');

%% ��ӡ������ʵ������ ���ݽṹΪ��2T4R��TX2���16����������
Re_Data=zeros(Doppler_Number,Range_Number); %��������洢���ݵĿվ���
Im_Data=zeros(Doppler_Number,Range_Number); %��������洢���ݵĿվ���

for chirp=1:Doppler_Number %�鿴����chirp������ 
for j=1:1:Tx_Number
    for k=1:1:Rx_Number
        for i=1:1:Range_Number
            Re_Data(chirp,i) = ADC_Data(j,chirp,k,(i-1)*2+2);
            Im_Data(chirp,i) = ADC_Data(j,chirp,k,(i-1)*2+1);
        end
    end 
end
end

%% �鲿+ʵ����������õ����ź� 
ReIm_Data = complex(Re_Data,Im_Data); %����ֻ���������ߵ����һ�����ݡ�ԭ�����ݴ�СӦ����16*256*8=32768������ֻ��16*256*1=4096��
ReIm_Data_All =complex(Re_Data_All,Im_Data_All);

ReIm_Data_all1 = zeros(Range_Number,Doppler_Number,4);
ReIm_Data_all2 = zeros(Range_Number,Doppler_Number,4);

%% ������������ 4ͨ��->8ͨ��
for nn=1:4
    for mm=1:Range_Number       
            ReIm_Data_all1(mm,:,nn) = ReIm_Data_All((nn-1)*Range_Number+ ((mm-1)*4*Range_Number+1):((mm-1)*4*Range_Number+Range_Number)+(nn-1)*Range_Number  );          
            ReIm_Data_all2(mm,:,nn) = ReIm_Data_All((nn-1)*Range_Number+131072/2+((mm-1)*4*Range_Number+1):131072/2+((mm-1)*4*Range_Number+Range_Number) +(nn-1)*Range_Number );
    end
end

ReIm_Data_All = cat(3,ReIm_Data_all1(:,:,1:4), ReIm_Data_all2(:,:,1:4));

%% 1D FFT
fft1d= zeros(Doppler_Number,Range_Number,8);
for qq =1:8
    for chirp_fft=1:Doppler_Number 
        fft1d(chirp_fft,:,qq) = fft((ReIm_Data_All(chirp_fft,:,qq)));
    end
end

FFT1_mag=abs(fft1d(:,:,1));
figure(1);
mesh(FFT1_mag);
xlabel('��������');ylabel('������');zlabel('����');
title('����άFFT���');

%%  MTI ��Ŀ����ʾ
fft1d_MTI= zeros(Range_Number,Doppler_Number,8);
for cc =1:8
    for ii =1:Doppler_Number-1
        fft1d_MTI (ii,:,cc) = fft1d(ii+1,:,cc)-fft1d(ii,:,cc);
    end
end
% 
% 
% mesh(abs(ReIm_Data_All(:,:,1)));

%%  ������ֵ�����㷨-��̬�Ӳ��˳�
fft1d_avg = zeros(128,128,8);
for n=1:8
    avg = sum(fft1d(:,:,n))/128;

    for chirp=1:128
        fft1d_avg(chirp,:,n) = fft1d(chirp,:,n)-avg;
    end
end

% figure;
% mesh(X,Y,abs(fft1d_jingtai(:,:,1)));
% fft1d =fft1d_jingtai;
%%

%% 2D FFT 

fft2d= zeros(Doppler_Number,Range_Number,8);
for kk=1:8
    for chirp_fft=1:Range_Number 
        fft2d(:,chirp_fft,kk)     =fftshift( fft((fft1d(:,chirp_fft,kk)))); %δ������̬�Ӳ��˳�
        fft2d_MTI(:,chirp_fft,kk) =fftshift( fft((fft1d_MTI(:,chirp_fft,kk)))); %MTI
        fft2d_avg(:,chirp_fft,kk) =fftshift( fft((fft1d_avg(:,chirp_fft,kk)))); %������ֵ��
    end
end

figure(2);
fft2d_0 =abs(fft2d(:,:,1));
mesh(X,Y,fft2d_0);
xlabel('����ά(m)');ylabel('�ٶ�ά(m/s)');zlabel('����');
title('R-V�׾���');

figure(3);
fft2d_0 =abs(fft2d(:,:,1));
fft2d_0(63:65,:)=0;
mesh(X,Y,fft2d_0);
xlabel('����ά(m)');ylabel('�ٶ�ά(m/s)');zlabel('����');
title('����ͨ��ֱ�����㷨');


figure(4);
mesh(X,Y,abs(fft2d_MTI(:,:,1)));
xlabel('����ά(m)');ylabel('�ٶ�ά(m/s)');zlabel('����');
title('MTI');

figure(5);
mesh(X,Y,abs(fft2d_avg(:,:,1)));
xlabel('����ά(m)');ylabel('�ٶ�ά(m/s)');zlabel('����');
title('������ֵ����');

%ֱ������ȥ�� ������ֵ�����㷨-��̬�Ӳ��˳�
figure(6)
fft2d_avg(:,1,1) =0;
mesh(X,Y,abs(fft2d_avg(:,:,1)));
xlabel('����ά(m)');ylabel('�ٶ�ά(m/s)');zlabel('����');
title('������ֵ����');

