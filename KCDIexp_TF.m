%% KCDIģ��ʵ�飬TFģ��
clc
close all
clear
%% ��λȫ���� �ı�CCD�Ĳ����߳����ı���λ��
% L1=2867e-6;      %side length(m)  28672΢��=28.672mmʧ��
L1=28672e-6; 
M=2048;       %number of samples
dx1=L1/M;    %src sample interval
x1=-L1/2:dx1:L1/2-dx1;    %src coords
y1=x1;
lambda=632e-9;     %wavelength(m)
k=2*pi/lambda;      %wavenumber
z1=0.601;             %��͸����������������ľ���474+127mm
z2=0.463;             %��ccd�ķ������������463mm
zf=0.474;           %������۵ľ��루���ࣩ474mm
wl=400;              %͸����ͷ�ھ��뾶(����Ϊ��λ)
% wl=400;    %����ֵ͸���ھ�(����Ϊ��λ)
[X1,Y1]=meshgrid(x1,x1);
% u1=circ(sqrt(X1.^2+Y1.^2)/wl);   %Բ�β�
[x_array,y_array] = meshgrid(1:M,1:M); 
x_array = x_array - floor(max(x_array(:))/2+1); % center of image to be zero 
y_array = y_array - floor(max(y_array(:))/2+1); % center of image to be zero 
u1=(x_array./wl).^2+(y_array./wl).^2 <= 1; 
% figure;imagesc(u1);axis square;title('͸����С');

%% ͸������ⳡ
uout=u1.*exp(-1i*k/(2*zf)*(X1.^2+Y1.^2));
I1=(abs(u1)).^2;
% figure;imagesc(x1,y1,I1);colormap('gray');axis square;title('͸��ǰ')
% figure;imagesc(x1,y1,angle(u1));colormap('gray');axis square;title('ͨ��͸��ǰ����λ');

Iout=(abs(uout)).^2;
figure;imagesc(x1,y1,Iout);colormap('gray');title('ͨ��͸����');axis square
figure;imagesc(x1,y1,angle(uout));colormap('gray');axis square;title('ͨ��͸�������λ');

u2=propTF(u1,L1,lambda,z1); %propagation, TF�����u2�ı߳�L2=L1
I2=abs(u2.^2);             %obs irrad
x2=x1;
y2=y1;
% figure;
%imagesc(x2,y2,I2);axis xy;axis square;colormap('gray');xlabel('x(m)');ylabel('y(m)');title(['û��ͨ��͸��z=',num2str(z1),'m']);
%% ͨ��͸�������������
%͸�������ķ���������propagation, TF�����u2�ı߳�L2=L1
u3=propTF(uout,L1,lambda,z1); 
% u3=propIR(u1,L1,lambda,z1);   %ȥ��͸��������
I3=abs(u3.^2);             %obs irrad
x3=x1;
y3=y1;
% figure;
%imagesc(x2,y2,I3);axis xy;axis square;colormap('gray');xlabel('x(m)');ylabel('y(m)');title(['����������z=',num2str(z1),'m','������Ʒǰ']);
%% ������Ʒ
obj1=imread('USAF1951B250_2048','png');   %read image file
obj1=rgb2gray(obj1);
obj=obj1;
Iobj=single(obj);                     %integer to floating ����ת��Ϊ�����ȸ�����
Iobj=Iobj/max(max(Iobj));               %normalize ideal image
uobj=sqrt(Iobj);                      %ideal image field
xobj=x1;yobj=y1;
figure%ԭͼ
imagesc(xobj,yobj,Iobj);title('ԭͼ');colormap('gray');xlabel('u(m)');ylabel('v(m)');axis square;
uobj=circshift(uobj,-100,1);
%% ƽ�沨ͨ����Ʒ�������������Ʒ��͸����
u4=u3.*uobj;  
% u4 = u3;   %û�������Ʒ��ֻ����͸���ķ���������
x4=x1;
y4=y1;
I4=(abs(u4)).^2;
figure('color',[1 1 1]);imagesc(x4,y4,I4);colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');%title('ͨ����Ʒ��Ĺ�ǿ');
figure;imagesc(x4,y4,(angle(u4)),[-pi,pi]);axis square;title('ͨ����Ʒ����λ');
%% ���������䵽��ccd
L5=L1;
M=2048;       %number of samples
dx5=L5/M;    %src sample interval
x5=-L5/2:dx5:L5/2-dx5;    %src coords
y5=x5;
u5=propTF(u4,L5,lambda,z2); %ԭ����TF
% u5=propIR(u4,L5,lambda,z2); 
I5=(abs(u5)).^2;
% figure;imagesc(x5,y5,I5);title('CCD intensity');colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');
figure('color',[1 1 1]);imagesc(x5,y5,I5);colormap('gray');axis square;xlabel('x(m)');ylabel('y(m)');
figure('color',[1 1 1]);imagesc(x5,y5,angle(u5));axis square;xlabel('x(m)');ylabel('y(m)');%title('CCD phase');
% figure;imagesc(x5,y5,I5);title('CCD');axis square
% figure;imagesc(x5,y5,angle(u5));xlabel('x(m)');ylabel('y(m)');title('CCD phase');axis square
figure('color',[1 1 1]);plot(x5,unwrap(angle(u5(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;%title('CCD phase M/2+10');

%% ��������
% Intensity = I5;
% save('KCDI_CCDsample_intensity_TF.mat','Intensity')
% lens_Intensity = I5;
% save('KCDI_CCDlens_intensity_TF.mat','lens_Intensity')
% % Angle = angle(u5);
% % save('KCDI_CCDsample_angle_TF.mat','Angle')
% lens_Angle = angle(u5);
% save('KCDI_CCDlens_angle_TF.mat','lens_Angle')
% % % save('CCD_field.mat','u5')
%% ���嵥������
% uobj_dif=propTF(uobj,L5,lambda,z2);
% figure
% imagesc(x5,y5,abs(uobj_dif));
% xlabel('x(m)');ylabel('y(rad)');
% title('���嵥������')
% axis square;colormap('gray')
% Angle_uobj = angle(uobj_dif);
% save('KCDI_CCD_Angle_uobj.mat','Angle_uobj')
%% inverse IR
% u42=propIR_inverse(u5,L5,lambda,z2); 
% figure;imagesc(x5,y5,abs(u42));xlabel('x(m)');ylabel('y(m)');title('IR��任');axis square
% figure;imagesc(x5,y5,abs(u5));xlabel('x(m)');ylabel('y(m)');title('IR���任');axis square
%% kk��ϵ
ws = 100;
ulens_2=propTF(u3,L1,lambda,z2); 
kr=angle(ulens_2);
kk_imaginary_part = kk3(I5,L1,lambda,M);
initial_diffraction=exp((log(I5)./2-1i*kk_imaginary_part)+1i*kr);   %krд��͸����������������֮�����λ
ph_p_kk=imag(log(initial_diffraction));
figure('color',[1 1 1]);imagesc(x5,y5,angle(initial_diffraction));axis square;xlabel('x(m)');ylabel('y(m)');%title('CCD phase');
figure('color',[1 1 1]);plot(x5,unwrap(angle(initial_diffraction(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');axis square;%title('CCD phase M/2+10');
figure
subplot(1,2,1);imagesc(x1,y1,ph_p_kk);xlabel('x(m)');ylabel('y(m)');title('kk E(r) angle','fontsize',22);axis square
subplot(1,2,2);plot(x1,unwrap(phase(initial_diffraction(M/2+10,:))));xlabel('x(m)');ylabel('y(rad)');title('kk E(r) angle M/2+10','fontsize',22);axis square