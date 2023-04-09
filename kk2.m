function [kk_imaginary_part] = kk2(I,L,lambda,M)
dx1=L/M;    %src sample interval
x1=-L/2:dx1:L/2-dx1;    %src coords
y1=x1;
k=2*pi/lambda;      %wavenumber
% [X1,Y1]=meshgrid(x1,y1);
%% ��CCD̽�⵽��ǿ��ͼ��ʹ��kk��ϵ����CCD��ⳡ��λ
%�Թ�ǿ��kk����ʹ��ϣ�����ر任
real_part=log(I)./2;
kk_imaginary_part = imag(hilbert(real(real_part).')).';
figure
imagesc(x1,y1,kk_imaginary_part)
title('�Թ�ǿ������kk�任')
xlabel('x(m)');ylabel('y(m)')
axis square
end