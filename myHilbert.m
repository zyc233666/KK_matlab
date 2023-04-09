function [y]=myHilbert(xn) 
    [xn,nshifts] = shiftdim(xn);
    n = length(xn); %����x�ĳ���
    x = fft(xn,n,1);   %��������������FFT�����źŴ�ʱ��ת����Ƶ��
    h = zeros(n,~isempty(x));%����ϵͳ����h
    %%%ϵͳ����
    if n > 0 && 2*fix(n/2) == n
      % even and nonempty
      h([1 n/2+1]) = 1;
      h(2:n/2) = 2;
    elseif n>0
      % odd and nonempty
      h(1) = 1;
      h(2:(n+1)/2) = 2;
    end
    yy = ifft(x.*h(:,ones(1,size(x,2))),[],1);
    % Convert back to the original shape.
    y = shiftdim(yy,-nshifts);
end

