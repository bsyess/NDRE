function [sigout]=VTB_DRE(signalin,h,pnob)
%% �������ġ�Low-Resolution Digital Pre-Compensation Enabled by Digital Resolution Enhancer���߼���д
%%   author��cc
%% PNOBĬ�ϴ���1
if size(signalin,1) > 1
    signalin = signalin.';
end

sigini=real(signalin);
% Qbit_num=floor(2^pnob);
Qbit_num=floor(2^pnob);
I_max = max(sigini);
I_min = min(sigini);

step = (I_max-I_min)/Qbit_num;
for ii = 1:Qbit_num-1
    ipartition(1,ii) = I_min+step*ii;
end
for jj = 1:Qbit_num
    icodebook(1,jj)  = I_min + step*(jj-1);
end
[Iorder,iroundsig]=quantiz(sigini,ipartition,icodebook);

Iorder = Iorder+1;

% figure;hist(Iorder,2000)
% scatterplot(iroundsig+1i*qroundsig,2)
% figure;plot(iroundsig-sigini)
% figure;hist(iroundsig-sigini,1000)
I_error = iroundsig-sigini;

%% ...

n=length(icodebook);
m=length(h);
nvtb=length(signalin);
I_error = [I_error,I_error(1:m)];
Iorder = [Iorder,Iorder(1:m)];
Qro1 = zeros(3,m);
Qsq1 = zeros(3,m);
Qro2 = zeros(2,m);
Qsq2 = zeros(2,m);
for ii = m+1:nvtb+m;
    if Iorder(ii) ~=1 && Iorder(ii) ~=Qbit_num  % �ж��Ƿ��Ǳ�֡����Ϊ��֡�Ļ�ֻ�����ֿ�����
        Qro1(1,:) = I_error(ii-m+1:ii)-[zeros(1,m-1),step];
        Qro1(2,:) = I_error(ii-m+1:ii);
        Qro1(3,:) = I_error(ii-m+1:ii)+[zeros(1,m-1),step];
        Qsq1(1,:) = conv(Qro1(1,:),h,'same');
        Qsq1(2,:) = conv(Qro1(2,:),h,'same');
        Qsq1(3,:) = conv(Qro1(3,:),h,'same');
        [~,num] = min(sum(abs(Qsq1),2));
        I_error(ii) = Qro1(num,m);
    else
        Qro2(1,:) = I_error(ii-m+1:ii)-[zeros(1,m-1),step].*sign(Iorder(ii)-n/2); %% ����if���ж���䣬�ӿ���������ٶ�
        % ��������λ��������ͺʹε�λ����������λ��������ߺʹθ�λ
        Qro2(2,:) = I_error(ii-m+1:ii);
        Qsq2(1,:) = conv(Qro2(1,:),h,'same');
        Qsq2(2,:) = conv(Qro2(2,:),h,'same');
        [~,num] = min(sum(abs(Qsq2),2));
        I_error(ii) = Qro2(num,m);
    end
end

I_error1 = [I_error(end-m+1:end),I_error(m+1:nvtb)];

sigini1 = I_error1+sigini;

sigout = sigini1;
% scatterplot(sigout)


if size(sigout,2) > 1
    sigout = sigout.';
end

end
