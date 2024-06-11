clc;
clear;
th=-90:.1:90;
t=(sind(th));
lambda=300/10;
d=15.4;

k=2*pi/lambda;
AF=0;
scan=0;
trunc=-60;
N=20;

% map=[];
% for j=1:1000
%     amp = rand(1, 20);
%     map(j)=amp;
% end
% 
amp=ones(1,N);
%  amp=taylorTappOdd((N-1)/2,35,5);  % for N odd
%amp=taylorTappEven(N/2,30,5);  %% for N even 
% mag=[];
% for j=1:1000
%     ph = 0:5.625:360;
%     rm = randi(length(ph), 1, 20);
%     phase = ph(rm);
%     mag(j)=phase;
% end


ph=zeros(1,N);
for k1=1:N
    AF=AF+amp(k1)*exp(-j*ph(k1))*exp(-j*k*d*(k1)*(t-sind(scan)));
end
    AFn=(abs(AF)./max(abs(AF)));%.*cos(th*pi/180);
     AFdB=20*log10(abs(AFn));
    indices=find(AFdB<trunc);
    AFdB(indices)=trunc;
plot(th,AFdB,'LineWidth',1.5,'Color',[.1 .1 0.1]);
hold on 
grid on
% mesh(th,pi,AFdB)
beamwidth=-th(min(find(AFdB>=-3)))+th(max(find(AFdB>=-3)))


%Peak Side Lobe Level
pks=findpeaks(AFn);
%pks1=findpeaks(AFdB);

d=sort(pks,'descend');
value=d(1);
value1=d(2);
log=20*log10(value);
log1=20*log10(value1);
PSLL=log+log1

%Average of peaks
d(1)=d(1)-d(1);

b=numel(d);
sum=0;
for i=1:b
    sum=sum+d(i);
end

avg=sum/(b-1);
Average=20*log10(avg)

d(1)=[];
rms_sll = sqrt(mean(d .^ 2));
RMS = 20 * log10(rms_sll)



%RMS side lobe level
% x=0;
% for j=1:b
%     x=x^2+d(i)^2;
% end
% 
% rms=sqrt(x/(b-1));
% RMS=20*log10(rms)
% Directivity = 10*log10(19000/beamwidth)