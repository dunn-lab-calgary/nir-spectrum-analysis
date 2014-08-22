function spc111 = spcfil(spc000,nn);

%construct the window, spc00 is the input data
%rising period nn: first nn1 data point, dropping period:the last nn2 data point
nn1=nn;
nn2=nn;
ww0=(nn1:-1:0).*pi/nn1;
aa1=0.5+0.5.*cos(ww0); %raised cosine function
ww1=(0:nn2).*pi/nn2;
aa2=0.5+0.5.*cos(ww1);
nn0=length(spc000);
xw=ones(1,nn0);
xw(1:(nn1+1))=aa1;
xw((nn0-nn2):end)=aa2;
spc001=spc000.*xw;


%Create a parks-mcclellan filter with passband ripple of 0.1db at 0.3pi,
%and stopband attenuation of 80db at 0.4pi.
As = 10^(-80/20); 
Ap = 1 - 10^(-0.1/20); 
[N,Fo,Ao,W] = firpmord ([0.3 0.4], [1 0], [Ap As]);

%adjust the order of filter N to meet the specification
N = N+4;
b = firpm (N, Fo, Ao, W);

%filter the spectrum
spc011=filter(b,1,spc001);
spc111=zeros(1,nn0);
spc111((nn+1):(end-N/2))= spc011((nn+1+N/2):end); % remove delay


