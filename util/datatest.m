%This program reads two sets of .ASP files, for example RHS and LHS time series
%Outputs are: optical pathlegth, and absolute [Hb] and [HbO2] changes (determined assuming 80% water in optical path), 
%fractional changes in water content (from 840 nm band), changes in offset (scattering). 
%The spectra are input as two sets of .ASP files that are numbered sequentially, as provided by GRAMS.
%The function 'match' is required.
%
%Changes in [Hb] and [HbO2] are evaluated relative to the average of a set of baseline spectra that 
%are specified in input.
%The wavelength range to fit component spectra to (and from whicn concentration changes are derived) is specified
%on input
%
% This version also fits by default the second derivative of the 810-870 nm range with three components (Hb, HbO2, water), with
% the specific aim of getting the intensity of the 840 nm water band. Assuming the tissue to be 80% water, the
% optical pathlength at 840 nm is then determined as follows:
%
% A(810-870 range)=a1 * epsilon(H2O)  +  a2 * epsilon(Hb)  +  a3 * epsilon(HbO2) + a4 (offset term)
% 
% but a1 is the product (concentration*pathlength) for water. We know the concentration is 44.4 molar (80% of the
% concentration of pure water), so we may divide a1 (=(concentration*pathlength)), as provided by the least squares
% fit, by the known water concentration of 44.4 M to yield the estimated optical pathlength at 840 nm. 
%
% The optical pathlength at 840 nm is used to recover absolute changes in [Hb] and [HbO2] and water from their fitting
% coefficients (=concentration*pathlength). 


fclose('all');
clear *

nspectra=input ('How many spectra (.ASP files) to read? ');

%Here we are inputting files as split from a GRAMS multifile.
lhaspname=input ('Base name for LHS input files: ','s');
rhaspname=input ('Base name for RHS input files: ','s');


%Changes in Hb and HbO2 are evaluated relative to baseline. The spectra that constitute the baseline
%are input here.
firstbaseline=input ('First in range of baseline spectra: ');
lastbaseline=input ('Last in range of baseline spectra: ');

%To create a time axis, input the interval between spectra.
%Note: The last of the 'baseline' spectra is taken to fall at t=0.

deltatime=input ('Time interval between spectra: ');

%Input wavelength range to include in fit

leftwavelength=input('Input low wavelength in range to be fit: ');
rightwavelength=input('Input high wavelength in range to be fit: ');

% Make time axis. Note that t=0 for the last of the 'baseline' spectra
for i=1:nspectra
   time(i)=deltatime * (i-lastbaseline);
end
   
%read .ASP spectra
for i=1:nspectra
   x=int2str(i);    
   y=(strcat(lhaspname,x,'.ASP'));
   xx=load(y);
   lhspectra(:,i)=xx;
   clear y;
   y=(strcat(rhaspname,x,'.ASP'));
   xx=load(y);
   rhspectra(:,i)=xx;
   clear y;
end

% Make LHS wavelength axis. Note, there are npts-1 axis segments between firstx and lastx. 
load(strcat(lhaspname,'1.ASP'));
xx=load(strcat(lhaspname,'1.ASP'));
lhnpts=xx(1);
firstx=xx(2);
lastx=xx(3);
clear xx;
range=lastx-firstx;
for i=1:lhnpts
   lhaxis(i)=firstx+((i-1)/(lhnpts-1))*range;
end
% Make RHS wavelength axis. Note, there are npts-1 axis segments between firstx and lastx. 
load(strcat(rhaspname,'1.ASP'));
xx=load(strcat(rhaspname,'1.ASP'));
rhnpts=xx(1);
firstx=xx(2);
lastx=xx(3);
clear xx;
range=lastx-firstx;
for i=1:rhnpts
   rhaxis(i)=firstx+((i-1)/(rhnpts-1))*range;
end

%read absorptivity spectra and create arrays for slope and baseline.
hbb=load('hbepsilon.prn');
hbO2=load('hbO2epsilon.prn');
water=load('waterinterpepsilon.prn');


%Zero matrix of component spectra.
component=zeros(4,1100);
%Put absorptivities into single dimension arrays. Note: Spacing must be exactly 1 nm.
lngth=length(hbb);
for i=1:lngth
   lambda=hbb(i,1);
   hbepsilon(lambda)=hbb(i,2);
   component(1,lambda)=hbb(i,2);
end

lngth=length(hbO2);
for i=1:lngth
   lambda=hbO2(i,1);
   hbO2epsilon(lambda)=hbO2(i,2);
   component(2,lambda)=hbO2(i,2);
end

lngth=length(water);
for i=1:lngth
   lambda=water(i,1);
   waterepsilon(lambda)=water(i,2);
   component(3,lambda)=water(i,2);
end

%The fourth component spectrumn is a baseline (an array of ones).
component(4,:)=ones(length(component),1)';

% Get rid of .ASP header lines
lhspectra=lhspectra(7:lhnpts+6,:);
rhspectra=rhspectra(7:rhnpts+6,:);

%Evaluate mean of baseline spectra
nbaseline=lastbaseline-firstbaseline+1;
lhbaselinespec=0;
rhbaselinespec=0;
for i=1:nbaseline
   lhbaselinespec=lhbaselinespec+lhspectra(:,firstbaseline+i-1);
   rhbaselinespec=rhbaselinespec+rhspectra(:,firstbaseline+i-1);
end
lhbaselinespec=lhbaselinespec/nbaseline;
rhbaselinespec=rhbaselinespec/nbaseline;

%Subtract baseline spectrum from all others
for i=1:nspectra
   lhdeltaspectra(:,i)=lhspectra(:,i)-lhbaselinespec;
   rhdeltaspectra(:,i)=rhspectra(:,i)-rhbaselinespec;
end

%Convert spectra to round axis. Note that two new arrays are created, one for the rounded spectra (roundspectra), 
%and one for the rounded difference spectra (rounddeltaspectra).
lhorigroundaxis=round(lhaxis);
lhfirstx=lhorigroundaxis(1);
lhlastx=lhorigroundaxis(length(lhorigroundaxis));
lhroundaxis=[lhfirstx:1:lhlastx];

rhorigroundaxis=round(rhaxis);
rhfirstx=rhorigroundaxis(1);
rhlastx=rhorigroundaxis(length(rhorigroundaxis));
rhroundaxis=[rhfirstx:1:rhlastx];


%average datapoints with the same rounded wavelength
lhrounddeltaspectra (lhroundaxis, :) = interp1 (lhaxis, lhdeltaspectra, lhroundaxis, 'linear');
lhroundspectra (lhroundaxis, :) = interp1 (lhaxis, lhspectra, lhroundaxis, 'linear');
rhrounddeltaspectra (rhroundaxis, :) = interp1 (rhaxis, rhdeltaspectra, rhroundaxis, 'linear');
rhroundspectra (rhroundaxis, :) = interp1 (rhaxis, rhspectra, rhroundaxis, 'linear');


% lhrounddeltaspectra=zeros(lhlastx,nspectra);
% lhroundspectra=zeros(lhlastx,nspectra);
% lhnpts = length (lhroundaxis); %lhnpts=length(data1(1,1,:));
% rhnpts = length (rhroundaxis); %rhnpts=length(data2(1,1,:));
% kk=0;
% k=0;
% for i=lhfirstx:lhlastx
%    kk=kk+k;
%    k=0;
%    for j=kk+1:lhnpts
%       if(lhorigroundaxis(j)==i)
%          lhrounddeltaspectra(i,:)=lhrounddeltaspectra(i,:) + lhdeltaspectra(j,:);
%          lhroundspectra(i,:)=lhroundspectra(i,:) + lhspectra(j,:);
%       	k=k+1;
%       end
%    end
%    lhrounddeltaspectra(i,:)=lhrounddeltaspectra(i,:)/k;
%    lhroundspectra(i,:)=lhroundspectra(i,:)/k;
% end
% 
% rhrounddeltaspectra=zeros(rhlastx,nspectra);
% rhroundspectra=zeros(rhlastx,nspectra);
% kk=0;
% k=0;
% for i=rhfirstx:rhlastx
%    kk=kk+k;
%    k=0;
%    for j=kk+1:rhnpts
%       if(rhorigroundaxis(j)==i)
%          rhrounddeltaspectra(i,:)=rhrounddeltaspectra(i,:) + rhdeltaspectra(j,:);
%          rhroundspectra(i,:)=rhroundspectra(i,:) + rhspectra(j,:);
%       	k=k+1;
%       end
%    end
%    rhrounddeltaspectra(i,:)=rhrounddeltaspectra(i,:)/k;
%    rhroundspectra(i,:)=rhroundspectra(i,:)/k;
% end


%Tranpose so that each row is a spectrum.
lhrounddeltaspectra=lhrounddeltaspectra';
lhroundspectra=lhroundspectra';
rhrounddeltaspectra=rhrounddeltaspectra';
rhroundspectra=rhroundspectra';


%Make matrices of component spectra and difference spectra that contain only the wavelength range of interest.
%Make matrix of component spectra in the 810-870 nm region that is used to get the pathlength from water band intensity
componentseg=component(:,leftwavelength:rightwavelength);
componentwaterfitseg=component(:,810:870);
lhrounddeltaspectraseg=lhrounddeltaspectra(:,leftwavelength:rightwavelength);
lhroundspectraseg=lhroundspectra(:,leftwavelength:rightwavelength);
rhrounddeltaspectraseg=rhrounddeltaspectra(:,leftwavelength:rightwavelength);
rhroundspectraseg=rhroundspectra(:,leftwavelength:rightwavelength);
lhroundspectra810_870=lhroundspectra(:,810:870);
rhroundspectra810_870=rhroundspectra(:,810:870);

%Take derivatives for water fit.
derlhroundspectra=savgol(lhroundspectra,45,2,2);
derrhroundspectra=savgol(rhroundspectra,45,2,2);
dercomponent=savgol(component,45,2,2);
dercomponentwaterfitseg=dercomponent(:,810:870);
derlhroundspectra810_870=derlhroundspectra(:,810:870);
derrhroundspectra810_870=derrhroundspectra(:,810:870);

%Do the least squares. A=SC'((CC')^-1) where A=coefficients of Hb, HbO2, water and baseline component spectra,
%S is the matrix of difference spectra S=rounddeltaspectraseg and C the component spectra for the range of interest
%C=componentseg.

%'diffcoeffs' are the coefficients to fit the difference spectra, and 'speccoeffs' are the coefficients to fit the
%original spectra.

%First, do the 810-870 fit and get the pathlength:
% lhwaterspeccoeffs=lhroundspectra810_870 * componentwaterfitseg' * (inv(componentwaterfitseg*componentwaterfitseg'));
% rhwaterspeccoeffs=rhroundspectra810_870 * componentwaterfitseg' * (inv(componentwaterfitseg*componentwaterfitseg'));
% lhwaterbaselinecoeff=mean(lhwaterspeccoeffs(firstbaseline:lastbaseline,3));
% rhwaterbaselinecoeff=mean(rhwaterspeccoeffs(firstbaseline:lastbaseline,3));
% lhpathlength=lhwaterbaselinecoeff/44.4
% rhpathlength=rhwaterbaselinecoeff/44.4

%Do the pathlength evaluation using the derivative spectra:
%Note: Must eliminate the 'scatter' (baseline shift) component, which is zero in the derivative. 
dercomponentwaterfitseg=dercomponentwaterfitseg(1:3,:);
derlhwaterspeccoeffs=derlhroundspectra810_870 * dercomponentwaterfitseg' * (inv(dercomponentwaterfitseg*dercomponentwaterfitseg'));
derrhwaterspeccoeffs=derrhroundspectra810_870 * dercomponentwaterfitseg' * (inv(dercomponentwaterfitseg*dercomponentwaterfitseg'));
derlhwaterbaselinecoeff=mean(derlhwaterspeccoeffs(firstbaseline:lastbaseline,3));
derrhwaterbaselinecoeff=mean(derrhwaterspeccoeffs(firstbaseline:lastbaseline,3));
derlhpathlength=derlhwaterbaselinecoeff/44.4
derrhpathlength=derrhwaterbaselinecoeff/44.4

%Get changes in water relative to zero baseline. Note: 810-870 derivatives used for water ALWAYS
deltaderlhwaterspeccoeffs=derlhwaterspeccoeffs(:,3)-derlhwaterbaselinecoeff;
deltaderrhwaterspeccoeffs=derrhwaterspeccoeffs(:,3)-derrhwaterbaselinecoeff;
deltaderlhwaterspeccoeffs=deltaderlhwaterspeccoeffs/derlhpathlength;
deltaderrhwaterspeccoeffs=deltaderrhwaterspeccoeffs/derrhpathlength;

%Now, do the fit to the chosen spectral region:
lhdiffcoeffs = lhrounddeltaspectraseg * componentseg' * (inv(componentseg*componentseg'));
rhdiffcoeffs = rhrounddeltaspectraseg * componentseg' * (inv(componentseg*componentseg'));
lhspeccoeffs = lhroundspectraseg * componentseg' * (inv(componentseg*componentseg'));
rhspeccoeffs = rhroundspectraseg * componentseg' * (inv(componentseg*componentseg'));
lhdiffcoeffs = lhdiffcoeffs/derlhpathlength;
rhdiffcoeffs = rhdiffcoeffs/derrhpathlength;
lhspeccoeffs = lhspeccoeffs/derlhpathlength;
rhspeccoeffs = rhspeccoeffs/derrhpathlength;


%plot the spectra
figure
subplot(3,1,1), plot(time,lhdiffcoeffs(:,1),'r',time,lhdiffcoeffs(:,2),'r',time,(lhdiffcoeffs(:,1)+lhdiffcoeffs(:,2)),'r', time,rhdiffcoeffs(:,1),'b',time,rhdiffcoeffs(:,2),'b',time,(rhdiffcoeffs(:,1)+rhdiffcoeffs(:,2)),'b');
title('LHS is red, RHS is blue')
ylabel('Hb and HbO2. mM')
subplot(3,1,2), plot(time,deltaderlhwaterspeccoeffs,'r',time,deltaderrhwaterspeccoeffs,'b');
ylabel('Water. M')
subplot(3,1,3), plot(time,lhdiffcoeffs(:,4),'r',time,rhdiffcoeffs(:,4),'b');
ylabel('Scatter. H2O scaled')


%Make "best fit" spectra.
lhfitspectra=lhdiffcoeffs*component;
rhfitspectra=rhdiffcoeffs*component;

fclose('all');
