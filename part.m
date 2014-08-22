function varargout = part(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @part_OpeningFcn, ...
                   'gui_OutputFcn',  @part_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% -----------------------------------------------------------------------

% --- Executes just before part is made visible.
function part_OpeningFcn(hObject, eventdata, handles, varargin)

global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

addpath(genpath(progpath),'-begin');

handles.img1 = imread([progpath '/img/logo.jpg']); % Read the image file logo.jpg
handles.img2 = imread([progpath '/img/spectrum2.tif']);
set(hObject, 'Position', [1 10 865 620]);  % control the show size 

axes(handles.axes1);
image(handles.img1);
set(handles.axes1,'visible', 'off')

axes(handles.axes2);
image(handles.img2);
set(handles.axes2,'visible', 'off')

set(handles.edit1,'String',[  ]);

load ([progpath '/ref/calfit']);
handles.calfit=calfit;

load([progpath '/ref/RefLight.mat']);
handles.Ir0=Ir0;

handles.Isys=0;
handles.Imea=0;

% Choose default command line output for part
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes part wait for user response (see UIRESUME)

% --- Outputs from this function are returned to the command line.
function varargout = part_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function edit1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==========================================================================
% --- Executes on button press in initializeButton.
% --- INITIALIZE BUTTON ----
function initializeButton_Callback(hObject, eventdata, handles)

% initialize
% load error code table
global ANDOR;% make the variable ANDOR can be used by multi program
global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

load([progpath '/ref/andor_errors.mat']);
ANDOR = andor_errno; 

% path to Andor files
addpath(fullfile(matlabroot,'toolbox','Andor'));
path = '';
installpath = fullfile(matlabroot,'toolbox','Andor','Camera Files');
cd (installpath);

i1=int8(1);
i2=int8(1);
i3=int8(1);

%================================
% Initialize system, call the function 'Initialize'
result = AndorInitialize(path);

if (result ~= ANDOR. DRV_SUCCESS)
    set(handles.edit1,'String',['Error during initialization of Andor system: ' result]);
    i1=0;
end
   
if i1==1
    Isys=1;
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb '--System is initialized ']);
else
    Isys=0;
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb '--System is not successfully initialized ']);
end

%============================
% Initializes the measurement
i4=1;
i5=1;
i6=1;
i7=1;
i8=1;
i9=1;
i10=1;
i11=1;
i12=1;
i13=1;
i14=1;
i15=1;

bb = get(handles.edit1,'string');
set(handles.edit1,'String',[bb '--Measurement is initializing...       ']);

acquisitionMode = 1; % 1-SingleScan, 2-Accumulate, 3-Kinetics, 4-FastKinetics,

result = SetAcquisitionMode(acquisitionMode);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1, [bb 'String','Error during SetAcquisitionMode: ' result]);
    i4=0
end

readMode = 3; % 0-FVB, 1-MultiTrack, 2-RandomTrack, 3-SingleTrack, 4-Image
result = SetReadMode(readMode);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during SetReadMode: ' result]);
    i5=0;
end

shuttertype = 1;                        % 0-TTL low, 1-TTL high
shuttermode = 1;                        % 0-Auto, 1-Open, 2-Close
shutterclosingtime = 0;                 % time for shutter to close
shutteropeningtime = 0;                 % time for shutter to open
result = SetShutter(shuttertype, shuttermode, shutterclosingtime, shutteropeningtime);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during SetShutter: ' result]);
    i6=0;
end

% int type 0-speed of data readout, 1-speed of binning, 2-speed of dumping pixels; int index
result = SetHSSpeed(0, 0);  
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during setHSSpeed: ' result]);
    i8=0
end

result = SetHSSpeed(1, 0);  % Binning at 1 mu s
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during setHSSpeed: ' result]);
    i9=0
end

result = SetHSSpeed(2, 0);  % Dumping at 1 mu s
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during setHSSpeed: ' result]);
    i10=0
end

result = SetVSSpeed(0);  % vertical shift speed
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during setHSSpeed: ' result]);
    i11=0;
end

% Initializes the Single track position
centre = 79;%formerly 79  
height = 56;   
result = SetSingleTrack(centre, height);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during SetSingleTrack:' result]);
    i12=0;
end

exposureTime = 0.002;
result = SetExposureTime(exposureTime);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during SetExposureTime: ' result]);
    i13=0;
end

triggerMode = 0; % 0-Internal, 1-External, 6-ExternalStart
result = SetTriggerMode(triggerMode);
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during SetTriggermode: ' result]);
    i14=0;
end

exposure = 0;
accumulate = 0;
kinetic = 0;
[result, exposure, accumulate, kinetic] = GetAcquisitionTimings;
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,[bb 'String','Error during GetAcquisitionTimings: ' result]);
    i15=0;
end

n5=double(50);

set(handles.text2,'String',['Exposure Time(s): ' num2str(exposure)]);
set(handles.text3,'String',['Accumulate Time(s): ' num2str(accumulate)]);
set(handles.text4,'String',['Average Number: ' num2str(n5)]);

ii=i11+i12+i13+i14+i15+i4+i5+i6+i7+i8+i9+i10;

if ii==12
    Imea=1;
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb '--Measure is initialized ']);
else
    Imea=0;
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String','Measure is not successfully initialized ');
end

%======================================    
% to do: Initalize cooler & Temperature
tem1=0;
tem2=0;
tem0=0;

[result,tem1,tem2] = GetTemperatureRange;
[result, tem0] = GetTemperature;% get temperature

remind={'Current Temperature: ','Temperature Range:', 'Input Target Temperature: '};
line=1;
range=[num2str(tem1) '~' num2str(tem2)];
default={num2str(tem0),range,' '}
tem3dlg=inputdlg(remind,'Temperature Set',line,default); % key in the temperature
tem3=str2num(tem3dlg{3,1});
while ((tem3 < tem1) | (tem3 > tem2))
    uiwait(msgbox(['Desired temperature must be within range of ', num2str(tem1), ' and ', num2str(tem2), '. Please re-enter.'],'Oops!','modal'));
    tem3dlg=inputdlg(remind,'Temperature Set',line,default); % key in the temperature
    tem3=str2num(tem3dlg{3,1});
end

result = SetTemperature(tem3);
if (result ~= ANDOR. DRV_SUCCESS)
    bb= get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during set the temperature:' result]);
    i2=0
end

result=CoolerON;
if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during cooler on:' result]);
    i3=0
end

ii23=i2+i3;
if ii23==2
    while tem0>(tem3+1) %the extra +1 is to allow for a +/- 1 degree error in temp
        pause(10);
        [result, tem0] = GetTemperature;% get temperature
        bb = get(handles.edit1,'string');
        set(handles.edit1,'String',[bb '  current temperature: ' num2str(tem0)]);
    end
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb ' --The Cooler is set up --All initialization done' ]);
else
     bb = get(handles.edit1,'string');
     set(handles.edit1,'String',[bb '  The Cooler is not successfully set up ']);
       
end


handles.Isys=Isys;
handles.Imea=Imea;
    
guidata(hObject, handles); %updata handles


%=========================================================================
% --- Executes on button press in calibrateButton.
% --- CALIBRATE BUTTON -----
function calibrateButton_Callback(hObject, eventdata, handles)

global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));
load([progpath '/ref/calorgstd716983.mat']); %load the standard spectrum for calibration

Isys=handles.Isys;
Imea=handles.Imea;

switch and (Isys==1,Imea==1)
    case 0
        set(handles.edit1,'String','Pls. Initialize System Firstly ');
    case 1
        uiwait(msgbox('Get ready for first background measurement.','Measurement','modal'));
        dark1 = Measure;
        uiwait(msgbox('Get ready for second background measurement.','Measurement','modal'));
        dark2 = Measure;
        uiwait(msgbox('Get ready for calibration measurement (orange light).','Measurement','modal'));
        % Measurement is the title, modal is the formal of the graph
        
        orange = Measure;
        noisespc = dark1 - dark2;
        calibspc = orange - dark2;

        % evaluate noise level
        noise.mean = mean (noisespc);
        noise.std = std (noisespc);
        
        % display the noise level
        set (handles.edit1, 'string', [' Noise_Mean: ' num2str(noise.mean) ' Noise_std: ' num2str(noise.std)]);
        
        level = 5; % note the noise level
        calibspc(calibspc < (level * noise.std + noise.mean))=0;

        % display the signals and noise in spectrum
        axes (handles.axes2);
        set(handles.axes2,'visible', 'on');
        
        plot (calibspc, '-r', 'LineWidth',2);
        hold on 
        plot (calorgstd(:,4), '-b', 'LineWidth',2);
        hold off
        
        % choose points for calibration (15X2 points is suggested)
        bb = get(handles.edit1,'string');
        set (handles.edit1, 'String', [bb '  Choose Peaks for calibration' 'Red-measure' 'Blue-standard']);
        
        [xx,yy]=ginput;
        x1=xx(1:2:end); % the pixels in chip
        x2=int16(xx(2:2:end));
        y1=calorgstd(x2,3); % standard wavelength
        
        % find the wrong points and remove it
        n1=find(y1==0);
        x1(n1)=[];
        y1(n1)=[];
        
        if not(isempty(n1));
            bb = get(handles.edit1,'string');
            set (handles.edit1, 'String', [bb 'remove a wrong point']);
        end
        
        load([progpath '/ref/calfit.mat']);
        bb = get(handles.edit1,'string');
        set (handles.edit1, 'String', [bb '   before calfit=' num2str(calfit)]);

        calfit=polyfit(x1,y1,2);
        handles.calfit=calfit;

        
        bb = get(handles.edit1,'string');
        b0 = num2str(calfit);
        set (handles.edit1, 'String', [bb ' calibration is finished, coefficents: ' b0]);
        
        ans2=inputdlg('replace the existed calfit file? y/n');
        ans1=ans2{1,1};
        if ans1=='y'
            save ([progpath '/ref/calfit'],'calfit');
            bb = get(handles.edit1,'string');
            b0 = num2str(calfit);
            set (handles.edit1, 'String', [bb ' Spectrum Calibration is finished, coefficents: ' b0]);
        end
end

guidata(hObject, handles);

% --- Executes on button press in pauseButton.
% --- PAUSE BUTTON ---
function pauseButton_Callback(hObject, eventdata, handles)

uiwait(msgbox('Program is in Pause, Click OK to Continue','Pause','modal'));
return;


%=========================================================================
% --- Executes on button press in measureButton.
% --- MEASURE BUTTON ---
function measureButton_Callback(hObject, eventdata, handles)
Isys=handles.Isys;
Imea=handles.Imea;

% Get calfit coefficient
calfit=handles.calfit;
% Get the reference light
Ir0=handles.Ir0;

global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

switch and (Isys==1,Imea==1)
    case 0
        set(handles.edit1,'String','Pls. Initialize System Firstly ');
    case 1
        % load data
        % the measure range: 716.53-983.86nm; 
        
        % load absorptivity spectra (based on natural lograthim ); 
        load([progpath '/ref/puredataln.mat']);
        % puredataln:4X351, wavelength: 650-1000nm 
        % pure(1,:)-wavelength, pure(2,:)-pure water OD/cm/55.5M, 
        % pure(3,:)-Hb & pure(4,:)-HbO OD/mm/uM,
        % change the units to OD/cm/uM
        pure(2,:)=pure(2,:)*10^-6/55.5;
        pure(3,:)=pure(3,:)*10;
        pure(4,:)=pure(4,:)*10;
       
        % calculate the differential spectra
        diffpure(:,1)=(savgol(pure(2,:),29,4,2))'; % 2nd diff pure H2O
        diffpure(:,2)=(savgol(pure(3,:),29,4,2))'; % 2nd diff pure Hb
        diffpure(:,3)=(savgol(pure(4,:),29,4,2))'; % 2nd diff pure HbO
        
        % transfer the pixel number to wavelength, 1:1024 to wl0=715.4477-985.5299 nm
        wl0=polyval(calfit,1:1024) 
        wl6=find(((716:985)>=725)&((716:985)<=965));    
        wl7=725:965;
        % pure maters and signal range for water significant peak: 900-965nm/810-850nm
        wl1=find((wl7>=900)&(wl7<=965));
        wl2=find(((650:1000)>=900)&((650:1000)<=965));
        % wl1=find((wl7>=810)&(wl7<=850));
        % wl2=find(((650:1000)>=810)&((650:1000)<=850));
        % pure maters and signal range for Hb significant peak, 740-810 nm
        wl3=find((wl7>=740)&(wl7<=810)); % for signal
       
        wl4=find(((650:1000)>=740)&((650:1000)<=810)); % for pure maters
        
        % Measuring Background
        uiwait(msgbox('Background Measurement','Measuring Background','modal'));
        ans1='y';
        
        a1=['c', 'm', 'y', 'r', 'g', 'b','k'];
        noise3=zeros(100,1024);
        noise31=zeros(1,1024);

        while ans1=='y'
            uiwait(msgbox('Get Ready for Measuring Background','Measuring Background','modal'));
            % measure 100 times
            
            bb = get(handles.edit1,'string');
            set(handles.edit1,'String',[bb 'MeasuringTheBackground...    ']);
            
            for n2=1:100
                noise3(n2,:)=Measure;
            end
            
            noise31=sum(noise3,1)./100;
            a2=round((rand(1)+0.1)*7);

            axes (handles.axes2);
            set(handles.axes2,'visible', 'on');
            plot(wl0,noise31,a1(a2));
            hold on; 
            xlabel('Wavelength (nm)');
            ylabel('Count');
            title('Background Measurement');
            ylim([1030,1050]);
            % xlim([700,1000]);

            ans2=inputdlg('continue meausring Background? y/n');
            ans1=ans2{1,1};
        end

        % Input water concentration
        ans1=inputdlg('Water Concentration (rat:0.8;human:0.85)');
        c1=str2double(ans1{1,1})*55.5*10^6; %water--uM
       
        % Mearsuring Signal
        uiwait(msgbox('Signal Measurement','Measuring Signal','modal'));

        % input the average times
        ans2=inputdlg('input sample times for average');
        avg=str2double(ans2{1,1});
        
        set (handles.edit1, 'String', 'TestSampleRate-50DataPoints  ');
        
        % test the real sampling period
        tt=datestr(now);
        tt1=tt(12:20);
        signal3=zeros(avg,1024);
        AveSig1=zeros(50,1024);
        c5=zeros(50,1); 
        tt3=zeros(50,9);
        tt3=char(tt3);
        
        axes (handles.axes2);
        set(handles.axes2,'visible', 'on');
        hold off;
  
        for n1=1:50
            for n2=1:avg
                signal3(n2,:)=Measure;
            end
            
            AveSig1(n1,:)=((sum(signal3,1))./avg)-noise31(1,:);

            tt=datestr(now);
            tt3(n1,1:9)=tt(12:20);

            % calculation
            Im22=AveSig1(n1,:);

            % do ABS to remove negative data
            % Im22=abs(Im22);
            % nn2=find(Im22==0);
            % Im22(1,nn2)=0.0001;

            % Caculate the attenuation spectra, Ir0-RefLight.mat
            A1=log(Ir0)-log(Im22);
            % do filter for the Attenuation spectra
            fA1=spcfil(A1,1);

            % do decimate the fA1 to round wavelength points
            dfA0=interp1(wl0,fA1,716:985,'linear');

            % signal range to remove the filter window and delay,choose 725-970 nm
            dfA1=dfA0(wl6);
            % 2nd differential algorithm, specifying the number of points in filter 
            % (width), the order of the polynomial (order), and the derivative (deriv).

            diffSig1=(savgol(dfA1,21,4,2))'; % 2nd diff attenuation

            % now, directly use 'lsqnonneg'do the multi-regression
            % fitting in wl1: 900-965nm/810-870nm for the H2O pathlength
            [cw4(:,1), ew1(:,1), ew2(:,1)]=lsqnonneg(diffpure(wl2,1:2),diffSig1(wl1,1));

            pl(n1,1)=cw4(1,1)/c1; % pl: pathlength
            fitSig=diffpure(wl2,1:2)*cw4;

            % ew3(n4)=ew1(:,n4)^0.5/norm(diffSig1(wl1,n4));
            [c4(:,1), e1(:,1), e2(:,1)]=lsqnonneg(diffpure(wl4,:),diffSig1(wl3,1));

            c5(n1,1)=(c4(2,1)/pl(n1,1))*0.78-1.5; %c5(n1,:) concentration of Hb
            % diffA1(:,n4)=diffpure(wl4,:)*c4(:,n4);
            % e3(n4)=e1(:,n4)^0.5/norm(diffSig1(wl3,n4));


            plot(c5','*b');
            xlim([0,50]);
            ylim([0,50]);
            xlabel('data number');
            ylabel('[Hb] concentration (uM)');
            title('[Hb] Measurement');

        end
        
        tt=datestr(now);
        tt2=tt(12:20);
        bb = get(handles.edit1,'string');
        set (handles.edit1, 'String', [bb '50Data: ' tt1 '--' tt2 '  ']);

        
        % Signal Measurement
        ans2=inputdlg('input Sampling Time for One Point(s)','Measure Signal',1);       
        redt1=str2double(ans2{1,1});
     
        % input the record time
        ans2=inputdlg('input Total Record Time(min)','Measure Signal',1);
        redt0=str2double(ans2{1,1});

        redt2=redt0*60/redt1;
        redt3= 0:(redt0/(redt2-1)):redt0;
        
        bb = get(handles.edit1,'string');
        set (handles.edit1, 'String', [bb 'TotalRecordTime: ' num2str(redt0) 'min  ' 'OneDatumSamplingTime: ' num2str(redt1) 's  ']);
        
        signal3=zeros(avg,1024);
        Sig1=zeros(redt2,1024);
        c5=zeros(redt2,1);
        pl=zeros(redt2,1);
        t1=zeros(redt2,9);
        t1=char(t1);

        % input the head of file name
        ans2=inputdlg('input the file name','Measure Signal',1);
        namehead=ans2{1,1};

        uiwait(msgbox('Get ready for Measuring Signal','Measurement','modal'));
        
        axes (handles.axes2);
        set(handles.axes2,'visible', 'on');
        
        tt=datestr(now);
        bb = get(handles.edit1,'string');
        set (handles.edit1, 'String', [bb '[Hb]MeasureBegin- ' tt(12:20) ' ']);
        bb = get(handles.edit1,'string');
        
        for n1=1:redt2
            % measure 10 times to improve SNR
            for n2=1:avg
                signal3(n2,:)=Measure;
        
            end
    
        Sig1(n1,:)=(sum(signal3,1))./avg -noise31(1,:);

        %% record the sample time, 'char()' change ASCII to string
        tt=datestr(now);
        t1(n1,1:9)=tt(12:20);
        set (handles.edit1, 'String', [bb 'CurrentSamplingTime ' tt(12:20) ' ']);
        
        % calculation
        Im22=Sig1(n1,:);
       
       
    
        % do ABS to remove negative data
        %Im22=abs(Im22);
        %nn2=find(Im22==0);
        %Im22(1,nn2)=0.0001;
        % Caculate the attenuation spectra, Ir0 is from "load('RefLight.mat')"
        A1=log(Ir0)-log(Im22);   
        
        % do filter for the Attenuation spectra
        fA1=spcfil(A1,1);    
   
        % do decimate the fA1 to round wavelength points
        dfA0=interp1(wl0,fA1,716:985,'linear');
        % signal range to remove the filter window and delay,choose 725-970 nm
        dfA1=dfA0(wl6);   
        
        % 2nd differential algorithm, specifying the number of points in filter 
        % (width), the order of the polynomial (order), and the derivative (deriv).
        diffSig1=(savgol(dfA1,21,4,2))'; % 2nd diff attenuation

        % calculate the [Hb], [HbO]
        % e1: norm of the residual: norm((pure''*c-A'')^2). 
        % e2: the residual, A''-pure''*c.
        % now, directly use 'lsqnonneg'do the multi-regression
        % fitting in wl1: 900-965nm/810-870nm for the H2O pathlength
        [cw4(:,1), ew1(:,1), ew2(:,1)]=lsqnonneg(diffpure(wl2,1:2),diffSig1(wl1,1));
       
        pl(n1,1)=cw4(1,1)/c1; % pl: pathlength
        fitSig=diffpure(wl2,1:2)*cw4;
        % ew3(n4)=ew1(:,n4)^0.5/norm(diffSig1(wl1,n4));
        % disp(['pl= ' num2str(pl),'    fitting errors= ', num2str(ew3*100), '%']);
        % now, directly use 'lsqnonneg'do the multi-regression     
        [c4(:,1), e1(:,1), e2(:,1)]=lsqnonneg(diffpure(wl4,:),diffSig1(wl3,1));
       
        
        
        % add correction parameter (0.78,-1.5)
        c5(n1,1)=(c4(2,1)/pl(n1,1))*0.78-1.5; %c5(n1,:) concentration of Hb
        % diffA1(:,n4)=diffpure(wl4,:)*c4(:,n4);
        % e3(n4)=e1(:,n4)^0.5/norm(diffSig1(wl3,n4));
        
        handles.c5=c5;
        handles.t1=t1;
        guidata(hObject, handles);


        % plot will impact the execution speed, don't use the 'hold on'
        plot(redt3, c5(1:redt2,1)','-b', 'LineWidth',1);
        
        xlim([0,redt0+1]);
        ylim([0,100]);
        xlabel('Time (min)');
        ylabel('[Hb] concentration (uM)');
        title('[Hb] Measurement');
        
        end

        % save the data
        dataname=[progpath '/data/' namehead tt(1,13:14) tt(1,16:17) tt(1,19:20)];
        % save (dataname,'AveSig1','noise3','c5', 't1');
        save (dataname,'c5', 't1','redt3');
end

%========================================================================
% --- Executes on button press in exitButton.
% --- EXIT BUTTON ---
function exitButton_Callback(hObject, eventdata, handles)

load andor_errors.mat;
ANDOR = andor_errno; 

result = CoolerOFF;

if (result ~= ANDOR. DRV_SUCCESS)
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb 'Error during cooler off:' result]);
end

tem4=0;
[result,tem4] = GetTemperature;% get temperature

while tem4<5
    pause(5);
    [result, tem4] = GetTemperature;% get temperature
    bb = get(handles.edit1,'string');
    set(handles.edit1,'String',[bb '  current temperature: ' num2str(tem4)]);
end

AndorShutDown;
quit;


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in axes3_menu.
function axes3_menu_Callback(hObject, eventdata, handles)
% hObject    handle to axes3_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns axes3_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from axes3_menu


% --- Executes on selection change in axes2_menu.
function axes2_menu_Callback(hObject, eventdata, handles)
% hObject    handle to axes2_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

%Value:     1.0 - Path Length time-series; no slider
%           2.0 - dHb concentration time-series; no slider
%           3.0 - Attenuation spectra
%           4.0 - 2nd differential spectra
val = get(hObject, 'Value');
switch(val);
case 1
    axes (handles.axes2);
    HbLenStruct = load([progpath '/resi.mat'], 'HbLen');
    HbLen=HbLenStruct.HbLen;
    plot (1:length(HbLen(1,:)),HbLen(1,:), '-bs','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',3);
    xlim([0 length(HbLen(1,:))]);
    xlabel('Measurement Number');
    ylabel('Path Length (cm)');
    set(handles.axes2_slider,'Enable','off');
            
case 2
    axes (handles.axes2);
    HbLenStruct = load([progpath '/resi.mat'], 'HbLen');
    HbLen=HbLenStruct.HbLen;
    plot (1:length(HbLen(2,:)),HbLen(2,:), '-r','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',3);
    xlim([0 length(HbLen(2,:))]);
    xlabel('Measurement Number');
    ylabel('dHb Concentration (M)');
    set(handles.axes2_slider,'Enable','off');
    
case 3
    atteStruct = load([progpath '/resi.mat'], 'atte');
    atte = atteStruct.atte;
    calfit=handles.calfit;
    wavelength = polyval (calfit, 1:length(atte(1,:)));
    
    sliderStatus = get(handles.axes2_slider,'Enable')
    if (~strcmp(sliderStatus, 'on'));    
        set(handles.axes2_slider,'Enable','on');
        set(handles.axes2_slider,'Value',1);
        set(handles.axes2_slider,'Min',1);
        set(handles.axes2_slider,'Max',length(atte(:,1)));
    end
    
    axes (handles.axes2);
    index = int32(get(handles.axes2_slider,'Value'));
    plot (wavelength, atte(index,:), '-r', 'LineWidth',1);
    xlabel('Wavelength nm');
    ylabel('Attenuation (OD)');
    title('Attenuation Spectra');

case 4
    der2Struct = load([progpath '/resi.mat'], 'der2');
    der2 = der2Struct.der2;
    calfit=handles.calfit;
    wavelength = polyval (calfit, 1:length(der2(1,:)));
    
    sliderStatus = get(handles.axes2_slider,'Enable')
    if (~strcmp(sliderStatus, 'on'));    
        set(handles.axes2_slider,'Enable','on');
        set(handles.axes2_slider,'Value',1);
        set(handles.axes2_slider,'Min',1);
        set(handles.axes2_slider,'Max',length(der2(:,1)));
    end
    
    axes (handles.axes2);
    index = int32(get(handles.axes2_slider,'Value'));
    plot (wavelength, der2(index,:), '-r', 'LineWidth',1);
    xlabel('Wavelength nm');
    ylabel('2nd Derivative)');
    title('2nd Derivative Spectra');
    set(handles.axes2_slider,'Enable','on');
end

% Hints: contents = cellstr(get(hObject,'String')) returns axes2_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from axes2_menu

% --- Executes on slider movement.
function axes2_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes2_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

%Value:     1.0 - Path Length time-series; no slider
%           2.0 - dHb concentration time-series; no slider
%           3.0 - Attenuation spectra
%           4.0 - 2nd differential spectra
val = get(handles.axes2_menu, 'Value');
switch(val);
case 1
    axes (handles.axes2);
    HbLenStruct = load([progpath '/resi.mat'], 'HbLen');
    HbLen=HbLenStruct.HbLen;
    plot (1:length(HbLen(1,:)),HbLen(1,:), '-bs','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',3);
    xlim([0 length(HbLen(1,:))]);
    xlabel('Measurement Number');
    ylabel('Path Length (cm)');
    set(handles.axes2_slider,'Enable','off');
            
case 2
    axes (handles.axes2);
    HbLenStruct = load([progpath '/resi.mat'], 'HbLen');
    HbLen=HbLenStruct.HbLen;
    plot (1:length(HbLen(2,:)),HbLen(2,:), '-r','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',3);
    xlim([0 length(HbLen(2,:))]);
    xlabel('Measurement Number');
    ylabel('dHb Concentration (M)');
    set(handles.axes2_slider,'Enable','off');
    
case 3
    atteStruct = load([progpath '/resi.mat'], 'atte');
    atte = atteStruct.atte;
    calfit=handles.calfit;
    wavelength = polyval (calfit, 1:length(atte(1,:)));
    
    axes (handles.axes2);
    index=int32(get(handles.axes2_slider,'Value'));
    plot (wavelength, atte(index,:), '-r', 'LineWidth',1);
    xlabel('Wavelength nm');
    ylabel('Attenuation (OD)');
    title('Attenuation Spectra');

case 4
    der2Struct = load([progpath '/resi.mat'], 'der2');
    der2 = der2Struct.der2;
    calfit=handles.calfit;
    wavelength = polyval (calfit, 1:length(der2(1,:)));
    
    axes (handles.axes2);
    index=int32(get(handles.axes2_slider,'Value'));
    plot (wavelength, der2(index,:), '-r', 'LineWidth',1);
    xlabel('Wavelength nm');
    ylabel('2nd Derivative)');
    title('2nd Derivative Spectra');
    set(handles.axes2_slider,'Enable','on');
end

% --- Executes on slider movement.
function axes3_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes3_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on button press in calibrateIntensityButton.
function calibrateIntensityButton_Callback(hObject, eventdata, handles)
% hObject    handle to calibrateIntensityButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calibrateIntensityButton
Isys=handles.Isys;
Imea=handles.Imea;
calfit=handles.calfit;

global progpath;
[progpath, ~, ~] = fileparts(mfilename('fullpath'));

switch and (Isys==1,Imea==1)
    case 0
        set(handles.edit1,'String','Pls. Initialize System Firstly ');
    case 1
        % Intensity Reference measurement
        % Measure Noise
        dark3=zeros(100,1024);
        uiwait(msgbox('Get ready for background measurement.','Intensity Reference','modal'));
        
        ans2=inputdlg('Are any neutral density attenuation filters being used? y/n?');
        ans1=ans2{1,1};
        if ans1=='y'
            ans2=inputdlg('What is the OD rating of the filter?');
            filterOD=str2double(ans2{1,1});
        else
            filterOD=0;
        end
        
        bb = get(handles.edit1,'string');
        set(handles.edit1,'String',[bb 'System is measuring the background... ']);
        
        for n1=1:100
            dark3(n1,:) = Measure;
            pause(1);
        end

        % display the noise spectrum
        axes (handles.axes2);
        set(handles.axes2,'visible', 'on');
        wl0=polyval(calfit,1:1024);
        plot(wl0, dark3(1:20:end,1:1024));

        load([progpath '/ref/RefLight.mat']);
        Ir1=Ir0;
        
        % Measure Intensity Reference
        I1=zeros(100,1024);
        uiwait(msgbox('Get ready for Reference Measurement-Empty & High Intensity','Intensity Reference','modal'));

        bb = get(handles.edit1,'string');
        set(handles.edit1,'String',[bb 'System is measuring the Reference Light... ']);
        
        for n1=1:100
            I1(n1,:) = Measure;    
            pause(1);
        end
        plot(wl0, I1(1:5:end,1:1024));
        xlabel('Wavelength (nm)');
        ylabel('Count');
        title('Reference Measurement');
        hold on;

        I1=I1-dark3;
        Ir0=sum(I1,1)./100; %for averaging
        Ir0=Ir0.*(10^filterOD); 

        plot(wl0, Ir0,'--r','LineWidth',2);
        plot(wl0, Ir1,'--b','LineWidth',1);
        hold off;
        
        ans2=inputdlg('red-new RefLight; blue-existed RefLight; replace the existed RefLight file? y/n');
        ans1=ans2{1,1};
        if ans1=='y'
            save ([progpath '/ref/RefLight'],'Ir0');
            handles.Ir0=Ir0;
            bb = get(handles.edit1,'string');
            set (handles.edit1, 'String', [bb ' Reference Measurement is finished ']);
        end
        
end

guidata(hObject, handles);


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveButton

c5=handles.c5;
t1=handles.t1;
tt=datestr(now);
dataname=['SN' tt(1,13:14) tt(1,16:17) tt(1,19:20)];
save (dataname,'c5', 't1');
