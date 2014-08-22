
if ~exist ('path', 'var')  %check the variable 'path', if not exist, then
    path = 'andor/';
end

configfilepath = path;

if ispc()
    librarypath = [path 'win32/'];
elseif isunix()
    librarypath = [path 'unix/'];
end

% Loads the functions defined in header file 'atmcd32d.h' and found in library 'atmcd32d.dll' into MATLAB.
if ispc()
    loadlibrary ([librarypath 'atmcd32d.dll'], [librarypath 'atmcd32d.h']);
elseif isunix()
    %insert unix library loading code here
end


global ANDOR;% make the variable ANDOR can be used by multi program
load andor_errors.mat;
ANDOR = andor_errno; 

result=calllib('atmcd32d','CoolerOFF');

if (result ~= ANDOR. DRV_SUCCESS)
    errordlg(['Error during cooler off:' result]);
end

tem4=libpointer('int32Ptr',int32(1));
result=calllib('atmcd32d', 'GetTemperature', tem4);% get temperature

disp('Waiting for the Temperature arrive 5' )

while tem4.value<5
    pause(5);
    result=calllib('atmcd32d', 'GetTemperature', tem4);% get temperature
    disp(['Current Temperature: ' num2str(tem4.value)]);
end
