function spectrum = Measure ()
global ANDOR

% execute the measurement
result = StartAcquisition;
if (result ~= ANDOR. DRV_SUCCESS)
    if (result == ANDOR. DRV_ACQUIRING)
        fprintf ('Still acquiring\n');
    else
        error (['Error during StartAcquisition: ' result]);
    end
end

status = 0;
[result, status] = AndorGetStatus;
fprintf ('Acquiring data ');
while (status  == ANDOR. DRV_ACQUIRING)
    [result, status] = AndorGetStatus;
    pause (0.01);        % change to appropriate time
    fprintf ('.');
end;
fprintf (' done.\n');

if (status ~= ANDOR. DRV_IDLE)
    error (sprintf ('Error during acquisition: %s (%i)', translate_errno (status), status));
end

size = 1024;   % change according to measurement parameters          
ResultArray = zeros (1, size);

[result, ResultArray]= GetAcquiredData( size);

if (result ~= ANDOR. DRV_SUCCESS)
    error (sprintf ('Error fetching data: %i', result));
end

ResultArrayRow = reshape(ResultArray, [1, size]);

spectrum = double(ResultArrayRow);