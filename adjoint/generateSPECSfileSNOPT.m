function snopt = generateSPECSfileSNOPT( oc, pulse, snopt )

% Determine if there is a savedir already
% timestr = strcat('__', char(datetime('now', 'format','yyMMddHHmmss')) );
if isfield( oc, 'saveDir' )
    saveDirTemp = oc.saveDir;
else
    currDir = fileparts( mfilename( "fullpath" ) );
    tld = strsplit( currDir, 'multiphoton' );
    tld = tld{ 1 };
    mpd = fullfile( tld, 'multiphoton' );
    saveDirTemp = fullfile( mpd, "data", "opt",...
        strcat( string( datetime('now', 'format','yyyy-MM-dd') ), "-SNOPT-opt-file" ) );
end
SPECSfile = fullfile( saveDirTemp, strcat( "SPECS", ".spc" ) );
if oc.saveResult
    PRINTfile = fullfile( saveDirTemp, strcat( "PRINT", ".out" ) );
    % SUMMARYfile = fullfile( saveDirTemp, strcat( "SUM", ".out" ) );
end

% Determine opt name
if isfield( snopt, 'name' )
    pulseName = snopt.name;
else
    pulseName = strcat( pulse.name, '-', pulse.type );
end

% Iterate over possible SPECS to include:
specStrArray = sprintf( "Begin\t%s", pulseName );

% Major Print Level
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "majorprintlevel", "Major Print level",...
    1, snopt );

% Minor Print Level
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "minorprintlevel", "Minor Print level",...
    0, snopt );

% Summary Frequency
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "summaryfrequency", "Summary frequency",...
    500, snopt );

% Print Frequency
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "printfrequency", "Print frequency",...
    500, snopt );

% Summary file
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "summaryfile", "Summary file",...
    1, snopt );

% Solution on print file
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "solution", "Solution",...
    'yes', snopt );

% Major feasibility tolerance
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "majorfeasibilitytolerance", "Major feasibility tolerance",...
    1e-9, snopt );

% Minor feasibility tolerance
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "minorfeasibilitytolerance", "Minor feasibility tolerance",...
    1e-6, snopt );

% Major optimality tolerance
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "majoroptimalitytolerance", "Major optimality tolerance",...
    1e-6, snopt );

% Minor optimality tolerance
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "minoroptimalitytolerance", "Minor optimality tolerance",...
    1e-6, snopt );

% Time limit
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "timelimit", "Time limit",...
    3600, snopt );

% Verify gradients
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "verifylevel", "Verify level",...
    -1, snopt );

% Scale option (scaling parameters)
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "scaleoption", "Scale option",...
    1, snopt );

% Major iterations limit
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "majoriterationslimit", "Major iterations limit",...
    500, snopt );

% Minor iterations limit
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "minoriterationslimit", "Minor iterations limit",...
    500, snopt );

% Major step limit
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "majorsteplimit", "Major step limit",...
    2.0, snopt );

% Derivative level: 3 means we have calculated gradients for all
% derivatives
[ specStrArray, snopt ] = addSpecsFile( specStrArray, "derivativelevel", "Derivative level",...
    3, snopt );

specStrArray = [ specStrArray; sprintf( "End\t%s", pulseName ) ];

if oc.saveResult
    snopt.printfile = char( PRINTfile );
    % snopt.summaryfile = char( SUMMARYfile );
end
snopt.specsfile = char( SPECSfile );

if ~isfolder( saveDirTemp )
    mkdir( saveDirTemp )
end
saveTxtFromStrArray( SPECSfile, specStrArray );

end

%% Helper Function
% ----------------------------------------------------------------------- %
function [ specStrArray, snopt ] = addSpecsFile( specStrArray, snoptStattribute, printfiletext, defaultval, snopt )
if ~isfield( snopt, snoptStattribute )
    snopt.(snoptStattribute) = defaultval;
end
specStrArray = [ specStrArray; sprintf( strcat("\t", printfiletext,"\t\t\t%s"), num2str( snopt.(snoptStattribute) ) ) ];

snopt = rmfield( snopt, snoptStattribute );
end
% ----------------------------------------------------------------------- %