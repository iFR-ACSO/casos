function [xs,ys,info] = sedumi2sdpam(As,bs,cs,K,opts)
% function [x,y,info] = sedumi2sdpam(A,b,c,K,opts);
%
% DESCRIPTION
%   This function solves a cone optimization specified in Sedumi format
%   using the SDPAM solver.  The function converts the inputs from
%   Sedumi format to SDPAM format, solves the problem with SDPAM, and
%   converts the SDPAM outputs back into Sedumi output format.
%
% INPUTS
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See
%       Sedumi documentation for more details.  SOCP constraints
%       are not handled by SDPAM.
%   opts: Options structure in SDPAM format.
%
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: SDPAM solver info.  See SDPAM documentation.
%
% SYNTAX
%   [x,y,info] = sedumi2sdpam(A,b,c,K,opts);
%
% EXAMPLE

% 3/6/08   PJS Initial Coding
% 11/18/10 PJS Update to 7.3.1 (Directly call executable, readspdaresult)

% Initalize fields of K
if ~isfield(K,'f')
    K.f = [];
end
if length(K.s)==1 && K.s==0
    K.s = [];
end

% Default print option
if ~isfield(opts,'print')
    opts.print = 'no';
end

%XXXX If opts.lambdaStar is too small then sdpam might only quit after
% one iteration. Need to increase lambdaStar if this happens.

% SDPAM does not handle free variables in the primal form. Split free
% variables into positive/negative parts.
[A2,b2,c2,K2] = splitfreevars(As,bs,cs,K);

if exist('mexSedumiWrap','file')
    % Call SDPAM Sedumi Wrapper if compiled mex file exists in path
    pars = [];
    [xs,ys,info]=sedumiwrap(A2,b2,c2,K2,pars,opts);
elseif ispc
    % Path to SDPA Executable
    SDPAFullPath = 'C:\Users\pseiler\Desktop\MatlabLocal\SDPSolvers\sdpam\sdpa.7.3.1.src\bin\';
    
    % Create temporary file names for SPDA input/output/parameters
    datfile = 'sedumi2sdpamtempfile.dat-s';
    if isfield(opts,'resultFile') && ~isempty(opts.resultFile)
        resultfile = opts.resultFile;
    else
        resultfile = 'sedumi2sdpamtempfile.result';
    end
    
    % Two ways to calling SDPA:
    % A) Call SDPA with its full path name from the current Matlab working
    %    directory. In this case, SDPA expects a param.sdpa to exist in
    %    the current working directory (it doesn't look in the folder
    %    where the SDPA exectuable exists). The .result file will be
    %    generated in the current working directory.
    % B) Change into the directory where the SDPA executable exists.
    %    Create the param.sdpa file in this folder, run SDPA (generating
    %    the result file), and then change back into the working directory.
    % The second option seems preferable since all input/output/param
    % files will be created in one location.  The first option might
    % result in these files cluttering user folders, e.g. the files will
    % not get deleted if the function exits early due to an error.
    p = pwd;
    cd(SDPAFullPath);
    
    try
        % Save Sedumi format data to a temporary SDPA sparse data file
        SedumiToSDPA(datfile,A2,b2,c2,K2)
        
        % Create options file
        LOCALwriteparam(SDPAFullPath,opts);
        
        % Solve the optimization with SDPA
        if strcmp(opts.print,'display')
            [status,result]=dos([SDPAFullPath 'sdpa ' datfile ' ' ...
                resultfile],'-echo');
        else
            [status,result]=dos([SDPAFullPath 'sdpa ' datfile ' ' ...
                resultfile]);            
        end
        
        % Read results from SDPA output file
        [xs,ys,info] = LOCALreadsdparesult(resultfile);
        
        % Delete created files
        ds1=dos(['del ' datfile]);
        ds2=dos(['del ' resultfile]);
    catch
        % Clean up and rethrow error
        ds1=dos(['del ' datfile]);
        ds2=dos(['del ' resultfile]);
        cd(p);
        rethrow(lasterr);
    end
    
    % Change back to user's working directory
    cd(p);
else
    error('Cannot find SDPA Mex Executable');
end

% Recombine positive/negative parts of the free variables
if ~isempty(K.f)
    xs = [xs(K.f+1 : 2*K.f) - xs(1:K.f); xs(2*K.f+1:end)];
end



%----------------------------------------------------------------------
% Local function to write spda.param options file
function LOCALwriteparam(SDPAFullPath,opts)

% Error checking on opts and fill in default values
% Call to maxNumCompThreads gives a warning in param file. 
warning off;  
defaultopts = param;
opts = param(opts);
warning on;

% Text for sdpa.param file
paramtext = {...
    '100',' unsigned int maxIteration;','maxIteration';...
    '1.0E-7',' double 0.0 < epsilonStar;','epsilonStar';...
    '1.0E2',' double 0.0 < lambdaStar;','lambdaStar';...
    '2.0',' double 1.0 < omegaStar;','omegaStar';...
    '-1.0E5',' double lowerBound;','lowerBound';...
    '1.0E5',' double upperBound;','upperBound';...
    '0.1',' double 0.0 <= betaStar <  1.0;','betaStar';...
    '0.2',' double 0.0 <= betaBar  <  1.0, betaStar <= betaBar;','betaBar';...
    '0.9',' double 0.0 < gammaStar  <  1.0;','gammaStar';...
    '1.0E-7',' double 0.0 < epsilonDash;','epsilonDash';...
    '%+8.3e',' char*  xPrint   (default %+8.3e,   NOPRINT skips printout)','xPrint';...
    '%+8.3e',' char*  XPrint   (default %+8.3e,   NOPRINT skips printout)','XPrint';...
    '%+8.3e',' char*  YPrint   (default %+8.3e,   NOPRINT skips printout)','YPrint';...
    '%+10.16e',' char*  infPrint (default %+10.16e, NOPRINT skips printout)','infPrint';...
};

% Some options are not available for direct exectuable call
if isfield(opts,'isSymmetric') && opts.isSymmetric ~= 0
    warning(['isSymmetric option only available with SDPA mex interface. ' ...
        'Setting isSymmetric to default value.']);
    opts.isSymmetric = defaultopts.isSymmetric;
end
if isfield(opts,'isDimacs') && opts.isDimacs ~= 0
    warning(['isDimacs option only available with SDPA mex interface. ' ...
        'Setting isDimacs to default value.']);
    opts.isDimacs = defaultopts.isDimacs;
end
if isfield(opts,'NumThreads') && opts.NumThreads ~= defaultopts.NumThreads
    warning(['NumThreads option only available with SDPA mex interface. ' ...
        'Setting NumThreads to default value.']);
    opts.NumThreads = defaultopts.NumThreads;
end

% Update text for sdpa.param file based on user options
Nopts = size(paramtext,1);
for i1=1:Nopts   
    if isfield(opts,paramtext{i1,3})
        val = opts.(paramtext{i1,3}); 
        if ischar(val)
            paramtext{i1,1} = val;
        else
            paramtext{i1,1} = num2str(val,'%0.16g');
        end
    end
end

% Write options text to sdpa.param file
valtxt = char(paramtext(:,1)); 
opttxt = char(paramtext(:,2)); 
fid = fopen([SDPAFullPath 'param.sdpa'],'wt');
%fid = fopen('junk.param.sdpa','wt'); 
for i1=1:Nopts, 
    cnt = fwrite(fid,valtxt(i1,:),'char'); 
    cnt = fwrite(fid,opttxt(i1,:),'char');
    cnt = fwrite(fid,char(10),'char'); 
end; 
fclose(fid);




%----------------------------------------------------------------------
% Local function to read in results from the SDPA output file
% and convert to Sedumi output format.
function [x,y,info] = LOCALreadsdparesult(filename)

% Open SDPA output file for reading
fid = fopen(filename,'r');

% Find start of optimization summary information
tline = LOCALstrsearch(fid,'phase');

% Create structure with basic summary information
info = [];
while 1
    idx = strfind(tline,'=');
    field = tline(1:idx-1);
    field( field=='.' ) = [];
    field( field==' ' ) = [];
    
    value = tline(idx+1:end);
    value( value==' ' ) = [];
    
    if strcmp( field, 'phasevalue')
        info.(field) = value;
    else
        info.(field) = eval(value);
    end
    
    tline = fgetl(fid);
    idx = findstr('** Para',tline);
    if ~isempty(idx)
        break
    end
end

% Find start of approximate optimal decision variable xVec
tline = LOCALstrsearch(fid,'xVec');

% Create Sedumi dual var y ( = SDPA primal var xVec)
tline = fgetl(fid);
tline(1) = '[';
tline(end) = ']';
y = eval(tline);
y = y(:);

% Find start of approximate optimal decision variable yMat
tline = LOCALstrsearch(fid,'yMat');

% Create Sedumi primal var x ( = SDPA dual var yMat)
x = [];
while 1
    tline = fgetl(fid);
    idx = strfind(tline,'loop time');
    if ~isempty(idx)
        break
    end
    
    sidx = strfind(tline,'{');
    eidx = strfind(tline,'}');
    if ~isempty(sidx) && ~isempty(eidx)
        tmp = tline(sidx(end):eidx(1));
        tmp(1) = '[';
        tmp(end) = ']';
        x = [x eval(tmp)];
    end
end
x = x(:);

% Close SDPA output file
fclose(fid);


%----------------------------------------------------------------------
% Local function to sequentially search lines of SDPA output file
% for a specified str.  The entire line containing str is returned.
function tline = LOCALstrsearch(fid,str)

while 1
    tline = fgetl(fid);
    idx = strfind(tline,str);
    if ~isempty(idx)
        break
    end
end
