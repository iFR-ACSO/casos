function disp_matrix(obj,delim,out)
% Print matrix to command line (display).

% Mostly taken from multipoly, version history:
% 6/7/2002: PJS  Initial Coding
% 6/9/2002: PJS  Use char conversion and display matrices

if nargin < 3
    % get String Representation of object
    out = str(obj);
end
if nargin < 2
    % use default delimiter
    delim = '[]';
end
sza = size(obj);

% Compute sizes of character arrays and max size
nchar = cellfun('size',out,2);
maxr = max(nchar,[],1);
numc = sum(maxr);
rootprops = fieldnames(get(0));
if any(strcmp(rootprops,'CommandWindowSize'))
    ws = get(0,'CommandWindowSize');
else
    % Older Versions of matlab don't allow access to the
    % window size.  Modify window size (ws) if matrix
    % line breaks are not in the right spot.
    ws = 80;
end
maxchar=ws(1)-6-3*(sza(2)-1);

if isempty(obj)
    fprintf('Empty %s: %dx%d\n',class(obj),sza(1),sza(2));
    
elseif max(numc) < maxchar
    % Display as matrix if all rows can fit on screen
    for i1 = 1:sza(1)
        if all(sza==[1 1])
            d = '';
        else
            d = delim(1);
        end
        for i2 = 1:sza(2)
            d = [d blanks(maxr(i2)-nchar(i1,i2)) out{i1,i2}];
            if i2~=sza(2)
                d = [d ', '];
            elseif ~all(sza==[1 1])
                d = [d delim(2)];
            end
        end
        disp(d);
    end
    
else
    % Else display entry-by-entry going down columns
    for i2 = 1:sza(2)
        if min(sza) > 1
            disp(['----  Column ' int2str(i2) ' ----------'])
            disp(' ')
        end
        for i1 = 1:sza(1)
            % Display the (i1,i2) entry
%             if ~all(sza==[1 1])
%                 disp([inputname(1) '(' int2str(i1) ',' int2str(i2) ')  = ']);
%             else
%                 disp([inputname(1) ' = ']);
%             end
            
            % Break lines at +,-, or *
            sr = out{i1,i2};
            while ~isempty(sr) %length(sr)>0
                if length(sr) < (ws(1)-6)
                    disp(['  ' sr]);
                    sr = [];
                else
                    idx1 = sort([strfind(sr,'-') strfind(sr,'+') strfind(sr,'*')]);
                    idx1 = [setdiff(idx1,1) length(sr)];
                    %idx2 = max(find(idx1< (ws(1)-6) ));
                    idx2 = find(idx1< (ws(1)-6), 1, 'last' );
                    if isempty(idx2)
                        disp(['  ' sr(1:idx1(1)-1)]);
                        sr = sr(idx1(1):end);
                    else
                        disp(['  ' sr(1:idx1(idx2)-1)]);
                        sr = sr(idx1(idx2):end);
                    end
                end
            end
            if ~all([i1 i2]==sza)
                disp(' ');
            end
        end
    end
    
end

% if ~compact
%     disp(' ');
% end
