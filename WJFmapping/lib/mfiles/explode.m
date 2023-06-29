% exploded = explode(separator, stringToExplode, itemWrapper)
%
% Import of a very cool php function.  Also optionally allows an itemWrapper,
% which is useful for removing quotes around the cellArray items.  Currently,
% itemWrapper can only be one character long.
%
% Easiest to use examples:
%   explode(',', {'one,two,three') produces {'one', 'two', 'three'}
%   explode(';', '"one";"two";"three"', '"') produces {'one', 'two', 'three'}
%
% RETURNS: a cell array that is the PHP-style explosion of the string
%
% This function is part of mrvista: http://white.stanford.edu/software/
% written by B. Dougherty, 2001
%  (c) Board of Trustees, Leland Stanford Junior University


function exploded = explode(separator, stringToExplode, itemWrapper)

if (~exist('separator','var') | ~exist('stringToExplode','var') ...
        | isempty(stringToExplode) | ...
        (exist('itemWrapper','var') & length(itemWrapper)>1))
    help mfilename;
    return;
end
if (~exist('itemWrapper','var'))
    itemWrapper = '';
end

sep = findstr(separator, stringToExplode);
if (isempty(sep))
    exploded{1} = stringToExplode;
else
    exploded{1} = stringToExplode(1:sep(1)-1);
    for i = [1:length(sep)-1]
        exploded{i+1} = stringToExplode(sep(i)+1:sep(i+1)-1);
    end
    exploded{length(exploded)+1} = stringToExplode(sep(end)+1:end);
end

if (~isempty(itemWrapper))
    for i = [1:length(exploded)]
        if (exploded{i}(1) == itemWrapper)
            exploded{i} = exploded{i}(2:end);
        end
        if (exploded{i}(end) == itemWrapper)
            exploded{i} = exploded{i}(1:end-1);
        end
    end
end