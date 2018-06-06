function str = getOrder(i)
% Function: Return a string to represent the order of the current number
%   examples: 
%       1 -> 1st, 2 -> 2nd, 3 -> 3rd
%       4 -> 4th, 11 -> 11st ...        
% 
% Usage:
% 
%       str = getOrder(i)
%   where:
%       str - the order as a string
%       i - the number
% 
% Author: Zhen Zhang
% Last modified: 6 Jun. 2018

% Get the lowest digit
ldigit = mod(i, 10);

switch ldigit
    case 1
        str = [num2str(i) 'st'];
    case 2
        str = [num2str(i) 'nd'];
    case 3
        str = [num2str(i) 'rd'];
    otherwise
        str = [num2str(i) 'th'];
end

end