function p = plotInit
%% Line styles specifiers
p.line   = {'--' '-' '-.'}; % ':'

%% Marker specifiers
p.marker = {
    'n' % None
    's' % Square     
    'o' % Circle             
    '^' % Upward  -pointing triangle
    'v'	% Downward-pointing triangle    
    '>'	% Right   -pointing triangle
    '<'	% Left    -pointing triangle       
    'd' % Diamond
    'p' % Five-pointed star (pentagram)        
    }; 
%   '+' % Plus sign
%   '.' % Point
%   '*' % Asterisk
%   'x' % Cross
%   'h' % Six-pointed star (hexagram)
%% Color specifiers
p.color = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.25  0.25  0.25
    0.95  0.95  0.00
    0.75  0.75  0.75    
    0.00  1.00  0.00
    0.54  0.63  0.22
    0.34  0.57  0.92
    0.75  0.75  0.00
    0.88  0.75  0.73
    0.10  0.49  0.47
    0.66  0.34  0.65
    0.99  0.41  0.23
    1.00  0.10  0.60
    0.25  0.25  0.75    
    0.76  0.57  0.17    
    0.75  0.25  0.25    
    ];