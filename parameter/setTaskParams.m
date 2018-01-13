function [ pTask, pI ] = setTaskParams( sTaskName, sDimension, sExampleName, sPath, sBaseRunID )
%% Task parameters
pTask.name         = sTaskName;
pTask.sDimension   = sDimension;
pTask.sExampleName = sExampleName;
pTask.sPath        = [ sPath '\' sTaskName '\' sDimension '\' sExampleName ]; % [ sPath '\' 'local' '\' sDimension '\' sExampleName ]; % 
pTask.sBaseRunID   = sBaseRunID;
% addpath( [ cd '/parameter/function/' sDimension ] );
switch pTask.name
    case 'Schrodinger'
        pTask.hbar    = 1;
        pTask.rho_inf = 1;
        switch sDimension % set values for dimentional parameters
            case '1D'
                pTask.B_inf = 2;
            case '2D'
                pTask.B1_inf = 2; pTask.B2_inf = pTask.B1_inf;
        end
        switch pTask.sExampleName
            case 'EX01r'
                pTask.X      = 1.6;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 0.8;       % x0
                pI.exAlpha   = 1/120;     % alpha
                pI.exKoef    = 100;       % k
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        pI.exY0      = pI.exX0;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX01'
                pTask.X      = 1.5; %TODO: currently unused! Must to be checked: pTask.X <= pCalc.X_0
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 0.8;       % x0
                pI.exAlpha   = 1 / 120;   % alpha
                pI.exKoef    = 100;       % k
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        pI.exY0      = pI.exX0;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX01a'
                pTask.X      = 1.5;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 0.8;       % x0
                pI.exAlpha   = 1/60;      % alpha
                pI.exKoef    = 1;         % k
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        pI.exY0      = pI.exX0;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case {'EX02','EX02t'}
                pTask.X      = 1.75;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        % initial wave parameters
                        pI.exY0      = pI.exX0;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX02m'
                pTask.X      = 1.75;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = 2.8; % 2011-CMP: 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX02n' % Splitting-in-potential Ducomet, Zlotnik, Zlotnik - 2012 - 2D
                pTask.X      = 1.75;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1/120;
                pI.exKoef    = 15;
                switch sDimension
                    case '2D'
                        pTask.Y      = 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX02k'
                pTask.X      = 4;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0    = 0.7;
                pI.exAlpha = 1/120;
                pI.exKoef  = 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX02_VATRO'
                pTask.X      = 10;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 5;
                pI.exAlpha   = 1/4;
                pI.exKoef    = 0;
                switch sDimension
                    case '2D'
                        pTask.Y      = 20;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX03'
                pTask.X      = 1.75;
                pTask.V_Linf = 1000; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                pI.exLambda  = - pTask.V_Rinf / ( pTask.hbar * pTask.rho_inf );
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        % initial wave parameters
                        pI.exY0      = pI.exX0;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX04' %% Arnold,Ehrhardt,Schulte - Numerical Simulation of Quantum Waveguides (4.1)
                pTask.X      = 1;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 3/4;
                pI.exAlpha   = 1 / 480;
                pI.exKoef    = 120;
                switch sDimension
                    case '2D'
                        pTask.Y      = 1;
                        % initial wave parameters
                        pI.exY0      = 1/4;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX05'
                pTask.X      = 1.25;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 0.625;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = 3;
                        % initial wave parameters
                        pI.exY0      = 1;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = 10;
                end
            case 'EX10'
                pTask.X      = 1.75;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                % initial wave parameters
                pI.exX0      = 0.75;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = - 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case {'EX20','EX20r'} %% Review - 2008
                pTask.X      = 18; % 24; % for the reference solution
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x    = [ 15 15.5 16 16.5 17 ];
                pTask.V_Q    = [ 25/2 5/2 0 25/2 ];
                % initial wave parameters
                pI.exX0      = 9;
                pI.exAlpha   = 1;
                pI.exKoef    = sqrt(7);
                switch sDimension
                    case '1D'
                        pTask.B_inf  = 1;
                    case '2D'
                        pTask.B1_inf = 1; pTask.B2_inf = pTask.B1_inf;
                        pTask.Y      = 2.8;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX21' % delta-function
                pTask.X      = 1.75;
                % potential parameters
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_j0   = 256;
                pTask.V_h    = 1 / 160;
                pTask.V_x    = [ (pTask.V_j0 - 1) * pTask.V_h pTask.V_j0 * pTask.V_h ];
                pTask.V_Q    = - 100 * 0.1;
                pTask.V_Q    = [ pTask.V_Linf pTask.V_Q/pTask.V_h pTask.V_Rinf ];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '2D'
                        pTask.Y      = 2.8;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
            case 'EX22' % potential barier / well % Zlotnik - Vestnik MPEI - 2010
                pTask.X      = 2.8;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x    = [1.6 1.7];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.V_Q    = 1000 * 1.5; % 2D: 1, 1.5, 4
                        pTask.Y      = 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22a'
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x    = [1.6 1.7];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension % set values for dimentional parameters
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.Y      = 1*2.8;
                        pTask.V_Q    = 1000 * 1.5; % 2D: 1, 1.5, 4
                        % D            = pTask.Y / 2^5 / 2; % pTask.Y / 2^2; % 
                        pTask.V_y    = [ pTask.Y/4 3*pTask.Y/4 ]; % pTask.Y / 2 + [ -D D ]; % [ 0+D pTask.Y-D ]; % 
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22m' % potential barier / well % Zlotnik, Zlotnik - KRM - 2012
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x    = [1.6 1.7];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.Y      = 2.8; % 2.75;
                        pTask.V_Q    = 1000 * 1.5; % 2D: 1, 1.5, 4
                        pTask.V_y    = [0 pTask.Y];
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22r' % potential barier / well % Zlotnik, Zlotnik - KRM - 2012
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x    = 0.4+[1.6 1.7];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.Y      = 2.8; % 2.75;
                        pTask.V_Q    = 1000 * 1.5; % 2D: 1, 1.5, 4
                        pTask.V_y    = [0 pTask.Y];
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22b' % potential-well -Q\chi_{(a,b)}
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_Q    = - 5 * 1000;
                pTask.V_x    = [1.8 2.0];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX32' % potential-well -Q(\tanh(a*(x-xL))+1)(\tanh(a*(xR-x))+1)
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_Q    = 4000; % 1000; % 
                pTask.V_alpha= 50;
                pTask.V_x    = [1.9 2.0];
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
            case 'EX22n' % potential barier / well % Zlotnik, Zlotnik - KRM - 2012
                pTask.X      = 3;
                pTask.V_Linf = 0;
                pTask.V_x    = [2.3 5.3]; % [1.6 5.6];
                % initial wave parameters
                pI.exX0      = 0.7; % 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.Y      = 2.75;
                        pTask.V_Q    = 1000 * 1; % 2D: 1, 1.5, 4
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Rinf = pTask.V_Q; pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22k'
                pTask.X      = 4;
                pTask.V_Linf = 0;
                pTask.V_x    = [2.3 5.3]; % [1.6 5.6];
                % initial wave parameters
                pI.exX0      = 0.7; % 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                    case '2D'
                        pTask.Y      = 2.75;
                        pTask.V_Q    = 1000 * 1.5; % 2D: 1, 1.5, 4
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Rinf = pTask.V_Q; pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX22s' % based on EX22m
                pTask.X      = 3;
                pTask.V_Linf = 0;
                pTask.V_x    = pTask.X/2;
                % initial wave parameters
                pI.exX0      = pTask.X/4;
                pI.exAlpha   = 1/120;
                pI.exKoef    = 30;
                switch sDimension
                    case '1D'
                        pTask.V_Rinf = -3000;
                    case '2D'
                        pTask.Y      = 2.8;
                        pTask.V_Rinf = 1500;
                        pTask.V_y    = [0 pTask.Y];
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Rinf ];
            case 'EX23' % delta function
                pTask.X      = 3;
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_Q    = 55;
                pTask.V_x0   = sqrt(3)-0.1;
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension % set values for dimentional parameters
                    case '2D'
                        pTask.Y      = 2.75;
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            case 'EX24' % potential step
                pTask.X = 16;
                c_e = 1.602; c_p = 1.054; c_m = 9.109;
                pTask.V_Linf = 0; 
                pTask.V_Rinf = 100*(2*c_m*c_e)/c_p^2;
                pTask.V_x    = pTask.X/2;
                % initial wave parameters
                pI.exX0      = pTask.X/4;
                sigma = 0.1; % 0.100722
                pI.exAlpha   = sigma^2/2;
                pI.exKoef    = 2*pi/sigma;
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        pTask.V_y    = [0 pTask.Y];
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Rinf ];            
            case 'EX25' % potential step by Sullivan - 2012, Example 1.2.2, Fig. 1.6 (Se1_1.m) with J=400, M=4500, KE=0.171, T=90 fs
                pTask.X = 40;
                c_e = 1.602; c_p = 1.054; c_m = 9.109;
                pTask.V_Linf = 0; % -0.75*0.1*(2*c_m*c_e)/c_p^2; 
                pTask.V_Rinf = -0.1*(2*c_m*c_e)/c_p^2; %  0.75*0.1*(2*c_m*c_e)/c_p^2;
                pTask.V_x    = pTask.X/2;
                % initial wave parameters
                pI.exX0      = pTask.X/4;
                pI.exAlpha   =(3/sqrt(2))^2/2;
                pI.exKoef    = 2*pi/3;
                switch sDimension
                    case '2D'
                        pTask.Y      = pTask.X;
                        pTask.V_y    = [0 pTask.Y];
                        % initial wave parameters
                        pI.exY0      = pTask.Y / 2;
                        pI.exAlphaY  = pI.exAlpha;
                        pI.exKoefY   = pI.exKoef;
                end
                pTask.V_Q = [ pTask.V_Linf pTask.V_Rinf ];            
            case 'EX30' % Poeshl-Teller
                pTask.X      = 4;
                % potential parameters
                pTask.V_Linf = 0; pTask.V_Rinf = pTask.V_Linf;
                pTask.V_x0   = 2.2; % 2.5; %
                pTask.V_alpha= 6;
                pTask.V_lambda= 9; % WALL: 4/9/17
                % initial wave parameters
                pI.exX0      = 1;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 10;
                switch sDimension
                    case '2D'
                        pTask.Y     = 2.75;
                        % initial wave parameters
                        pI.exY0     = pTask.Y / 2;
                        pI.exAlphaY = pI.exAlpha;
                        pI.exKoefY  = pI.exKoef;
                end
            case 'EX31' % Potential step
                pTask.X      = 4;
                % potential parameters
                pTask.V_Linf = 0;
                pTask.V_x0   = 2.3;
                pTask.V_alpha= 5;
                pTask.V_Q    = 1000 * 0.8; % 1D: 0.2, 0.8, 2
                % initial wave parameters
                pI.exX0      = 0.7;
                pI.exAlpha   = 1 / 120;
                pI.exKoef    = 30;
                switch sDimension
                    case '2D'
                        pTask.Y     = 2.75;
                        % initial wave parameters
                        pI.exY0     = pTask.Y / 2;
                        pI.exAlphaY = pI.exAlpha;
                        pI.exKoefY  = pI.exKoef;
                end
                pTask.V_Rinf = pTask.V_Q; pTask.V_Q = [ pTask.V_Linf pTask.V_Q pTask.V_Rinf ];
            otherwise
                logMessage('error', 'UnknownExampleName', 'Unknown example name.\n Current value is %s', pTask.sExampleName, mfilename);
        end
        switch sDimension
            case '1D'
                if pTask.B_inf <= 0
                    logMessage('error', 'WrongPhysicalValue', 'B_inf must be greather than zero.\n Current value is %d', pTask.B_inf, mfilename);
                end
            case '2D'
                if pTask.B1_inf <= 0 || pTask.B2_inf <= 0
                    logMessage('error', 'WrongPhysicalValue'...
                        , 'B1_inf and B2_inf must be greather than zero.\n Current values are %d and %d consequently'...
                        , [pTask.B1_inf pTask.B2_inf], mfilename);
                end
        end
        if pTask.hbar <= 0
            logMessage('error', 'WrongPhysicalValue', 'hbar must be greather than zero.\n Current value is %d', pTask.hbar, mfilename);
        end
        if pTask.rho_inf <= 0
            logMessage('error', 'WrongPhysicalValue', 'rho_inf must be greather than zero.\n Current value is %d', pTask.rho_inf, mfilename);
        end
        pTask.pI = pI;
    otherwise
        logMessage('error', 'UnknownTaskName', 'Unknown task name.\n Current value is %s', pTask.name, mfilename);
end