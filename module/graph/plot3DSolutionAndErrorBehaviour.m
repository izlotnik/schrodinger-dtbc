function OutParams = plot3DSolutionAndErrorBehaviour(U, E,  Params, pTask)
%%
[X, Y] = meshgrid(Params.nMeshGridX, Params.nMeshGridY);
%% Choose style
switch Params.sStyle
    case 'phd-2012'
        ax = NaN * ones(1, 2);
        for i=1:2
            ax(i) = subplot(1, 2, i ...
                , 'FontSize'       , Params.nFontSize   ...
                , 'FontWeight'     , Params.sFontWeight ...
                , 'NextPlot'       , Params.sNextPlot   ...
                , 'Projection'     , 'Orthographic'     ...% 'Perspective'      ...%, 'CameraPosition' , Params.nCameraPosition
                );
            xlim(ax(i), Params.nLimitX); % xlabel(ax(i), Params.sLabelX); % , 'HorizontalAlignment', 'Right'
            ylim(ax(i), Params.nLimitY); % ylabel(ax(i), Params.sLabelY); % , 'HorizontalAlignment', 'Left'
            % zlim(ax(i), [-1 1]); % zlabel(Params.sLabelZ);
            switch i
                case 1
                    W = abs(U);
                    % clim = [-max(max(W)) max(max(W))];
                case 2
                    W = abs(E); % imag(U);
            end
            % Create surf % 'FaceColor':'texturemap'
            surfc(X, Y, W, 'Parent', ax(i), 'FaceLighting', 'Gouraud', 'FaceColor', 'Interp', 'EdgeColor', 'None');
            % contourf(Params.nMeshGridX, Params.nMeshGridY, abs(U), 'FaceColor', 'interp');
            % shading flat % shading interp % shading faceted %
            % axis tight; daspect([1 1 1]); camproj perspective; camlight; lighting gouraud; alpha(0.75); rotate3d on;
            if Params.bPlotPotential % Plot potential \chi_I(x)*\chi_J(y)
                if isfield(pTask, 'V_x')
                    ind_x = pTask.V_x(1) <= Params.nMeshGridX & Params.nMeshGridX <= pTask.V_x(2);
                else
                    ind_x = [];
                end
                if isfield(pTask, 'V_y')
                    ind_y = pTask.V_y(1) <= Params.nMeshGridY & Params.nMeshGridY <= pTask.V_y(2);
                else
                    ind_y = [];
                end
                if length(ind_x) + length(ind_y) > 0
                    set(ax(i), 'NextPlot', 'Add');
                    P = zeros(length(Params.nMeshGridY), length(Params.nMeshGridX)); P(ind_y, ind_x) = max( max(max(W)), eps);
                    surf(X, Y, P, 'Parent', ax(i), 'FaceColor', [0 1 0], 'FaceAlpha', 0.10, 'EdgeAlpha', 0); % 'FaceAlpha', 0.15
                    set(ax(i), 'NextPlot', 'ReplaceChildren')
                end
            end
            % Create light
            light('Parent', ax(i), 'Style', 'Local'); % 'Position', [4 -12 14]
            box(ax(i), 'on'); grid(ax(i), 'on'); % hold(ax(i), 'off');
            view(2); % view(ax(i), [15 20]); % [65 10]); % [40 30]); %
            colorbar
        end
end

if Params.bPlayPause % Play pause
    keyboard % pause(1*0.01)
end

%% Movie control
Params.pMovie.nFrameCount = Params.pMovie.nFrameCount + 1; Params.pMovie.pFrame(Params.pMovie.nFrameCount) = getframe(Params.hFigure);

%% Set output values
OutParams = Params;