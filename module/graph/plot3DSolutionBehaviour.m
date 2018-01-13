function OutParams = plot3DSolutionBehaviour(U, Params, pTask)
%%
[X, Y] = meshgrid(Params.nMeshGridX, Params.nMeshGridY);
if Params.nDenominator == 0
    Params.nDenominator = 1;
end
%% Choose style
switch Params.sStyle
    case 'phd-2012'
        I = 2;
        ax = NaN * ones(1, I);
        for i=1:I
            ax(i) = subplot(1, I, i ...
                , 'FontSize'       , Params.nFontSize   ...
                , 'FontWeight'     , Params.sFontWeight ...
                , 'NextPlot'       , Params.sNextPlot   ...
                , 'Projection'     , 'Orthographic'     ...% 'Perspective'      ...%, 'CameraPosition' , Params.nCameraPosition
                );
            xlim(ax(i), Params.nLimitX); % xlabel(ax(i), Params.sLabelX); % , 'HorizontalAlignment', 'Right'
            ylim(ax(i), Params.nLimitY); % ylabel(ax(i), Params.sLabelY); % , 'HorizontalAlignment', 'Left'
            zlim(ax(i), [-1 1]); % zlabel(Params.sLabelZ);
            if I > 1
                switch i
                    case 1
                        W = abs(U);
                        clim = [-max(max(W)) max(max(W))];
                    case 2
                        W = real(U); % imag(U);
                end
            else
                switch Params.sType
                    case 'modulus'
                        W = abs(U);
                    case 'real'
                        W = real(U);
                    case 'imaginary'
                        W = imag(U);
                end
                clim = [-max(max(W)) max(max(W))];
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
                    P = zeros(length(Params.nMeshGridY), length(Params.nMeshGridX)); 
                    P(ind_y, ind_x) = 1; % P(ind_y, ind_x) = max(clim(2),eps); % P(:, ind_x) = 1; % 
                    surf(X, Y, P, 'Parent', ax(i), 'FaceColor', [0 1 0], 'FaceAlpha', 0.2, 'EdgeAlpha', 0); % 'FaceAlpha', 0.1
                    set(ax(i), 'NextPlot', 'ReplaceChildren')
                end
            end
            % Create light
            light('Parent', ax(i), 'Style', 'Local'); % 'Position', [4 -12 14]
            box(ax(i), 'on'); grid(ax(i), 'on'); % hold(ax(i), 'off');
            view(ax(i), [15 20]); % [60 10]); % [40 30]); % 
        end
        %% Update coloring
        if sum(abs(clim))>0
            caxis(clim); set(ax, 'CLim', clim);
        end
    case 'phd-2012-d'
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
            zlim(ax(i), [-1 1]); % zlabel(Params.sLabelZ);
            switch i
                case 1
                    W = abs(U);
                    clim = [-max(max(W)) max(max(W))];
                case 2
                    W = real(U); % imag(U);
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
                    P = zeros(length(Params.nMeshGridY), length(Params.nMeshGridX)); 
                    P(ind_y, ind_x) = 1; % P(ind_y, ind_x) = max(clim(2),eps); % P(:, ind_x) = 1; % 
                    surf(X, Y, P, 'Parent', ax(i), 'FaceColor', [0 1 0], 'FaceAlpha', 0.10, 'EdgeAlpha', 0); % 'FaceAlpha', 0.15
                    set(ax(i), 'NextPlot', 'ReplaceChildren')
                end
            end
            % Create light
            light('Parent', ax(i), 'Style', 'Local'); % 'Position', [4 -12 14]
            box(ax(i), 'on'); grid(ax(i), 'on'); % hold(ax(i), 'off');
            view(ax(i), [15 20]); % [60 10]); % [40 30]); % 
        end
        %% Update coloring
        if sum(abs(clim))>0
            caxis(clim); set(ax, 'CLim', clim);
        end
        % colorbar('FontWeight', Params.sFontWeight, 'FontSize', Params.nFontSize, 'Location', 'SouthOutside');
    case 'vestnik-mpei-2010'
        subplot(2, 2, 1);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontName'       , 'Times New Roman'...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left' , 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        surf(X, Y, abs(U), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(0, 90);
        grid('off');
        clear zto zt ztm ztl yto yt ytm ztl xt xtl xtm;
        yto = get(gca, 'YTick');
        if yto(1) == 0
            yt = yto(2:end);
        else
            yt = yto;
        end
        ytl = num2str(yt', '%1.1f'); % ytl = num2str(yt', '%3.2f'); 1D
        for i=1:length(yt)
            ytm(i,:) = strrep(ytl(i,:), '.', ',');
        end
        set(gca, 'YTick', yt, 'YTickLabel', ytm);
        xto  = get(gca, 'XTick');
        if xto(1) == 0
            xt = xto(2:end);
        else
            xt = xto;
        end
        xtl = num2str(xt', '%1.1f');
        for i=1:length(xt)
            xtm(i,:) = strrep(xtl(i,:), '.', ',');
        end
        set(gca, 'XTick', xto, 'XTickLabel', ['0  ' ; xtm ]);
        set(get(gca, 'XLabel'), 'pos', [2.8 -0.1]);                 % set(get(gca, 'XLabel'), 'pos', [2.8 -0.0015]);
        set(get(gca, 'YLabel'), 'Rotation', 0.0, 'pos', [-0.1 2.7]);% set(get(gca, 'YLabel'), 'Rotation', 0.0, 'pos', [-0.1 0.0425]);
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask),length(Params.nMeshGridY),1);
            contourf(X, Y, P/Params.nDenominator, 10);
        end
        hc = colorbar('FontWeight', 'bold', 'FontSize', 14);
        zto = get(hc, 'YTick');
        if zto(1) == 0
            zt = zto(2:end);
        else
            zt = zto;
        end
        ztl = num2str(zt', '%1.1f');
        for i=1:length(zt)
            ztm(i,:) = strrep(ztl(i,:), '.', ',');
        end
        set(hc, 'YTick', zto, 'YTickLabel', ['0  ' ; ztm ]);
        %
        subplot(2, 2, 2);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontName'       , 'Times New Roman'...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left' , 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        zlim([0 1.2])
        surf(X, Y, abs(U), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        view([ 5 20 ]);
        grid('on');
        clear zto zt ztm ztl yto yt ytm ztl xt xtl xtm;
        zto = get(gca, 'ZTick');
        if zto(1) == 0
            zt = zto(2:end);
        else
            zt = zto;
        end
        ztl = num2str(zt', '%1.1f');
        for i=1:length(zt)
            ztm(i,:) = strrep(ztl(i,:), '.', ',');
        end
        set(gca, 'ZTick', zt, 'ZTickLabel', ztm);
        yto = get(gca, 'YTick');
        if yto(1) == 0
            yt = yto(2:end);
        else
            yt = yto;
        end
        ytl = num2str(yt', '%1.0f'); % ytl = num2str(yt', '%3.2f'); 1D
        for i=1:length(yt)
            ytm(i,:) = strrep(ytl(i,:), '.', ',');
        end
        set(gca, 'YTick', yt, 'YTickLabel', ytm);
        xto = get(gca, 'XTick');
        if xto(1) == 0
            xt = xto(2:end);
        else
            xt = xto;
        end
        xtl = num2str(xt', '%1.1f');
        for i=1:length(xt)
            xtm(i,:) = strrep(xtl(i,:), '.', ',');
        end
        set(gca, 'XTick', xto, 'XTickLabel', ['0  ' ; xtm ]);
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask),length(Params.nMeshGridY), 1);
            surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat');
        end
        %
        subplot(2, 2, 3);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontName'       , 'Times New Roman'...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        set(get(gca, 'YLabel'), 'Rotation', 0.0, 'pos', [0 0.045]);
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left' , 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        surf(X, Y, real(U), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(0, 90);
        grid('off');
        clear zto zt ztm yto yt ytm ytl xt xtl xtm;
        yto = get(gca, 'YTick');
        if yto(1) == 0
            yt = yto(2:end);
        else
            yt = yto;
        end
        ytl = num2str(yt', '%1.1f'); % ytl = num2str(yt', '%3.2f'); 1D
        for i=1:length(yt)
            ytm(i,:) = strrep(ytl(i,:), '.', ',');
        end
        set(gca, 'YTick', yt, 'YTickLabel', ytm);
        xto = get(gca, 'XTick');
        if xto(1) == 0
            xt = xto(2:end);
        else
            xt = xto;
        end
        xtl = num2str(xt', '%1.1f');
        for i=1:length(xt)
            xtm(i,:) = strrep(xtl(i,:), '.', ',');
        end
        set(gca, 'XTick', xto, 'XTickLabel', ['0  ' ; xtm ]);
        set(get(gca, 'XLabel'), 'pos', [2.8 -0.1]);                 % set(get(gca, 'XLabel'), 'pos', [2.8 -0.0015]);
        set(get(gca, 'YLabel'), 'Rotation', 0.0, 'pos', [-0.1 2.7]);% set(get(gca, 'YLabel'), 'Rotation', 0.0, 'pos', [-0.1 0.0425]);
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            contourf(X, Y, P/Params.nDenominator, 10);
        end
        hc = colorbar('FontWeight', 'bold', 'FontSize', 14);
        zto = get(hc, 'YTick');
        if zto(1) == 0
            zt = zto(2:end);
        else
            zt = zto;
        end
        ztl = num2str(zt', '%1.1f');
        for i=1:length(zt)
            ztm(i,:) = strrep(ztl(i,:), '.', ',');
        end
        if sum(zt<=0)>0
            set(hc, 'YTick', zto, 'YTickLabel', ztm);
        else
            set(hc, 'YTick', zto, 'YTickLabel', ['0  ' ; ztm ]);
        end
        %
        subplot(2, 2, 4);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontName'       , 'Times New Roman'...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left' , 'FontName', 'Times New Roman', 'FontAngle', 'italic');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        zlim([-1.2 1.2])
        surf(X, Y, real(U), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        view([ 5 20 ]);
        grid('on');
        yt = get(gca, 'YTick');
        clear zto zt ztm ztl yto yt ytm ztl xt xtl xtm;
        zto = get(gca, 'ZTick');
        if zto(1) == 0
            zt = zto(2:end);
        else
            zt = zto;
        end
        ztl = num2str(zt', '%1.1f');
        for i=1:length(zt)
            ztm(i,:) = strrep(ztl(i,:), '.', ',');
        end
        set(gca, 'ZTick', zt, 'ZTickLabel', ztm);
        yto = get(gca, 'YTick');
        if yto(1) == 0
            yt = yto(2:end);
        else
            yt = yto;
        end
        ytl = num2str(yt', '%1.0f'); % ytl = num2str(yt', '%3.2f'); 1D
        for i=1:length(yt)
            ytm(i,:) = strrep(ytl(i,:), '.', ',');
        end
        set(gca, 'YTick', yt, 'YTickLabel', ytm);
        xto = get(gca, 'XTick');
        if xto(1) == 0
            xt = xto(2:end);
        else
            xt = xto;
        end
        xtl = num2str(xt', '%1.1f');
        for i=1:length(xt)
            xtm(i,:) = strrep(xtl(i,:), '.', ',');
        end
        set(gca, 'XTick', xto, 'XTickLabel', ['0  ' ; xtm ]);
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            hp = surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat');
        end
        Params.hAxes = gca;
    case 'master-2010'
        subplot(2, 2, 1);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        % contourf(Params.nMeshGridX, Params.nMeshGridY, abs(U), 'FaceColor', 'interp');
        surf(X, Y, abs(U), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(0, 90);
        grid('off');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            contourf(X, Y, P/Params.nDenominator, 10);
        end
        colorbar('FontWeight', 'bold', 'FontSize', 14)
        %
        subplot(2, 2, 2);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        surf(X, Y, abs(U), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        % zlabel(Params.sLabelZ);
        view([ 5 20 ]);
        grid('on');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat');
        end
        %
        subplot(2, 2, 3);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        % contourf(Params.nMeshGridX, Params.nMeshGridY, abs(U), 'FaceColor', 'interp');
        surf(X, Y, real(U), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(0, 90);
        grid('off');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            contourf(X, Y, P/Params.nDenominator, 10);
        end
        colorbar('FontWeight', 'bold', 'FontSize', 14)
        %
        subplot(2, 2, 4);
        set(gca, 'CameraPosition' , Params.nCameraPosition...
            , 'FontSize'       , Params.nFontSize...
            , 'FontWeight'     , Params.sFontWeight...
            , 'nextplot'       , Params.sNextPlot...
            );
        xlabel(Params.sLabelX, 'HorizontalAlignment', 'Right');
        ylabel(Params.sLabelY, 'HorizontalAlignment', 'Left');
        xlim(Params.nLimitX);
        ylim(Params.nLimitY);
        surf(X, Y, real(U), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        view([ 5 20 ]);
        grid('on');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask), length(Params.nMeshGridY), 1);
            hp = surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat');
            % alpha(h, 0.5)
        end
        Params.hAxes = gca;
    case 'report-2010'
        [XI, YI] = meshgrid(Params.nMeshGridXI, Params.nMeshGridYI);
        switch Params.sType
            case 'modulus'
                W = abs(U);
            case 'real'
                W = real(U);
            case 'imaginary'
                W = imag(U);
        end
        WI = interp2(X, Y, W, XI, YI);
        surf(XI, YI, WI, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        zlabel(Params.sLabelZ);
        view([ 40 20 ]);
        grid('on');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask),length(Params.nMeshGridY),1);
            hp = surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat'); % alpha(h, 0.5)
            % meshz(X, Y, P/Params.nDenominator);
            set(gca, 'nextplot', Params.sNextPlot)
        end
    case 'krm-2009'
        switch Params.sType
            case 'modulus'
                W = abs(U);
            case 'real'
                W = real(U);
            case 'imaginary'
                W = imag(U);
        end
        meshc(X, Y, W);
        zlabel(Params.sLabelZ);
        view([ 45 30 ]);
        grid('off');
        if Params.bPlotPotential % Plot potential function
            set(gca, 'nextplot', 'add');
            P = repmat(V(Params.nMeshGridX, pTask),length(Params.nMeshGridY),1);
            hp = surf(X, Y, P/Params.nDenominator, 'FaceColor', 'texturemap', 'EdgeColor', 'flat'); % alpha(h, 0.5)
            % meshz(X, Y, P/Params.nDenominator);
            set(gca, 'nextplot', Params.sNextPlot)
        end
end

%% Other styles
% surf(X, Y, W);
% surf(X, Y, W, 'CData', abs(W), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'EdgeColor', 'none');
% surf(X, Y, W, 'CDataMapping', 'scaled', 'EdgeColor', 'none');
% surf(X, Y, W, 'CDataMapping', 'direct');
% surf(X, Y, W, 'CData', abs(W), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceColor', 'interp');
% surf(X, Y, W, 'CData', abs(W), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'EdgeColor', [0 0 1]);
% contourf(W, 10), colorbar;
% view(Params.nViewAngle);

if Params.bPlayPause % Play pause
    keyboard % pause(1*0.01)
end

%% Movie control
Params.pMovie.nFrameCount = Params.pMovie.nFrameCount + 1; Params.pMovie.pFrame(Params.pMovie.nFrameCount) = getframe(Params.hFigure);

%% Set output values
OutParams = Params;