function h = plot_acc_gyro(acc, gyr, fig)

    if nargin < 3
        h = figure('Position', [40,40,900,700]);
        xcol = [0, 0.2, 0.8];
        ycol = [0, 0.8, 0.8];
        zcol = [0.8, 0, 0];
        axh.accxh = subplot(321);
        title('Acc x')
        ylabel('m/s^2')
        axh.accyh = subplot(323);
        title('Acc y')
        ylabel('m/s^2')
        axh.acczh = subplot(325);
        title('Acc z')
        ylabel('m/s^2')
        axh.gyrxh = subplot(322);
        title('Gyro x')
        ylabel('rad/s')
        axh.gyryh = subplot(324);
        title('Gyro y')
        ylabel('rad/s')
        axh.gyrzh = subplot(326);
        title('Gyro z')
        ylabel('rad/s')
        set(h, 'UserData', axh);
        %set(get(h, 'Children'), 'NextPlot', 'add');
    else
        % Assume plots already exists. Add new plots
        axh = get(fig, 'UserData');
        
        if isempty(axh) | ~isstruct(axh) | length(fieldnames(axh)) ~= 6
            warning('Expected to find 6 subplots. Creating new figure instead.')
            h = plot_acc_gyro(acc,gyr);
            return
        end
        h = fig;
        set(fig,'Position', [40,40,900,700]);
        xcol = [0, 0.1, 0.5];
        ycol = [0, 0.5, 0.5];
        zcol = [0.5, 0.1, 0];
    end
    
    hold(axh.accxh, 'on')
    plot(axh.accxh, acc.t, acc.dta(1,:), 'Color', xcol);

    hold(axh.accyh, 'on')
    plot(axh.accyh, acc.t, acc.dta(2,:), 'Color', ycol);
    
    hold(axh.acczh, 'on')
    plot(axh.acczh, acc.t, acc.dta(3,:), 'Color', zcol);

    hold(axh.gyrxh, 'on')
    plot(axh.gyrxh, gyr.t, gyr.dta(1,:), 'Color', xcol);
    
    hold(axh.gyryh, 'on')
    plot(axh.gyryh, gyr.t, gyr.dta(2,:), 'Color', ycol);

    hold(axh.gyrzh, 'on')
    plot(axh.gyrzh, gyr.t, gyr.dta(3,:), 'Color', zcol);
end
