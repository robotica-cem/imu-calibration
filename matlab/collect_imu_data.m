function [t, acc, gyro] = collect_imu_data(N, plotdata, accCalib, gyrCalib)
% Will collect N data samples streamed by the Wireless IMU app
% If accCalib and gyrCalib are given, data are calibrated as they are
% received
%

% Kjartan Halvorsen
% 2019-11-20

%% Setup necessary infrastructure
  import('se.hendeby.sensordata.*');  % Used to receive data.


  try
    %% Create data link
    server = StreamSensorDataReader(3400);
    % Makes sure to resources are returned.
    sentinel = onCleanup(@() server.stop());

    server.start();  % Start data reception.
  catch e
    fprintf(['Unsuccessful connecting to client!\n' ...
      'Make sure to start streaming from the phone *after* '...
             'running this function!']);
    return;
  end


if nargin < 3
    accCalib.Dhat = eye(3);
    accCalib.bhat = zeros(3,1);
    gyrCalib.bhat = zeros(3,1);
elseif nargin < 4
    gyrCalib.bhat = zeros(3,1);
end
accDinv = inv(accCalib.Dhat);

gravAcc = 9.779; % Mexico city
%gravAcc = 9.818; % Stockholm

plotUpdatePeriod = 6; % Update plot with a period > 1 to avoid delays 

if nargin < 2
    plotdata = 0;
end

acc = nan(N,3);
gyro = nan(N,3);
t = nan(N,1);
accmagn = nan(1,N);

if plotdata
    acchandles = cell(3,1);
    gyrhandles = cell(3,1);
    h = figure('Position', [40,40,900,700]);

    subplot(321)
    acchandles{1} = plot(acc(:,1));
    set(acchandles{1}, 'Color', [0, 0.2, 0.7]);
    hold on 
    plot(1:N, repmat(gravAcc, 1, N), 'k:')
    plot(1:N, repmat(-gravAcc, 1, N), 'k:')
    title('Acc x')
    ylabel('m/s^2')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-20 20])
    
    subplot(323)
    acchandles{2} = plot(acc(:,2));
    set(acchandles{2}, 'Color', [0, 0.7, 0.7]);
    hold on
    plot(1:N, repmat(gravAcc, 1, N), 'k:')
    plot(1:N, repmat(-gravAcc, 1, N), 'k:')
    title('Acc y')
    ylabel('m/s^2')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-20 20])
    
    subplot(325)
    acchandles{3} = plot(acc(:,3));
    set(acchandles{3}, 'Color', [0.6, 0, 0]);
    hold on
    plot(1:N, repmat(gravAcc, 1, N), 'k:')
    plot(1:N, repmat(-gravAcc, 1, N), 'k:')
    title('Acc z')
    ylabel('m/s^2')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-20 20])

    subplot(322)
    gyrhandles{1} = plot(gyro(:,1));
    set(gyrhandles{1}, 'Color', [0, 0.2, 0.7]);
    title('Gyro x')
    ylabel('rad/s')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-4 4])
    
    subplot(324)
    gyrhandles{2} = plot(gyro(:,2));
    set(gyrhandles{2}, 'Color', [0, 0.7, 0.7]);
    title('Gyro y')
    ylabel('rad/s')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-4 4])

    subplot(326)
    gyrhandles{3} = plot(gyro(:,3));
    set(gyrhandles{3}, 'Color', [0.6, 0, 0]);
    title('Gyro z')
    ylabel('rad/s')
    set(gca, 'XLim', [0 N])
    set(gca, 'YLim', [-4 4])

    h2 = figure('Position', [940,700,500,300]);
    magnhandle = plot(accmagn, 'Color', [0.6, 0, 0]);
    hold on 
    plot(1:N, repmat(gravAcc, 1, N), 'k:');
    title('Magnitude of acceleration vector')
    ylabel('m/s^2')
    set(gca, 'XLim', [0 N])
    set(gca, 'Ylim', [0.5 1.5]*gravAcc)
    th = text(N/8, 1.3*gravAcc, ' ');
    
    drawnow
end
    
for k=1:N

    
    data = server.getNext(10);

    t_k = data(1,1);
    acc_k = data(1, 2:4)';
    gyro_k = data(1, 5:7)';
    if ~isnan(t_k)
        acc(k,:) = accDinv*(acc_k - accCalib.bhat);
        gyro(k,:) = gyro_k - gyrCalib.bhat;
        t(k) = t_k;
        accmagn(k) = sqrt(sum(acc(k,:).^2));
        if plotdata
            % Update plot only every mth sample
            if ~mod(k,plotUpdatePeriod)
                start(timer('StartDelay', 0.004, 'TimerFcn', @(~,~)update_plot(acchandles, gyrhandles, acc, gyro, magnhandle, accmagn, th, k)));
            end
        end
    end
end
t = t-t(1);
end

function update_plot(acchandles, gyrhandles, acc, gyro, magnhandle, accmagn, th, k)
    % Function meant for asynchrounous execution
    for i=1:3
        set(acchandles{i}, 'YData', acc(:,i));
        set(gyrhandles{i}, 'YData', gyro(:,i));
    end
    set(magnhandle, 'YData', accmagn);
    set(th, 'String', sprintf('%5.3f', accmagn(k)))
    
    drawnow update
end

function am = accMagn(acc)
    am = sqrt(sum(acc.^2));
end