function PlotTwoBiases(t, gt, est,var)
    figure
    ylab={'x';'y';'z'};
    sig3=var.^(0.5)*3;
    for i=1:3
        subplot(3,2,2*(i-1)+1),hold on, grid on
        plot(t, est(:,i),'linewidth',3)
%         plot(t, gt(:,i),'linewidth',2)
        plot(t, gt(:,i)+sig3(:,i),'r:','linewidth',3)
        plot(t, gt(:,i)-sig3(:,i),'r:','linewidth',3)
        if i==1
            title('Follower Gyro Bias (rad/s)')
        end
        if i==3
            xlabel('Time (s)')
        end
        ylabel(ylab{i})
        set(gca,'fontsize',30), grid on; 
    end
    for i=1:3
        subplot(3,2,2*i),hold on, grid on
        plot(t, est(:,i+3),'linewidth',3)
%         plot(t, gt(:,i+3),'linewidth',2)
        plot(t, gt(:,i+3)+sig3(:,i+3),'r:','linewidth',3)
        plot(t, gt(:,i+3)-sig3(:,i+3),'r:','linewidth',3)
        if i==1
            title('Leader Gyro Bias (rad/s)')
        end
%         if i==2
%             legend('Groundtruth', 'Estimation')
%         end
        if i==3
            xlabel('Time (s)')
        end
        ylabel(ylab{i})
        set(gca,'fontsize',30), grid on; 
    end
end