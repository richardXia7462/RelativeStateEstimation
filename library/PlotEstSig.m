function [p1, p2, p3] = PlotEstSig(t,errXm,var,ylab,tit,type,sig3type,varargin)
% t=tLen;true=zeros(nS,3);est=erre;var=p_cov(:,1:3)*(180/pi)^2;

    if size(var,1)==1
        var=repmat(var,size(errXm,1),1);
    end
    sig3=var.^(0.5)*3;
    
    if nargin==8
        if varargin{1}==1
            figure
        end 
    elseif nargin==7
        figure
    end
    subplot(311), hold on
    plot(t,sig3(:,1),sig3type,'linewidth',2); plot(t,-sig3(:,1),sig3type,'linewidth',2);
    p1=plot(t,errXm(:,1),type,'linewidth',2);
%     plot(t,sig3(:,1),sig3type,t,errXm(:,1),type,t,-sig3(:,1),sig3type)
    set(gca,'fontsize',20), grid on; 
    legend(p1,ylab(1))
%     ylabel(ylab(1))
    % title(tit)
    
    subplot(312), hold on
    plot(t,sig3(:,2),sig3type,'linewidth',2); plot(t,-sig3(:,2),sig3type,'linewidth',2);
    p2=plot(t,errXm(:,2),type,'linewidth',2);
%     plot(t,sig3(:,2),sig3type,t,errXm(:,2),type,t,-sig3(:,2),sig3type)
    set(gca,'fontsize',20), grid on
    legend(p2,ylab(2))
%     ylabel(ylab(2))
    
    subplot(313), hold on
    plot(t,sig3(:,3),sig3type,'linewidth',2); plot(t,-sig3(:,3),sig3type,'linewidth',2);
    p3=plot(t,errXm(:,3),type,'linewidth',2);
%     plot(t,sig3(:,3),sig3type,t,errXm(:,3),type,t,-sig3(:,3),sig3type)
    set(gca,'fontsize',20), xlabel('Time (sec)'), grid on
    legend(p3,ylab(3))
%     ylabel(ylab(3))
end