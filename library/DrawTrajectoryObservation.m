function DrawTrajectoryObservation(transF, transL, AllTrajF, AllTrajL, imgSample, kSam, nFace, bPixel, plotTag, bOneStep)

    figure
    Axis = 2*[0.1 0 0; 0 0.1 0; 0 0 0.1];    
    plotTrajF = transF(imgSample,:);
    plotTrajL = transL(imgSample,:);
    nSamFea = length(imgSample);
    for i=1:nSamFea
        tic
        plot3(plotTrajF(1:i,1),plotTrajF(1:i,2),plotTrajF(1:i,3),'b','linewidth',1)
        grid on, hold on, axis equal
        axis([-1.3 1.3 -1.3 1.3 -0.5 3])
        plot3(plotTrajL(1:i,1),plotTrajL(1:i,2),plotTrajL(1:i,3),'r','linewidth',1)

        Tf = AllTrajF(:,:,(i-1)*kSam+1);
        temp = Tf(1:3,:)*[Axis; ones(1,3)];
        plot3([Tf(1,4) temp(1,1)],[Tf(2,4) temp(2,1)],[Tf(3,4) temp(3,1)],'b','linewidth',1)
        plot3([Tf(1,4) temp(1,2)],[Tf(2,4) temp(2,2)],[Tf(3,4) temp(3,2)],'r','linewidth',1)
        plot3([Tf(1,4) temp(1,3)],[Tf(2,4) temp(2,3)],[Tf(3,4) temp(3,3)],'g','linewidth',1)
        Tl = AllTrajL(:,:,(i-1)*kSam+1);
        temp = Tl(1:3,:)*[Axis; ones(1,3)];
        plot3([Tl(1,4) temp(1,1)],[Tl(2,4) temp(2,1)],[Tl(3,4) temp(3,1)],'b','linewidth',1)
        plot3([Tl(1,4) temp(1,2)],[Tl(2,4) temp(2,2)],[Tl(3,4) temp(3,2)],'r','linewidth',1)
        plot3([Tl(1,4) temp(1,3)],[Tl(2,4) temp(2,3)],[Tl(3,4) temp(3,3)],'g','linewidth',1)

        for j=1:nFace
            if bPixel(i,j)==1
                temp = reshape(plotTag(i,(j-1)*3*4+1:j*3*4),3,4);
                p4 = plot3(temp(1,:), temp(2,:), temp(3,:),'mp','linewidth',1);
            end
        end
        hold off
        drawnow
        drwaTime = toc;
        if bOneStep
            pause;
        end
    end
end