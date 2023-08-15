function [AllTrajF,AllTrajL] = Trajectory(acc, kPer, height, radius, dt, kTraj)

AllTrajF = []; AllTrajL = [];
for i = 1:kPer
    traF = [radius 0 0; 
            0 radius height;
            -radius 0 2*height;
            0 -radius 3*height;
            radius 0 4*height];

    traF(:,3) = traF(:,3)+4*height*(i-1);
    
    traL = [kTraj*radius 0 0; 
            0 -kTraj*radius height;
            -kTraj*radius 0 2*height;
            0 kTraj*radius 3*height;
            kTraj*radius 0 4*height];

    traL(:,3) = traL(:,3)+4*height*(i-1)+2*height;
    eulL = -pi/2*[0 1 2 3 0];
    axixF = [1 0 0; 0 1 0; 1 1 1; 1 0 0; 0 1 0; 1 1 1];
    flag = 0;
    for j = 1:4
        Rf0 = axang2tform([0 0 1 pi/2])*axang2tform([axixF(i,:) pi/2*(j-1)]); 
%         Rf0 = axang2tform([0 0 1 pi/2])*axang2tform([axixF pi/2*(j-1)]); 
        tf0=trvec2tform(traF(j,:)); 
        Rl0 = axang2tform([0 0 1 pi/2])*axang2tform([0 0 1 eulL(j)])*axang2tform([0 1 0 eulL(j)])*axang2tform([1 0 0 flag*pi/6]); tl0=trvec2tform(traL(j,:)); 
        Tf0 = tf0*Rf0; 
        Tl0 = tl0*Rl0;
        
        flag = -flag;
        Rf1 = axang2tform([0 0 1 pi/2])*axang2tform([axixF(i,:) 1*pi/2*j]); 
%         Rf1 = axang2tform([0 0 1 pi/2])*axang2tform([axixF 1*pi/2*j]); 
        tf1=trvec2tform(traF(j+1,:));
        Rl1 = axang2tform([0 0 1 pi/2])*axang2tform([0 0 1 eulL(j+1)])*axang2tform([0 1 0 eulL(j+1)])*axang2tform([1 0 0 flag*pi/6]); tl1=trvec2tform(traL(j+1,:)); 
        Tf1=tf1*Rf1; 
        Tl1=tl1*Rl1;

        s = norm(Tf1(1:3,4)-Tf0(1:3,4));
        if s~=0
            n = round(3*sqrt(s/2/acc)/dt);
        else
            n = 100;
        end
        
        AllTrajF = cat(3,AllTrajF,ctraj(Tf0,Tf1,n));
        AllTrajL = cat(3,AllTrajL,ctraj(Tl0,Tl1,n));
    end
end
end

