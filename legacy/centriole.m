function [success overlap b_deg] = centriole(SYM,LK,CW_RADIUS,r,L,MTn,MT_RADIUS,GAMMA,CORRECT_OVERLAP,show_result)

% PARAMETERS --------------------------------------------------------------
    if nargin < 1
       clear all; clc; close all;
       SYM = 9; 
       CW_RADIUS = SYM*7.7/(2*pi); % central hub radius
       r = 27; % SAS-6 coiled coil length
       L_single = 20; % triplet length (from the end of A to the end of C)
       L_ratio = [1 1.75 2.65];
       MTn = 1; % number of MT (singlet, doublet,triplet)
       MT_RADIUS = 10; % Microtubule radius 
       LK = 28.5; % A-C linker length
       GAMMA = 60; % angle for LK anchoring
       CORRECT_OVERLAP = 0;
       show_result = 1;
       
       L = L_single*L_ratio(MTn); % real length
    end
%--------------------------------------------------------------------------

    disp(' ');
    disp(['-------------- desired LK: ',num2str(LK)]);
    
    if show_result
        figure;
        % MT colors
        cmap = colormap(gray);
        idxi = linspace(1,size(cmap,1),ceil(sqrt(SYM)));
        idxi = round(idxi);
        col = cmap(idxi,:);
    end
   
    t = linspace(0,2*pi,SYM+1);
    tt = linspace(0,2*pi,100);
    a = 2*pi/SYM;
    R = (r+CW_RADIUS-MT_RADIUS)/cos(a/2);
    g = pi*GAMMA/180; % convert GAMMA to radians
    
    bb = linspace(0,pi/2,181*2); % 0.5? precision

    % calculate LK and angle b (with respect to tangent)
    i = 1;
    for j = 1:length(bb) 
        % anchor points on A and C MT
        if MTn > 1 % doublet-triplet
            anchC(j,:) = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i))-(L-2*MT_RADIUS)*cos(pi/2-t(i)+bb(j)) - MT_RADIUS*cos(pi/2-t(i)+bb(j)-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i))+(L-2*MT_RADIUS)*sin(pi/2-t(i)+bb(j)) + MT_RADIUS*sin(pi/2-t(i)+bb(j)-g)];
            anchA(j,:) = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)) + MT_RADIUS*cos(pi/2-t(i+1)+bb(j)-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i+1)) - MT_RADIUS*sin(pi/2-t(i+1)+bb(j)-g)];
        else % singlet
            anchC(j,:) = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i)) - MT_RADIUS*cos(pi/2-t(i)-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i)) + MT_RADIUS*sin(pi/2-t(i)-g)];
            anchA(j,:) = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)) + MT_RADIUS*cos(pi/2-t(i+1)-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i+1)) - MT_RADIUS*sin(pi/2-t(i+1)-g)];
        end
        
        % distance between the two anchor points (i.e. LK)
        d = [anchC(j,:);anchA(j,:)];
        LKmV(j) = sqrt(sum((d(1,:)-d(2,:)).^2));
    end
    idx = find(abs(LKmV-LK)<0.15);
    if ~isempty(idx)
        didx = find(diff(idx)>1);
        if ~isempty(didx)
            idx = idx(didx(end)+1:end);
        end
        [mini idxi] = min(abs(LKmV(idx)-LK));
        idx = idx(idxi);
        success = 1;
    else
        [mini idx] = min(abs(LKmV-LK));
        success = 0;
        disp('It is not possible to reach the desired A-C LINKER.');
    end
    LKm = round2xdigit(LKmV(idx));
    b = bb(idx);
    b_deg = round2xdigit(b*180/pi);
    
    anchC = anchC(idx,:);
    anchA = anchA(idx,:);
    
    % check if overlap
    overlap = 0; % true, unless proven it is not
    if success
        r1 = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i)),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))];
        r2 = [(r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)),(r+CW_RADIUS+MT_RADIUS)*sin(t(i+1))];
        d2 = sqrt( (r1(1)-r2(1))^2 + (r1(2)-r2(2))^2 );
        if d2 < 2*MT_RADIUS
           disp('It is not possible to avoid MT overlap.. Increase SAS-6 length (r)');
           overlap = 1;
        else
            x0V = [(r+CW_RADIUS+MT_RADIUS)*cos(t(1)) (r+CW_RADIUS+MT_RADIUS)*cos(t(1))-(L-2*MT_RADIUS)*cos(pi/2-t(1)+b) (r+CW_RADIUS+MT_RADIUS)*cos(t(1))-0.5*(L-2*MT_RADIUS)*cos(pi/2-t(1)+b)];
            y0V = [(r+CW_RADIUS+MT_RADIUS)*sin(t(1)) (r+CW_RADIUS+MT_RADIUS)*sin(t(1))+(L-2*MT_RADIUS)*sin(pi/2-t(1)+b) (r+CW_RADIUS+MT_RADIUS)*sin(t(1))+0.5*(L-2*MT_RADIUS)*sin(pi/2-t(1)+b)];
            x0pV = [(r+CW_RADIUS+MT_RADIUS)*cos(t(2)) (r+CW_RADIUS+MT_RADIUS)*cos(t(2))-(L-2*MT_RADIUS)*cos(pi/2-t(2)+b) (r+CW_RADIUS+MT_RADIUS)*cos(t(2))-0.5*(L-2*MT_RADIUS)*cos(pi/2-t(2)+b)];
            y0pV = [(r+CW_RADIUS+MT_RADIUS)*sin(t(2)) (r+CW_RADIUS+MT_RADIUS)*sin(t(2))+(L-2*MT_RADIUS)*sin(pi/2-t(2)+b) (r+CW_RADIUS+MT_RADIUS)*sin(t(2))+0.5*(L-2*MT_RADIUS)*sin(pi/2-t(2)+b)];
            for j = 1:MTn
               for k = 1:MTn
                   x0 = x0V(j); y0 = y0V(j);
                   x0p = x0pV(k); y0p = y0pV(k); 
                   S = solve('(x-x0)^2+(y-y0)^2=MT_RADIUS^2','(x-x0p)^2+(y-y0p)^2=MT_RADIUS^2','x','y');
                   x = eval(S.x);
                   y = eval(S.y);
                   if isreal(x)
                       overlap = 1;
                       break; 
                   end
               end
            end
            if overlap && CORRECT_OVERLAP
               disp('It is not possible to reach the specified LK without overlap (showing closest solution)... Calculating (wait)...'); 
               for j = idx+1:length(bb)
                    x0V = [(r+CW_RADIUS+MT_RADIUS)*cos(t(1)) (r+CW_RADIUS+MT_RADIUS)*cos(t(1))-(L-2*MT_RADIUS)*cos(pi/2-t(1)+bb(j)) (r+CW_RADIUS+MT_RADIUS)*cos(t(1))-0.5*(L-2*MT_RADIUS)*cos(pi/2-t(1)+bb(j))];
                    y0V = [(r+CW_RADIUS+MT_RADIUS)*sin(t(1)) (r+CW_RADIUS+MT_RADIUS)*sin(t(1))+(L-2*MT_RADIUS)*sin(pi/2-t(1)+bb(j)) (r+CW_RADIUS+MT_RADIUS)*sin(t(1))+0.5*(L-2*MT_RADIUS)*sin(pi/2-t(1)+bb(j))];
                    x0pV = [(r+CW_RADIUS+MT_RADIUS)*cos(t(2)) (r+CW_RADIUS+MT_RADIUS)*cos(t(2))-(L-2*MT_RADIUS)*cos(pi/2-t(2)+bb(j)) (r+CW_RADIUS+MT_RADIUS)*cos(t(2))-0.5*(L-2*MT_RADIUS)*cos(pi/2-t(2)+bb(j))];
                    y0pV = [(r+CW_RADIUS+MT_RADIUS)*sin(t(2)) (r+CW_RADIUS+MT_RADIUS)*sin(t(2))+(L-2*MT_RADIUS)*sin(pi/2-t(2)+bb(j)) (r+CW_RADIUS+MT_RADIUS)*sin(t(2))+0.5*(L-2*MT_RADIUS)*sin(pi/2-t(2)+bb(j))];
                    x = []; y = [];
                    for l = 1:MTn
                       for k = 1:MTn
                           x0 = x0V(l); y0 = y0V(l);
                           x0p = x0pV(k); y0p = y0pV(k); 
                           S = solve('(x-x0)^2+(y-y0)^2=MT_RADIUS^2','(x-x0p)^2+(y-y0p)^2=MT_RADIUS^2','x','y');
                           if isreal(eval(S.x))
                               x(end+1,:) = eval(S.x);
                               y(end+1,:) = eval(S.y);
                           end
                       end
                    end
                    if isempty(x)
                       overlap = 0;
                       idx = j; % j or j+1?
                       b = bb(idx);
                       b_deg = round2xdigit(b*180/pi);
                       LKm = round2xdigit(LKmV(idx));
                       break; 
                    end
               end
            end
            if overlap && CORRECT_OVERLAP
                disp('No better solution was found... Showing with overlap.'); 
            elseif overlap && ~CORRECT_OVERLAP
                disp('Showing with overlap!');
            end
        end
    end

    % plot figure
    if show_result
        j = 1;
        for i = 1:SYM
            line([0 (r+CW_RADIUS)*cos(t(i))],[0 (r+CW_RADIUS)*sin(t(i))],'Color',[0 0 0]);
            if MTn > 2 % triplet
                line([(r+CW_RADIUS+MT_RADIUS)*cos(t(i))-(L-2*MT_RADIUS)*cos(pi/2-t(i)+b) - MT_RADIUS*cos(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)) + MT_RADIUS*cos(pi/2-t(i+1)+b-g)],[(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+(L-2*MT_RADIUS)*sin(pi/2-t(i)+b) + MT_RADIUS*sin(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i+1)) - MT_RADIUS*sin(pi/2-t(i+1)+b-g)],'Color',[0 0 0]);
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))-(L-2*MT_RADIUS)*cos(pi/2-t(i)+b)+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+(L-2*MT_RADIUS)*sin(pi/2-t(i)+b)+MT_RADIUS*(sin(tt)),col(j,:)); % C
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))-0.5*(L-2*MT_RADIUS)*cos(pi/2-t(i)+b)+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+0.5*(L-2*MT_RADIUS)*sin(pi/2-t(i)+b)+MT_RADIUS*(sin(tt)),col(j,:)); % B
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+MT_RADIUS*(sin(tt)),col(j,:)); % A
            elseif MTn > 1 % doublet
                line([(r+CW_RADIUS+MT_RADIUS)*cos(t(i))-(L-2*MT_RADIUS)*cos(pi/2-t(i)+b) - MT_RADIUS*cos(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)) + MT_RADIUS*cos(pi/2-t(i+1)+b-g)],[(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+(L-2*MT_RADIUS)*sin(pi/2-t(i)+b) + MT_RADIUS*sin(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i+1)) - MT_RADIUS*sin(pi/2-t(i+1)+b-g)],'Color',[0 0 0]);
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))-(L-2*MT_RADIUS)*cos(pi/2-t(i)+b)+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+(L-2*MT_RADIUS)*sin(pi/2-t(i)+b)+MT_RADIUS*(sin(tt)),col(j,:)); % B
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+MT_RADIUS*(sin(tt)),col(j,:)); % A
            elseif MTn == 1 % singlet
                line([(r+CW_RADIUS+MT_RADIUS)*cos(t(i)) - MT_RADIUS*cos(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*cos(t(i+1)) + MT_RADIUS*cos(pi/2-t(i+1)+b-g)],[(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+ MT_RADIUS*sin(pi/2-t(i)+b-g), (r+CW_RADIUS+MT_RADIUS)*sin(t(i+1)) - MT_RADIUS*sin(pi/2-t(i+1)+b-g)],'Color',[0 0 0]);
                patch((r+CW_RADIUS+MT_RADIUS)*cos(t(i))+MT_RADIUS*cos(tt),(r+CW_RADIUS+MT_RADIUS)*sin(t(i))+MT_RADIUS*(sin(tt)),col(j,:)); % A
            end
            j = j + 1; % MT color index
            if j > sqrt(SYM)
                j = 1;
            end
        end
        patch(CW_RADIUS*sin(tt),CW_RADIUS*cos(tt),[1 1 1]);
        
        %hold on
        %plot(x0V(2),y0V(2),'gx');
        %plot(x0pV(1),y0pV(1),'cx');
        
        xlim([-1.1*(r+CW_RADIUS+6*MT_RADIUS) 1.1*(r+CW_RADIUS+6*MT_RADIUS)]);
        ylim([-1.1*(r+CW_RADIUS+6*MT_RADIUS) 1.1*(r+CW_RADIUS+6*MT_RADIUS)]);
        axis equal; box on;
        %hold on; plot((r+CW_RADIUS)*cos(tt),(r+CW_RADIUS)*sin(tt),'c-');
        title(['SYM: ',num2str(SYM),'   LK: ',num2str(LK),'   CW: ',num2str(CW_RADIUS),'   SAS-6: ',num2str(r),'  MTl: ',num2str(L),'   MTr: ',num2str(MT_RADIUS)]);
        saveas(gcf,['SYM',num2str(SYM),'_LK',num2str(LK),'_CW',num2str(CW_RADIUS),'_SAS',num2str(r),'_MTl',num2str(L),'_MTr',num2str(MT_RADIUS),'.pdf']);
        saveas(gcf,['SYM',num2str(SYM),'_LK',num2str(LK),'_CW',num2str(CW_RADIUS),'_SAS',num2str(r),'_MTl',num2str(L),'_MTr',num2str(MT_RADIUS),'.fig']);

        % figure with LK(angle)
        figure; plot(bb*180/pi,LKmV);
        min(LKmV);
        xlabel('Angle (degrees)'); ylabel('LK');
    end
    
    % display the result
    disp(['Obtained LK: ',num2str(LKm),' with an angle beta = ',num2str(b_deg+90),'°.']);
    
end
    