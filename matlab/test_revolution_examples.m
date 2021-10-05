ex = struct('a',-20,'b',[10;3],'alpha0',1.5701,'total_length',2.628,'gamma0',[1.3844;0.9776],'num_spokes',20);
ex(2) = struct('a',-1,'b',[0.2;1],'alpha0',pi*0.4,'total_length',8,'gamma0',[1;-0.5],'num_spokes',9);
ex(3) = struct('a',-1*0.4,'b',[-0.5;1]*0.4,'alpha0',pi*0.9,'total_length',6.2,'gamma0',[4;-0.8],'num_spokes',5);
ex(4) = struct('a',-1*1.5,'b',[0;1]*1.5,'alpha0',pi*0.7,'total_length',6.9,'gamma0',[3;-0.23],'num_spokes',15);
ex(5) = struct('a',-1,'b',[1;-0.3],'alpha0',pi*0.7,'total_length',9,'gamma0',[3;0],'num_spokes',6);
ex(6) = struct('a',-12*0.1,'b',[1;-0.3]*0.1,'alpha0',0,'total_length',0.96,'gamma0',[0.2;0],'num_spokes',8);
% ex(7) = struct('a',-1,'b',[1;0],'alpha0',pi*1.15,'total_length',5.8,'gamma0',[2.5;0],'num_spokes',4);

extra_plot = false;

gray = [43 48 43]/255;
colors = [60 135 60;
    214 124 83;
    138 183 138;
    94 46 117;
    158 130 172;
    204 91 40;
    30 68 30]/255;

exi_range = 1:length(ex);
% exi_range = length(ex);
for exi=exi_range
    a = ex(exi).a;
    v = [1;0];
    b = ex(exi).b;
    alpha0 = ex(exi).alpha0;
    % solve for (k*kappa - 1)*<gamma,e_1> = a

    total_length = ex(exi).total_length;
    ds = 2.e-3;
    num_steps = ceil(total_length / ds);

    alpha = zeros(1,num_steps+1);
    gamma = zeros(2,num_steps+1);
    K = zeros(1,num_steps);
    curv = zeros(1,num_steps);

    gamma(:,1) = ex(exi).gamma0;
    alpha(:,1) = alpha0;

    for i=1:num_steps
        curv(i) = (a + dot(gamma(:,i), b)) / dot(gamma(:,i), v);
        alpha(i+1) = alpha(i) + ds * curv(i);
        gamma(:,i+1) = gamma(:,i) + ds * [cos(alpha(i+1)); sin(alpha(i+1))];
        K(:,i) = dot(gamma(:,i), v);
    end
    s = ds*(0:num_steps-1);

    num_spokes = ex(exi).num_spokes;
    half_angle = pi/num_spokes;
    hw0 = tan(half_angle) * gamma(1,:);

    
    p = [[gamma(1,:); -hw0; gamma(2,:)],[gamma(1,:); hw0; gamma(2,:)]];
    num_samples = size(gamma,2);
    F = GeometryGenerator.generateStripFaceList(num_samples);
    figure('Color','white');
    for si=1:num_spokes
        angle = 2*pi*(si-1)/num_spokes;
        R = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
        Rp = R*p;
        patch('Vertices',Rp','Faces',F,'EdgeColor','none','FaceColor',colors(exi,:),...
            'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',1);
        patch('Vertices',Rp(:,[1 num_samples+1 num_samples 2*num_samples])', 'Faces',[1 2;3 4],...
            'EdgeColor',gray,'LineWidth',3);
    end
    light('Color',[0.5 0.5 0.5],'Position',[0.5 0 1]);
    light('Color',[0.5 0.5 0.5],'Position',[0 0.5 1]);
    light('Color',[0.2 0.2 0.2],'Position',[0.5 0 -1]);
    light('Color',[0.2 0.2 0.2],'Position',[0 0.5 -1]);
    view(75,16);

    axis tight equal off;
    title(sprintf('Fig. 22 (%u/6)', exi));
    
    if extra_plot && exi==1
        figure('Color','white');
        
        hw0 = tan(pi/8) * gamma(1,:);
        p = [[gamma(1,:); -hw0; gamma(2,:)],[gamma(1,:); hw0; gamma(2,:)]];
        
        axis_points = [[0;0;p(3,1)], [0;0;p(3,end)]];
        
        patch('Vertices',p','Faces',F,'EdgeColor','none','FaceColor',colors(exi,:),...
            'SpecularStrength',0.3,'AmbientStrength',0.5,'DiffuseStrength',1);
        patch('Vertices',[axis_points, p(:,[1, num_samples+1, num_samples, 2*num_samples])]', ...
            'Faces',[1 2;1 3;1 4;2 5;2 6]);
        
        [~,i_cut] = max(gamma(1,:));
        patch('Vertices',[[0;0;p(3,i_cut)], p(:,[i_cut, num_samples+i_cut])]', ...
            'Faces',[1 2;1 3;2 3]);
        
        light('Color',[0.5 0.5 0.5],'Position',[0.5 0 1]);
        light('Color',[0.5 0.5 0.5],'Position',[0 0.5 1]);
        light('Color',[0.2 0.2 0.2],'Position',[0.5 0 -1]);
        light('Color',[0.2 0.2 0.2],'Position',[0 0.5 -1]);
        view(44,16);
        
        axis tight equal off;
        
        return;
    end
    
%     figure;
%     ax=axes();
%     plot(gamma(1,:), gamma(2,:));
%     axis tight equal;
%     xl = ax.XLim;
%     yl = ax.YLim;
%     
%     line_tangent = [-b(2);b(1)];
%     point_on_line =-a*b/dot(b,b);
%     patch('Vertices',[point_on_line-line_tangent*100, point_on_line+line_tangent*100]', 'Faces',[1 2]);
%     ax.XLim = [0 xl(2)];
%     ax.YLim = yl;
end


