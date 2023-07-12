clear
tic
%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.

% Copyright 2011 The MathWorks, Inc.


%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv = stlread('轮毂抛光面.STL');


%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.
patch_fv.faces = fv.faces;
patch_fv.vertices = fv.vertices;

patch(patch_fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([70 20]);

clear patch_fv;

%% 储存点相关信息  point


[ point.xyz,~,fv.sort_tag ] = unique(fv.vertices,'rows'); % point.xyz 为不重复坐标新点 fv.sort_tag 为旧点在新点中所处的位置

% point.tag = zeros(length( point.xyz ),1);
% for i = 1:1:length( point.xyz )
% point.tag(i) = i;
% end

%% 遍历每个面
pos_xyz = zeros( 20000,3) ;

fv.combination = zeros( length(fv.sort_tag) ,2);  % 点组合的可能性
[extremum.x_max,extremum.x_min,extremum.y_max,extremum.y_min,extremum.z_max,extremum.z_min] =...
    deal( max(point.xyz(:,1)),min(point.xyz(:,1)),max(point.xyz(:,2)),min(point.xyz(:,2)),max(point.xyz(:,3)),min(point.xyz(:,3)) ); % 极值extremum

length_d.x = 1;length_d.y = 80;length_d.z = 40;
length_d.subsection_z = linspace(extremum.z_min,extremum.z_max,length_d.z+1); %Z方向的分段情况 subsection 分段

for i = 1:1:size(fv.faces,1)
    
    [X1,Y1,Z1] = deal( fv.vertices( fv.faces( i ,1) ,1) , fv.vertices( fv.faces( i ,1) ,2) ,fv.vertices( fv.faces( i ,1) ,3) );
    [X2,Y2,Z2] = deal( fv.vertices( fv.faces( i ,2) ,1) , fv.vertices( fv.faces( i ,2) ,2) ,fv.vertices( fv.faces( i ,2) ,3) );
    [X3,Y3,Z3] = deal( fv.vertices( fv.faces( i ,3) ,1) , fv.vertices( fv.faces( i ,3) ,2) ,fv.vertices( fv.faces( i ,3) ,3) );
       
    for i_side =1:1:3 %一个面对应三条边
        switch i_side
            
            case 1  % i_side = 1
                j = sum(fv.combination( :,1) ~= 0) + 1; %非0元素的个数
                a = min( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,2) ) );
                b = max( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,2) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %去除重复交线
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z1 Z2],length_d.subsection_z); % Section区间  discretize 左闭右开区间
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z1,Z2]) || Z < min([Z1,Z2])
                            % 无交点
                        else
                            
                            sum_pos_i = sum(pos_xyz( :,1) ~= 0) + 1;
                            X = (X1-X2)*(Z-Z2)/(Z1-Z2) + X2;
                            Y = (Y1-Y2)*(Z-Z2)/(Z1-Z2) + Y2;
                            
                            pos_xyz ( sum_pos_i,1) = X;
                            pos_xyz ( sum_pos_i,2) = Y;
                            pos_xyz ( sum_pos_i,3) = Z;
                            
                        end
                    end
                                     
                else
                end

           case 2  % i_side = 2
                j = sum(fv.combination( :,1) ~= 0) + 1; %非0元素的个数
                a = min( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,3) ) );
                b = max( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,3) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %去除重复交线
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z1 Z3],length_d.subsection_z); % Section区间  discretize 左闭右开区间
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z1,Z3]) || Z < min([Z1,Z3])
                            % 无交点
                        else
                            
                            sum_pos_i = sum(pos_xyz( :,1) ~= 0) + 1;
                            X = (X1-X3)*(Z-Z3)/(Z1-Z3) + X3;
                            Y = (Y1-Y3)*(Z-Z3)/(Z1-Z3) + Y3;
                            
                            pos_xyz ( sum_pos_i,1) = X;
                            pos_xyz ( sum_pos_i,2) = Y;
                            pos_xyz ( sum_pos_i,3) = Z;
                            
                        end
                    end
                                     
                else
                end


           case 3  % i_side = 3
                j = sum(fv.combination( :,1) ~= 0) + 1; %非0元素的个数
                a = min( fv.sort_tag( fv.faces(i,2) ) ,fv.sort_tag( fv.faces(i,3) ) );
                b = max( fv.sort_tag( fv.faces(i,2) ) ,fv.sort_tag( fv.faces(i,3) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %去除重复交线
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z2 Z3],length_d.subsection_z); % Section区间  discretize 左闭右开区间
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z3,Z2]) || Z < min([Z3,Z2])
                            % 无交点
                        else
                            
                            sum_pos_i = sum(pos_xyz( :,1) ~= 0) + 1;
                            X = (X3-X2)*(Z-Z2)/(Z3-Z2) + X2;
                            Y = (Y3-Y2)*(Z-Z2)/(Z3-Z2) + Y2;
                            
                            pos_xyz ( sum_pos_i,1) = X;
                            pos_xyz ( sum_pos_i,2) = Y;
                            pos_xyz ( sum_pos_i,3) = Z;
                            
                        end
                    end
                                     
                else
                end

        end
    end
end

%% 轨迹规划


traj_plan.pos_xyz = unique([-pos_xyz(:,3),pos_xyz],'rows');

traj_plan.sort_xyz = zeros(length(traj_plan.pos_xyz)-1,5);  % temporary
sum_i = 1; traj_plan.sort_xyz(1,1) = 1;
for i = 2:1:length(traj_plan.pos_xyz)-1
    if traj_plan.pos_xyz(i,1) == traj_plan.pos_xyz(i-1,1)       
    else
        sum_i = sum_i + 1;
    end
    traj_plan.sort_xyz(i,1) = sum_i;    
end

left_i = 0;right_i = 0;
for i = 1:1:41
    left_right = sum(traj_plan.sort_xyz(:,1)==i);
    left_i = right_i + 1;
    right_i = right_i + left_right;
    clear temp_xyz temp_xyz1
    temp_xyz = zeros( left_right,4);
    temp_xyz(1:left_right,1) = -traj_plan.pos_xyz(left_i:right_i,3);
    temp_xyz(1:left_right,2:4) = traj_plan.pos_xyz(left_i:right_i,2:4);
    temp_xyz1 = unique(temp_xyz,'rows');
    traj_plan.sort_xyz(left_i:right_i,2:4) = temp_xyz1(1:left_right,2:4);
    hold on;
    line( temp_xyz1(:,2),temp_xyz1(:,3),temp_xyz1(:,4) );

    
end
%%
% plot3(traj_plan.pos_xyz(:,2),traj_plan.pos_xyz(:,3),traj_plan.pos_xyz(:,4),'o')
plot3(traj_plan.sort_xyz(:,2),traj_plan.sort_xyz(:,3),traj_plan.sort_xyz(:,4),'-')
% [aaa,ccc,bbb] = unique(-pos_xyz(:,3),'rows');

%%  点分片
for i = 2:1:length(traj_plan.sort_xyz)
    traj_plan.sort_xyz(i,5) = traj_plan.sort_xyz(i,3) - traj_plan.sort_xyz(i-1,3);
end

toc
