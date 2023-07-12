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

%% 旧版本代码

% %% 储存点相关信息  point
% [c,~,ic] = unique([-fv.vertices(:,3),fv.vertices],'rows'); % c为不重复坐标点 ic为旧点在新点中所处的位置
% 
% point.tag = zeros(length(c),1);
% for i = 1:1:length(c)
% point.tag(i) = i;
% end
% 
% point.XYZ = c(:,2:4);
% 
% point.face_tag = zeros(length(c),100); % 创建一个足够大的矩阵，用于存放面片位置
% point.face_position_tag = zeros(length(c),100); % 创建一个足够大的矩阵，用于存放点在面上所处的位置
% for i = 1:1:length(ic)
%     j = sum(point.face_tag( ic(i),:) ~= 0) + 1; %非0元素的个数
%     point.face_tag( ic(i),j ) = ceil(i/3); %节点所包含的面所在的位置
%     
%     compare_a = fv.vertices( fv.faces( ceil(i/3) ,1) ,3); %原始面所对应的指定轴坐标
%     compare_b = fv.vertices( fv.faces( ceil(i/3) ,2) ,3);
%     compare_c = fv.vertices( fv.faces( ceil(i/3) ,3) ,3);
%     
%    compare_abc = point.XYZ( ic(i) ,3) ; % 排序后的面所对应的指定轴坐标
%    
%    if compare_abc == max([compare_a,compare_b,compare_c]) %节点在面上所处的位置，最大为1，最小为3，中间为2
%        point.face_position_tag( ic(i),j ) = 1;
%    else
%        if compare_abc == min([compare_a,compare_b,compare_c])
%            point.face_position_tag( ic(i),j ) = 3;
%        else
%            point.face_position_tag( ic(i),j ) = 2;
%        end
%    end  
% end
% 
% assign_face_line = all(point.face_tag==0,1) ;%指定行所在位置
% point.face_tag(:,assign_face_line)= []; %删除空白列
% point.face_position_tag(:,assign_face_line)= []; %删除空白列
% 
% clear i j c ic compare_abc compare_a compare_b compare_c assign_face_line
% 
% % line(point.XYZ(:,1),point.XYZ(:,2),point.XYZ(:,3));
% 
% %% 创建包围盒 bounding_volume
% % 包围盒边界
% [extremum.x_max,extremum.x_min,extremum.y_max,extremum.y_min,extremum.z_max,extremum.z_min] =...
%     deal( max(point.XYZ(:,1)),min(point.XYZ(:,1)),max(point.XYZ(:,2)),min(point.XYZ(:,2)),max(point.XYZ(:,3)),min(point.XYZ(:,3)) );
% 
% length_d.x = 1;length_d.y = 50;length_d.z = 20;
% bounding_volume.tag = zeros( length_d.y*length_d.z,3 ); % X轴方向不进行包容盒划分,Y轴方向划分为50个，Z轴方向划分为20个
% sum_ij = 0 ; bounding_volume.tag(:,1) = 1;
% for i = 1:1:length_d.z
%     for j = 1:1:length_d.y
%         sum_ij = sum_ij +1;
%         bounding_volume.tag(sum_ij,2) = j;
%         bounding_volume.tag(sum_ij,3) = i;
%     end
% end
% 
% bounding_volume.value = zeros( length_d.y*length_d.z,6 );
% 
% sum_ij = 0;
% for i = 1:1:length_d.z
%     for j = 1:1:length_d.y
%         sum_ij = sum_ij +1;
%         bounding_volume.value(sum_ij,1) = extremum.x_min;
%         bounding_volume.value(sum_ij,2) = extremum.x_max;
%         bounding_volume.value(sum_ij,3) = extremum.y_min + (j-1)*(extremum.y_max-extremum.y_min)/length_d.y;
%         bounding_volume.value(sum_ij,4) = extremum.y_min + j * (extremum.y_max-extremum.y_min)/length_d.y;
%         bounding_volume.value(sum_ij,5) = extremum.z_min + (i-1)*(extremum.z_max-extremum.z_min)/length_d.z;
%         bounding_volume.value(sum_ij,6) = extremum.z_min + i * (extremum.z_max-extremum.z_min)/length_d.z;
%     end
% end
% 
% bounding_volume.point_tag = zeros( length_d.y*length_d.z,100 ); % 创建一个足够大的矩阵，用于存放点的位置
% length_d.dx = linspace(extremum.x_min,extremum.x_max,length_d.x+1);
% length_d.dy = linspace(extremum.y_min,extremum.y_max,length_d.y+1);
% length_d.dz = linspace(extremum.z_min,extremum.z_max,length_d.z+1); %检索区间
% 
% point.length_d_tag(:,1) = discretize(point.XYZ(:,1),length_d.dx); %所处区间位置 左闭右开区间，最后一个区间为左闭右闭区间
% point.length_d_tag(:,2) = discretize(point.XYZ(:,2),length_d.dy); %所处区间位置
% point.length_d_tag(:,3) = discretize(point.XYZ(:,3),length_d.dz); %所处区间位置
% 
% 
% for i = 1:1:length(point.tag)
%     i_tag = (point.length_d_tag(i,3)-1)*50 + point.length_d_tag(i,2); %所处位置
%     sum_i = sum(bounding_volume.point_tag( i_tag,:) ~= 0) + 1; %非0元素的个数
%     bounding_volume.point_tag( i_tag,sum_i ) = point.tag(i); %所包含的节点所在的位置
% end
% 
% bounding_volume.point_tag(:,all(bounding_volume.point_tag==0,1))= []; %删除空白列
% assign_line = all(bounding_volume.point_tag==0,2) ;%指定行所在位置
% bounding_volume.point_tag(assign_line,:)= []; %删除空白行
% bounding_volume.tag(assign_line,:)= []; %删除空白行
% bounding_volume.value(assign_line,:)= []; %删除空白行
% 
% clear  i j sum_i sum_ij i_tag;
% 
% %% 轨迹规划——生成交点  trajectory planning
% 
% segment_num = 40; %沿所需规划的方向划分的段数
% traj_plan.tag = zeros( segment_num+1 ,2);
% traj_plan.volume = zeros( segment_num+1 ,length_d.z);
% traj_plan.face_tag = zeros( segment_num+1 ,5000);
% traj_plan.face_sort = zeros( segment_num+1 ,5000);
% [traj_plan.pos_X ,traj_plan.pos_Y ,traj_plan.pos_Z ] = deal( zeros( segment_num+1 ,2000) );
% 
% for i = 1:1:segment_num+1 %该循环是为了移动截面的位置
%     traj_plan.tag(i,1) = i; %标签，从最上方往最下方走
%     traj_plan.tag(i,2) = extremum.z_max - (i-1)*(extremum.z_max-extremum.z_min)/segment_num; %截面标签位置
%     if i ==  1
%         tag_z = 20;
%     else
%         tag_z = ceil( (segment_num+2-i)/2 );
%     end
%     tag_volume = flip( find(bounding_volume.tag(:,3) == tag_z ) ); %对于一个行向量，使用flip(A)默认可以将列表左右翻转
%     traj_plan.volume(i , 1:length(tag_volume) ) = tag_volume'; %将对应的盒子放入对应的位置
%     
%     for i_volume = 1:1:sum(traj_plan.volume( i,:) ~= 0) %包含多少个盒子，循环就有多少次
%         for i_point = 1:1:sum(bounding_volume.point_tag ( traj_plan.volume(i,i_volume) ,:) ~= 0) %盒子里包含多少个点，循环就有多少次
%             i_point_line = bounding_volume.point_tag ( traj_plan.volume(i,i_volume) , i_point);  %点在盒子中的值，对应的也是点所被包含面中的位置
%             for i_face = 1:1:sum( point.face_position_tag( i_point_line ,:) ~= 0) %每个点被包含于多少个面中，循环就有多少次
%                 
%                 sum_i = sum(traj_plan.face_tag( i,:) ~= 0) + 1; %非0元素的个数
%                 traj_plan.face_tag(i,sum_i) = point.face_tag( i_point_line ,i_face);
%       
%             end
%             
%         end
% 
%     end
%     % 三次循环已经读取了所有的面片
%     
%     AAA = unique( traj_plan.face_tag(i,:)); % traj_plan.face_tag 矩阵的列维度要大，不然第三行会报错——索引超出矩阵维度
%     length_face_tag = sum( AAA  ~= 0);
%     traj_plan.face_sort(i,1:length_face_tag) = AAA(2:end) ; %第一个值为零，所以从第二个开始
%     
%     % 最重要的循环，主要是求轨迹线的交点
%     for i_traj = 1:1:length_face_tag % 截面与每个面片相交
%         face_pos = traj_plan.face_sort(i,i_traj); %排序后面所在的位置
%  
%         compare_a = fv.vertices( fv.faces( face_pos ,1) ,3); %原始面所对应的指定轴坐标
%         compare_b = fv.vertices( fv.faces( face_pos ,2) ,3);
%         compare_c = fv.vertices( fv.faces( face_pos ,3) ,3);
%         if traj_plan.tag(i,2) > max([compare_a,compare_b,compare_c]) %面片中要有点的坐标大于指定轴坐标，才有交点
%         else
%             Z = traj_plan.tag(i,2);
%             [X1,Y1,Z1] = deal( fv.vertices( fv.faces( face_pos ,1) ,1) , fv.vertices( fv.faces( face_pos ,1) ,2) ,fv.vertices( fv.faces( face_pos ,1) ,3) );
%             [X2,Y2,Z2] = deal( fv.vertices( fv.faces( face_pos ,2) ,1) , fv.vertices( fv.faces( face_pos ,2) ,2) ,fv.vertices( fv.faces( face_pos ,2) ,3) );
%             [X3,Y3,Z3] = deal( fv.vertices( fv.faces( face_pos ,3) ,1) , fv.vertices( fv.faces( face_pos ,3) ,2) ,fv.vertices( fv.faces( face_pos ,3) ,3) );
%   
%             for i_side =1:1:3 %一个面对应三条边
%                 switch i_side
%                     
%                     case 1
%                         X = (X1-X2)*(Z-Z2)/(Z1-Z2) + X2;
%                         Y = (Y1-Y2)*(Z-Z2)/(Z1-Z2) + Y2;
%                         if Z > max([Z1,Z2]) || Z < min([Z1,Z2])
%                             % 无交点
%                         else
%                             sum_pos_i = sum(traj_plan.pos_X( i,:) ~= 0) + 1;
%                             traj_plan.pos_X ( i,sum_pos_i) = X;
%                             traj_plan.pos_Y ( i,sum_pos_i) = Y;
%                             traj_plan.pos_Z ( i,sum_pos_i) = Z;
%                             
%                         end
%                         
%                     case 2
%                         X = (X1-X3)*(Z-Z3)/(Z1-Z3) + X3;
%                         Y = (Y1-Y3)*(Z-Z3)/(Z1-Z3) + Y3;
%                         if Z > max([Z1,Z3]) || Z < min([Z1,Z3])
%                             % 无交点
%                         else
%                             sum_pos_i = sum(traj_plan.pos_X( i,:) ~= 0) + 1;
%                             traj_plan.pos_X ( i,sum_pos_i) = X;
%                             traj_plan.pos_Y ( i,sum_pos_i) = Y;
%                             traj_plan.pos_Z ( i,sum_pos_i) = Z;
%                             
%                         end
%                     case 3
%                         X = (X3-X2)*(Z-Z2)/(Z3-Z2) + X2;
%                         Y = (Y3-Y2)*(Z-Z2)/(Z3-Z2) + Y2;
%                         if Z > max([Z3,Z2]) || Z < min([Z3,Z2])
%                             % 无交点
%                         else
%                             sum_pos_i = sum(traj_plan.pos_X( i,:) ~= 0) + 1;
%                             traj_plan.pos_X ( i,sum_pos_i) = X;
%                             traj_plan.pos_Y ( i,sum_pos_i) = Y;
%                             traj_plan.pos_Z ( i,sum_pos_i) = Z;
%                             
%                         end
%                 end
%                 
%             end
%         end
%     end
% end
%     
% clear  X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 X Y Z
% clear i tag_z i_face i_point i_point_line i_volume sum_i AAA length_face_tag   
% 
% %% 轨迹规划——排序
% length_point = 1;
% traj_plan_XYZ = zeros(20000,4);
% for i = 1:1:segment_num+1
% 
%     XYZ(:,1) = (traj_plan.pos_Y(i,:))' ;
%     XYZ(:,2) = (traj_plan.pos_X(i,:))' ;
%     XYZ(:,3) = (traj_plan.pos_Y(i,:))' ;
%     XYZ(:,4) = (traj_plan.pos_Z(i,:))' ;
%     AAA = unique ( XYZ,'rows' );
%     
%     length_xyz = size( AAA,1)-1;
%     traj_plan_XYZ(length_point:length_point+length_xyz-1 ,1) = i;
%     traj_plan_XYZ(length_point:length_point+length_xyz-1 ,2:4) = AAA(2:end,2:4);
%     length_point = length_point+length_xyz;
%     
%     
% end
% 
% traj_plan_XYZ(all(traj_plan_XYZ==0,2) ,:)= []; %删除空白行
% 
% line(traj_plan_XYZ(:,2),traj_plan_XYZ(:,3),traj_plan_XYZ(:,4));
toc