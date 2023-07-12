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
fv = stlread('����׹���.STL');


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

%% ����������Ϣ  point


[ point.xyz,~,fv.sort_tag ] = unique(fv.vertices,'rows'); % point.xyz Ϊ���ظ������µ� fv.sort_tag Ϊ�ɵ����µ���������λ��

% point.tag = zeros(length( point.xyz ),1);
% for i = 1:1:length( point.xyz )
% point.tag(i) = i;
% end

%% ����ÿ����
pos_xyz = zeros( 20000,3) ;

fv.combination = zeros( length(fv.sort_tag) ,2);  % ����ϵĿ�����
[extremum.x_max,extremum.x_min,extremum.y_max,extremum.y_min,extremum.z_max,extremum.z_min] =...
    deal( max(point.xyz(:,1)),min(point.xyz(:,1)),max(point.xyz(:,2)),min(point.xyz(:,2)),max(point.xyz(:,3)),min(point.xyz(:,3)) ); % ��ֵextremum

length_d.x = 1;length_d.y = 80;length_d.z = 40;
length_d.subsection_z = linspace(extremum.z_min,extremum.z_max,length_d.z+1); %Z����ķֶ���� subsection �ֶ�

for i = 1:1:size(fv.faces,1)
    
    [X1,Y1,Z1] = deal( fv.vertices( fv.faces( i ,1) ,1) , fv.vertices( fv.faces( i ,1) ,2) ,fv.vertices( fv.faces( i ,1) ,3) );
    [X2,Y2,Z2] = deal( fv.vertices( fv.faces( i ,2) ,1) , fv.vertices( fv.faces( i ,2) ,2) ,fv.vertices( fv.faces( i ,2) ,3) );
    [X3,Y3,Z3] = deal( fv.vertices( fv.faces( i ,3) ,1) , fv.vertices( fv.faces( i ,3) ,2) ,fv.vertices( fv.faces( i ,3) ,3) );
       
    for i_side =1:1:3 %һ�����Ӧ������
        switch i_side
            
            case 1  % i_side = 1
                j = sum(fv.combination( :,1) ~= 0) + 1; %��0Ԫ�صĸ���
                a = min( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,2) ) );
                b = max( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,2) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %ȥ���ظ�����
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z1 Z2],length_d.subsection_z); % Section����  discretize ����ҿ�����
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z1,Z2]) || Z < min([Z1,Z2])
                            % �޽���
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
                j = sum(fv.combination( :,1) ~= 0) + 1; %��0Ԫ�صĸ���
                a = min( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,3) ) );
                b = max( fv.sort_tag( fv.faces(i,1) ) ,fv.sort_tag( fv.faces(i,3) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %ȥ���ظ�����
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z1 Z3],length_d.subsection_z); % Section����  discretize ����ҿ�����
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z1,Z3]) || Z < min([Z1,Z3])
                            % �޽���
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
                j = sum(fv.combination( :,1) ~= 0) + 1; %��0Ԫ�صĸ���
                a = min( fv.sort_tag( fv.faces(i,2) ) ,fv.sort_tag( fv.faces(i,3) ) );
                b = max( fv.sort_tag( fv.faces(i,2) ) ,fv.sort_tag( fv.faces(i,3) ) );
                
                if ismember(fv.combination, [a b],'rows') == 0 %ȥ���ظ�����
                    fv.combination(j,1) = a; fv.combination(j,2) = b;

                    section_z = discretize([Z2 Z3],length_d.subsection_z); % Section����  discretize ����ҿ�����
                    section_z_a = min( section_z(1),section_z(2) ); section_z_b = max( section_z(1),section_z(2) );
                    
                    if section_z_b == length_d.z
                        section_z_b = section_z_b + 1;
                    end
                    
                    for tag_subsection = section_z_a:1:section_z_b
                        Z = length_d.subsection_z(tag_subsection);
                        
                        if Z > max([Z3,Z2]) || Z < min([Z3,Z2])
                            % �޽���
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

%% �켣�滮


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

%%  ���Ƭ
for i = 2:1:length(traj_plan.sort_xyz)
    traj_plan.sort_xyz(i,5) = traj_plan.sort_xyz(i,3) - traj_plan.sort_xyz(i-1,3);
end

%% �ɰ汾����

% %% ����������Ϣ  point
% [c,~,ic] = unique([-fv.vertices(:,3),fv.vertices],'rows'); % cΪ���ظ������ icΪ�ɵ����µ���������λ��
% 
% point.tag = zeros(length(c),1);
% for i = 1:1:length(c)
% point.tag(i) = i;
% end
% 
% point.XYZ = c(:,2:4);
% 
% point.face_tag = zeros(length(c),100); % ����һ���㹻��ľ������ڴ����Ƭλ��
% point.face_position_tag = zeros(length(c),100); % ����һ���㹻��ľ������ڴ�ŵ�������������λ��
% for i = 1:1:length(ic)
%     j = sum(point.face_tag( ic(i),:) ~= 0) + 1; %��0Ԫ�صĸ���
%     point.face_tag( ic(i),j ) = ceil(i/3); %�ڵ��������������ڵ�λ��
%     
%     compare_a = fv.vertices( fv.faces( ceil(i/3) ,1) ,3); %ԭʼ������Ӧ��ָ��������
%     compare_b = fv.vertices( fv.faces( ceil(i/3) ,2) ,3);
%     compare_c = fv.vertices( fv.faces( ceil(i/3) ,3) ,3);
%     
%    compare_abc = point.XYZ( ic(i) ,3) ; % ������������Ӧ��ָ��������
%    
%    if compare_abc == max([compare_a,compare_b,compare_c]) %�ڵ�������������λ�ã����Ϊ1����СΪ3���м�Ϊ2
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
% assign_face_line = all(point.face_tag==0,1) ;%ָ��������λ��
% point.face_tag(:,assign_face_line)= []; %ɾ���հ���
% point.face_position_tag(:,assign_face_line)= []; %ɾ���հ���
% 
% clear i j c ic compare_abc compare_a compare_b compare_c assign_face_line
% 
% % line(point.XYZ(:,1),point.XYZ(:,2),point.XYZ(:,3));
% 
% %% ������Χ�� bounding_volume
% % ��Χ�б߽�
% [extremum.x_max,extremum.x_min,extremum.y_max,extremum.y_min,extremum.z_max,extremum.z_min] =...
%     deal( max(point.XYZ(:,1)),min(point.XYZ(:,1)),max(point.XYZ(:,2)),min(point.XYZ(:,2)),max(point.XYZ(:,3)),min(point.XYZ(:,3)) );
% 
% length_d.x = 1;length_d.y = 50;length_d.z = 20;
% bounding_volume.tag = zeros( length_d.y*length_d.z,3 ); % X�᷽�򲻽��а��ݺл���,Y�᷽�򻮷�Ϊ50����Z�᷽�򻮷�Ϊ20��
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
% bounding_volume.point_tag = zeros( length_d.y*length_d.z,100 ); % ����һ���㹻��ľ������ڴ�ŵ��λ��
% length_d.dx = linspace(extremum.x_min,extremum.x_max,length_d.x+1);
% length_d.dy = linspace(extremum.y_min,extremum.y_max,length_d.y+1);
% length_d.dz = linspace(extremum.z_min,extremum.z_max,length_d.z+1); %��������
% 
% point.length_d_tag(:,1) = discretize(point.XYZ(:,1),length_d.dx); %��������λ�� ����ҿ����䣬���һ������Ϊ����ұ�����
% point.length_d_tag(:,2) = discretize(point.XYZ(:,2),length_d.dy); %��������λ��
% point.length_d_tag(:,3) = discretize(point.XYZ(:,3),length_d.dz); %��������λ��
% 
% 
% for i = 1:1:length(point.tag)
%     i_tag = (point.length_d_tag(i,3)-1)*50 + point.length_d_tag(i,2); %����λ��
%     sum_i = sum(bounding_volume.point_tag( i_tag,:) ~= 0) + 1; %��0Ԫ�صĸ���
%     bounding_volume.point_tag( i_tag,sum_i ) = point.tag(i); %�������Ľڵ����ڵ�λ��
% end
% 
% bounding_volume.point_tag(:,all(bounding_volume.point_tag==0,1))= []; %ɾ���հ���
% assign_line = all(bounding_volume.point_tag==0,2) ;%ָ��������λ��
% bounding_volume.point_tag(assign_line,:)= []; %ɾ���հ���
% bounding_volume.tag(assign_line,:)= []; %ɾ���հ���
% bounding_volume.value(assign_line,:)= []; %ɾ���հ���
% 
% clear  i j sum_i sum_ij i_tag;
% 
% %% �켣�滮�������ɽ���  trajectory planning
% 
% segment_num = 40; %������滮�ķ��򻮷ֵĶ���
% traj_plan.tag = zeros( segment_num+1 ,2);
% traj_plan.volume = zeros( segment_num+1 ,length_d.z);
% traj_plan.face_tag = zeros( segment_num+1 ,5000);
% traj_plan.face_sort = zeros( segment_num+1 ,5000);
% [traj_plan.pos_X ,traj_plan.pos_Y ,traj_plan.pos_Z ] = deal( zeros( segment_num+1 ,2000) );
% 
% for i = 1:1:segment_num+1 %��ѭ����Ϊ���ƶ������λ��
%     traj_plan.tag(i,1) = i; %��ǩ�������Ϸ������·���
%     traj_plan.tag(i,2) = extremum.z_max - (i-1)*(extremum.z_max-extremum.z_min)/segment_num; %�����ǩλ��
%     if i ==  1
%         tag_z = 20;
%     else
%         tag_z = ceil( (segment_num+2-i)/2 );
%     end
%     tag_volume = flip( find(bounding_volume.tag(:,3) == tag_z ) ); %����һ����������ʹ��flip(A)Ĭ�Ͽ��Խ��б����ҷ�ת
%     traj_plan.volume(i , 1:length(tag_volume) ) = tag_volume'; %����Ӧ�ĺ��ӷ����Ӧ��λ��
%     
%     for i_volume = 1:1:sum(traj_plan.volume( i,:) ~= 0) %�������ٸ����ӣ�ѭ�����ж��ٴ�
%         for i_point = 1:1:sum(bounding_volume.point_tag ( traj_plan.volume(i,i_volume) ,:) ~= 0) %������������ٸ��㣬ѭ�����ж��ٴ�
%             i_point_line = bounding_volume.point_tag ( traj_plan.volume(i,i_volume) , i_point);  %���ں����е�ֵ����Ӧ��Ҳ�ǵ������������е�λ��
%             for i_face = 1:1:sum( point.face_position_tag( i_point_line ,:) ~= 0) %ÿ���㱻�����ڶ��ٸ����У�ѭ�����ж��ٴ�
%                 
%                 sum_i = sum(traj_plan.face_tag( i,:) ~= 0) + 1; %��0Ԫ�صĸ���
%                 traj_plan.face_tag(i,sum_i) = point.face_tag( i_point_line ,i_face);
%       
%             end
%             
%         end
% 
%     end
%     % ����ѭ���Ѿ���ȡ�����е���Ƭ
%     
%     AAA = unique( traj_plan.face_tag(i,:)); % traj_plan.face_tag �������ά��Ҫ�󣬲�Ȼ�����лᱨ����������������ά��
%     length_face_tag = sum( AAA  ~= 0);
%     traj_plan.face_sort(i,1:length_face_tag) = AAA(2:end) ; %��һ��ֵΪ�㣬���Դӵڶ�����ʼ
%     
%     % ����Ҫ��ѭ������Ҫ����켣�ߵĽ���
%     for i_traj = 1:1:length_face_tag % ������ÿ����Ƭ�ཻ
%         face_pos = traj_plan.face_sort(i,i_traj); %����������ڵ�λ��
%  
%         compare_a = fv.vertices( fv.faces( face_pos ,1) ,3); %ԭʼ������Ӧ��ָ��������
%         compare_b = fv.vertices( fv.faces( face_pos ,2) ,3);
%         compare_c = fv.vertices( fv.faces( face_pos ,3) ,3);
%         if traj_plan.tag(i,2) > max([compare_a,compare_b,compare_c]) %��Ƭ��Ҫ�е���������ָ�������꣬���н���
%         else
%             Z = traj_plan.tag(i,2);
%             [X1,Y1,Z1] = deal( fv.vertices( fv.faces( face_pos ,1) ,1) , fv.vertices( fv.faces( face_pos ,1) ,2) ,fv.vertices( fv.faces( face_pos ,1) ,3) );
%             [X2,Y2,Z2] = deal( fv.vertices( fv.faces( face_pos ,2) ,1) , fv.vertices( fv.faces( face_pos ,2) ,2) ,fv.vertices( fv.faces( face_pos ,2) ,3) );
%             [X3,Y3,Z3] = deal( fv.vertices( fv.faces( face_pos ,3) ,1) , fv.vertices( fv.faces( face_pos ,3) ,2) ,fv.vertices( fv.faces( face_pos ,3) ,3) );
%   
%             for i_side =1:1:3 %һ�����Ӧ������
%                 switch i_side
%                     
%                     case 1
%                         X = (X1-X2)*(Z-Z2)/(Z1-Z2) + X2;
%                         Y = (Y1-Y2)*(Z-Z2)/(Z1-Z2) + Y2;
%                         if Z > max([Z1,Z2]) || Z < min([Z1,Z2])
%                             % �޽���
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
%                             % �޽���
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
%                             % �޽���
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
% %% �켣�滮��������
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
% traj_plan_XYZ(all(traj_plan_XYZ==0,2) ,:)= []; %ɾ���հ���
% 
% line(traj_plan_XYZ(:,2),traj_plan_XYZ(:,3),traj_plan_XYZ(:,4));
toc